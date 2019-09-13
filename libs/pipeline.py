#!/usr/env/bin python

# import off-the-shelf libraries
import os
import sys
import GPy, scipy
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.signal import savgol_filter
from scipy.stats import spearmanr

print '%s\t\t %s' % ('pandas',pd.__version__)
print '%s\t\t %s' % ('numpy', np.__version__)
print '%s\t\t %s' % ('scipy', scipy.__version__)
print '%s\t\t %s' % ('GPy', GPy.__version__)
print '%s\t\t %s' % ('seaborn', sns.__version__)
print '%s\t %s' % ('matplotlib', mpl.__version__)
print 

# set global parameters
sns.set_style('whitegrid')

#import in-house library
sys.path.append("/data/davidalb/users/fsm/biolog/phenotypic-characterization/")
from libs import classical,growth,plates

#sys.path.append("/data/davidalb/users/fsm/biolog/metadata")
#from flags import *

def readPlateReaderFolder(folderpath,save=False,save_dirname='../data_formatted',interval=600):
    
    df_dict = {};

    filepaths = plates.findPlateReaderFiles(folderpath)

    for filepath in sorted(filepaths):

        _, filebase, newfile = plates.breakDownFilePath(filepath,save_dirname);

        #df = plates.readPlateReaderData(filepath,save=save,save_dirname=save_dirname);
        df = plates.readPlateReaderData(filepath,save=save,save_dirname=save_dirname);
        #df = formatPlateData(df,interval=interval);
        #df = df.T.reset_index(drop=False); # now we can reset index to a random number 
        
        # Time A1 A2 ... H12
        # 0 0.1 0.102 ... 0.09 
        # ...
        # 6000 1.1 1.05 ... 1.0

        df_dict[filebase] = df;

        print df.head()

        if save:
            df.to_csv(newfile, sep='\t',header=True) # does not save header index name (i.e. Time)

    return df_dict

def readPlateReaderMapping(filepath,filter_plates=True):
    
    mapping_df = pd.read_csv(filepath,sep='\t',header=0,index_col=0)
    
    if filter_plates:
        to_remove = mapping_df[mapping_df.Ignore=='Yes'].index
        mapping_df = mapping_df.drop(labels=to_remove,axis=0)

    return mapping_df

def initializePlateMapping(data_dict,flagged_plates):
    
    # initialize list of plates based on available data
    plate_list_df = pd.DataFrame(index=data_dict.keys(),
                                 columns=['Isolate','PM','Rep']);

    for key in data_dict.keys():

        basename = key.split('_')[0];
        iso = basename;
        pmn = int(key.split('PM')[1][0]);
        rep = [int(key.split('-')[-1]) if '-' in key else 1][0];
        plate_list_df.loc[key] = [iso,pmn,rep]#,ribo];
    
    plate_list_df = plate_list_df.sort_values(['Isolate','PM','Rep'],ascending=True);
    plate_list_df = plate_list_df.dropna(axis=0)

    plate_list_df = plate_list_df.drop(flagged_plates)
    
    return plate_list_df

def expandPlateMapping(data_dict,mapping_df,flagged_plates):
    
    # initialize list of plates based on available data
    plate_list_df = pd.DataFrame(index=data_dict.keys(),
                                 columns=['Isolate','PM','Rep','Ribotype'])

    for key in data_dict.keys():

    	# if plate id is not in mapping data, ignore
        basename = key.split('_')[0];
        if basename not in mapping_df.index:
            continue
        iso = mapping_df.loc[basename,'Isolate'];
        pmn = int(key.split('PM')[1][0]);fl
        rep = [int(key.split('-')[-1]) if '-' in key else 1][0];
        ribo = mapping_df.loc[basename,'Ribotype'];
        plate_list_df.loc[key] = [iso,pmn,rep,ribo];
    
    plate_list_df = plate_list_df.sort_values(['Isolate','PM','Rep','Ribotype'],ascending=True);
    plate_list_df = plate_list_df.dropna(axis=0)

    plate_list_df = plate_list_df.drop(flagged_plates)
    
    return plate_list_df

def formatPlateData(df,interval=600):
    
    # convert headers from strings to integers
    df.columns = plates.listTimePoints(interval=interval,numTimePoints=df.shape[1])
    
    # maintain name for index column and headers row
    df.index.name = 'Well'
    df.T.index.name = 'Time'
    
    # remove columns (time points) with only NA vlaues
    df = df.iloc[:,np.where(~df.isna().all(0))[0]]
    df = df.astype(float)

    # overrides time to intervals of 600
    # sets index names to Well and Time respectively
    # gets raid of NA values
    # confirms that it is float

    return df

def summarizePlateData(plate_id,plate_data):

    summary = plates.summarizeGrowthData(plate_data,expand_well_id=True)
    
    plate_key = plates.initializePlateKey(plate_id)
    plate_key = plate_key.join(summary)

    return plate_key

def summarizeMultiplePlates(data_dict,plate_list):
    
    summary_dict = {};
    
    for plate_id in plate_list:
        
        print 'Summarizing %s' % plate_id
        
        summary_dict[plate_id] = summarizePlateData(plate_id,data_dict[plate_id]);
        
    return summary_dict

def preModellingDataFormatting(data_dict,summary_dict,plate_list):
    
    new_data_dict = {};
    new_summary_dict = {};
    
    for plate_id in plate_list:
        
        plate_data = data_dict[plate_id]
        plate_key = summary_dict[plate_id]
        
        plate_data = plate_data.T # transpose
        plate_data = plate_data.loc[:,plate_key.index] # make sure columns in same order as plate key rows
        plate_data = plate_data.reset_index(drop=False); # now we can reset index to a random number 
        
        new_data_dict[plate_id] = plate_data
        new_summary_dict[plate_id] = plate_key
        
    return new_data_dict, new_summary_dict
        
    
def modelPlateData(plate_data,plate_summary):
    
    params = ['GP_r','GP_td','GP_d','GP_K','GP_AUC','GP_max'];
   
    pred_data = plate_data.copy() 
    plate_data = growth.GrowthPlate(data=plate_data,key=plate_summary);
 
    plate_data.convertTimeUnits()
    plate_data.logData()
    plate_data.subtractBaseline()
    
    growth_summary_df = pd.DataFrame(index=plate_summary.index,columns=params);
    
    counter = 0;
    row_count = 1; print '%02d' % row_count,
    for well in growth_summary_df.index:
        
        if counter<12:
            print '.',
            counter += 1;
        else:
            row_count += 1;
            print '\n%02d .' % row_count,
            counter = 1;
        
        substrate = plate_summary.loc[well,'Substrate']
        well_data = plate_data.extractGrowthData({'Substrate':[substrate]});
        
        well_data = growth.GrowthMetrics(well_data);
        well_data.fitGP();
        well_data.inferGPDynamics();
        well_data.inferDoublingTime(mtype='GP');
        well_data.predictGP();
        
        to_header = params;
        to_index = well_data.key.index[0];
        
        growth_summary_df.loc[to_index,to_header] = [float(ii) for ii in well_data.key.loc[to_index,to_header].values];
        
        pred_data.loc[:,well] = well_data.pred;

    print 
    plate_summary = plate_summary.join(growth_summary_df.astype(float));

    return plate_summary,pred_data

def modelMultiplePlates(data_dict,summary_dict,plate_list):
        
    new_summary_dict = {};
    pred_data_dict = {};
    
    for plate_id in plate_list:
        
        print 'Modelling %s' % plate_id
        
        plate_data = data_dict[plate_id]
        plate_summary = summary_dict[plate_id]
        
        new_summary_dict[plate_id],pred_data_dict[plate_id] = modelPlateData(plate_data,plate_summary)
        
        filepath = "../results/%s/%s.txt" % (plate_summary.Isolate[0],
                                              plate_summary.Plate[0])
        
        plates.createFolder('../results/%s/' % plates.parsePlateName(plate_id)[0])
        
        if not os.path.isfile(filepath):
            new_summary_dict[plate_id].to_csv(filepath,sep='\t',header=True,index=True)
            
    return new_summary_dict,pred_data_dict


def visualCheck(data_dict,summary_dict,plate_list,save_dirname="../figures"):
    
    for plate_id in plate_list.index:
        
        print 'Plotting %s' % plate_id
        
        plate_data = data_dict[plate_id];
        plate_summary = summary_dict[plate_id]
    
        plate_obj = growth.GrowthPlate(data=plate_data,key=plate_summary)
        
        filepath = "%s/%s/%s.pdf" % (save_dirname,
        							plate_summary.Isolate[0],
                                    plate_summary.Plate[0])
        
        plates.createFolder('%s/%s/' % (save_dirname,plates.parsePlateName(plate_id)[0]))
        
        if not os.path.isfile(filepath):

            fig,ax = plate_obj.plot(savefig=True,filepath=filepath);

