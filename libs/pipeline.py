#!/usr/env/bin python

# import off-the-shelf libraries
import re
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

def initializeParameter(arg_files,arg_in,arg_value,sep=',',integerize=False):
    '''

    ARGs:
        arg_files (dictionary) -- keys are parameters, values are location of paramter files
        arg_in  (str) -- argument parsed by script
        arg_value (str) -- parameter, e.g. META, FLAG, HYPO, SUBSET, or INTERVAL
        sep (str) -- delimiter
        integerize (boolean) -- whether to conver values in a list to integers
    '''

    if arg_in is None:
        arg_dict = checkDictTxt(arg_files[arg_value],verbose=True);
    elif len(arg_in)>0:
        arg_dict = checkArgText(arg_in,sep=sep);
    else:
        arg_dict = {};

    if integerize:
        return integerizeDictValues(arg_dict);
    else:
        return arg_dict


def checkDictTxt(path,verbose=False,spliton=','):

    exists = os.path.exists(path);

    args = {};

    if exists:

        fid = open(path,'r');
        
        for line in fid:

            key,values = line.split(':');
            values = values.strip('\n').split(spliton);
            values = [ii.strip() for ii in values];
            values = [float(ii) if ii.isdigit() else ii for ii in values]
            args[key] = values

    return args

def checkMetaTxt(path,verbose=False):

    BOOL_meta = os.path.exists(path)

    if not BOOL_meta:
        df_meta = pd.DataFrame
        df_meta_plates = [];
    else:
        df_meta = pd.read_csv(path,sep='\t',header=0,index_col=None);

    try: 
        df_meta_plates = df_meta.Plate_ID.values
    except:
        df_meta_plates = [];

    if verbose and BOOL_meta: 
        print 
        print 'Meta-Data is'
        print df_meta
        print 
    elif verbose:
        print
        print 'No meta.txt file found.'
        print

    return df_meta,df_meta_plates

def checkArgText(command,sep=','):
    '''Parses command-line text and formats as a dictionary.

    * text is a list of items that are separated by semicolons (;)
    * each item is a variable name separated from a list of values with a colon (:)
    * each list has values separated by commmas (,)

    Args:
        command (str)

    Returns:
        dictionary with values as lists

    example ARGUMENT 'Isolate:CD630,PRB599;Substrate:Negative Control,D-Trehalose'
    example OUTPUT {'Isolate':['CD630','PRB599],'Substrate':['Negative Control','D-Trehalose]}
    '''
    if command is None:
        return None

    # strip flanking semicolons or whitespaces then split by semicolon
    lines = command.strip(';').strip(' ').split(';'); 

    # get names of variables
    lines_keys = [ii.split(':')[0] for ii in lines];

    # get list of values for all variables
    lines_values = [re.split(sep,ii.split(':')[1]) for ii in lines]; # why did I use regex here?
    
    # re-package variables and their list of values into a dictionary
    lines_dict = {ii:jj for ii,jj in zip(lines_keys,lines_values)};

    return lines_dict

def integerizeDictValues(dict):
    '''
    converts values in a dictionary into integers. this works if values are iterables (e.g. list)

    Args:
        command (str)

    Returns:
        dictionary with values as lists of integers

    example ARGUMENT {'CD630_PM1-1':str(500)}
    example OUTPUT {'CD630_PM1-1':int(500)}
    '''
    if dict is None:
        return None

    return {key:[int(vv) for vv in value] for key,value in dict.iteritems()}

def checkDirectoryExists(directory,generic_name='directory',verbose=False,sys_exit=False,initialize=False):
    '''
    Checks if a directory exists. If directory does not exist, it could be initialized. 

    Args:
        directory (str) -- path
        generic_name (str) -- used in communication with users
        verbose (boolean)
        sys_exit (boolean) -- used to exit in case of error or severe WARNING
        initialize (boolean)

    Returns:
        boolean
    '''
    exists = os.path.exists(directory);

    if exists and verbose:
        print('%s is %s\n' % (generic_name,directory))
        return True

    elif exists: 
        return True

    elif sys_exit:
        sys.exit('USER ERROR: %s %s does not exist.\n' % (generic_name,directory))

    elif initialize and verbose:
        os.makedirs(directory)
        print('WARNING: %s did not exist but was created.\n' % (directory))
        return True

    elif initialize:
        os.makedirs(directory)
        return True

    elif verbose:
        print('WARNING: %s does not exist.\n' % (directory))
        return False

    else:
        return False

def checkDirectoryNotEmpty(directory,generic_name='Data directory',verbose=False):
    '''
    checks that directory is not empty

    Args:
        directory (str) -- path
        generic_name (str) -- used in communication with users
        verbose (boolean)
    '''
    if len(os.listdir(directory)) == 0:
        sys.exit('USER ERROR: %s %s is empty.\n' % (generic_name,directory))

    elif verbose:
        print('%s %s has %s files:' % (generic_name,directory,len(os.listdir(directory))))
        printItemizedList(directory)
        print('\n')

    return True

def dropFlaggedWells(df,flags,verbose=False):
    '''

    ARGS:
        df (pandas.DataFrame) must have Plate_IDs and Well as columns
        flags (dictionary) with Plate_IDs (str) as keys and Wells (stR) as vlaues
        verbose (boolean)

    Returns: 
        df (pandas.DataFrame)
    '''

    if (len(flags)==0) and (verbose):
        print 'No flag.txt was found.'
        return df

    elif len(flags)==0:
        return df

    for plate_id,wells in flags.iteritems():

        mapping_dict = {'Plate_ID':[plate_id],'Well':wells};
        to_drop = df[df.isin(mapping_dict).sum(1)== len(mapping_dict)].index;
        df = df.drop(labels=to_drop)

    if verbose:
        print 'Following flags were detected'
        print flags

    return df

def subsetWells(df,criteria,verbose=False):
    ''' subsets form mapping file'''

    if (len(criteria)==0) and (verbose):
        print 'No subsetting was requested.'
        return df

    elif len(criteria)==0:
        return df

    df = df[df.isin(criteria).sum(1)==len(criteria)];

    if verbose:
        print 'The following criteria was used to subset the data'
        print criteria

    return df

def smartmapping(filebase,mapping_path,well_ids,df_meta,df_meta_plates):

    if os.path.exists(mapping_path):

        print '%s: Reading %s' % (filebase,mapping_path)
        df_mapping = pd.read_csv(mapping_path,sep='\t',header=0,index_col=0);

        if 'Plate_ID' not in df_mapping.index:

            df_mapping.loc[:,'Plate_ID'] = [filebase]*df_mapping.shape[0];

    elif filebase in df_meta_plates:
        
        print '%s: Found meta-data in meta.txt' % (filebase),

        metadata = df_meta[df_meta.Plate_ID==filebase] # pd.DataFrame
        biolog = plates.isBiologFromMeta(metadata);

        if biolog:
            df_mapping = plates.expandBiologMetaData(metadata)
            print 'and seems to be a BIOLOG plate' 

        else:
            df_mapping = plates.initKeyFromMeta(metadata,well_ids)
            print 'and does not seem to be a BIOLOG plate'

    elif plates.isBiologFromName(filebase):

        df_mapping = plates.initializeBiologPlateKey(filebase)
        print '%s: Did not find file or meta-data but seems to be a BIOLOG plate' % filebase
       
    else:
        df_mapping = pd.DataFrame(index=well_ids,columns=['Plate_ID'],
                                  data = [filebase]*len(well_ids))
        
        print '%s Did not find file or meta-data and does not seem to be a BIOLOG plate' % filebase

    df_mapping = df_mapping.reset_index(drop=False)

    return df_mapping

def meltData(df,plate_id):
    '''
    no index columns
    first column is wells 
    rest of header is time-points so ['Well',0,600,...]
    values are numeric
    '''

    #df.T.index.name = 'Time'

    #if df.index.name == 'Well':
    #    df = df.reset_index();

    #df = pd.melt(df,id_vars='Well',value_name='OD');
    df = df.melt(id_vars='Time',var_name='Well',value_name='OD')
    df.Time = df.Time.astype(float);
    df.loc[:,'Plate_ID'] = [plate_id]*df.shape[0];

    # Well Time OD Plate_ID
    # ...

    # one plate is 96 wells x 100 time points = 9,600 rows x 4 columns 

    return df

def grabCurve(df,Well_ID,Plate_ID):

    criteria = {'Well':[Well_ID],'Plate_ID':[Plate_ID]};
    df = df[df.isin(criteria).sum(1)==2];
    df = df.sort_values(['Plate_ID','Well','Time']);

    #df_time = df.Time;
    #df_OD = df.OD;

    return df

def grabVariable(df,varb):

    return pd.DataFrame(df.loc[:,varb])

def packGrowthData(df,key):
    
    time = grabVariable(df,'Time')
    od = grabVariable(df,'OD') 

    return growth.GrowthData(time,od,key)

def subsetCombineData(data,master_mapping):

    wide_data_dict = {}
    key_dict= {}

    # for each plate_id
    for pid in data.keys():

        # sub_mapping is wells by meta-data variables (including WEll, Plate_ID)
        # can be used for GrowthData as key
        
        # grab all wells and their info
        sub_mapping_df = master_mapping[master_mapping.Plate_ID == pid]; 

        # only continue if mapping is not empty
        if sub_mapping_df.shape[0] > 0:
            
            # list all well IDs (usually A1 through H12)
            wells = list(sub_mapping_df.Well.values)
            
            # recall that data[pid] is time (p) x wells (n)
            # so here you are selecting for wells in sub_mapping
            key_dict[pid] = sub_mapping_df;
            wide_data_dict[pid] = data[pid].loc[:,['Time'] + wells]; 

    return wide_data_dict,key_dict

def packageGrowthPlate(data_dict,key_dict):

    ## ll is short for left and rr is short for right
    gplate_data = reduce(lambda ll,rr: pd.merge(ll,rr,on='Time',how='outer'),data_dict.values());
    gplate_data = gplate_data.sort_values(['Time']).reset_index(drop=True);
    gplate_key = pd.concat(key_dict.values())#.set_index(['Sample_ID']).reset_index(drop=True);
    gplate_data.columns = ['Time'] + list(gplate_key.Sample_ID.values);

    gplate = growth.GrowthPlate(data=gplate_data,key=gplate_key);

    return gplate

def readPlateReaderFolder(folderpath,save=False,save_dirname='../data_formatted',interval=600,interval_dict={}):
    
    df_dict = {};

    filepaths = plates.findPlateReaderFiles(folderpath)

    for filepath in sorted(filepaths):

        _, filebase, newfile = plates.breakDownFilePath(filepath,save_dirname);

        if filebase in interval_dict.keys():
            plate_interval = interval_dict[filebase][0]
        else:
            plate_interval = interval

        df = plates.readPlateReaderData(filepath,save=save,save_dirname=save_dirname,interval=plate_interval);

        # Time A1 A2 ... H12
        # 0 0.1 0.102 ... 0.09 
        # ...
        # 6000 1.1 1.05 ... 1.0

        df_dict[filebase] = df;

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


def printItemizedList(directory,sort=True):
    '''
    prints a list of the contents in a directory

    ARGs:
        directory (str)
        sort (boolean)

    Returns:
        None
    '''
    if sort: 
        items = sorted(os.listdir(directory));

    for item in items:
        print item

    return None

def visualCheck(data_dict,summary_dict,plate_list,save_dirname="../figures"):
    
    for plate_id in plate_list.index:
        
        print 'Plotting %s' % plate_id
        
        plate_data = data_dict[plate_id];
        plate_summary = summary_dict[plate_id]
    
        plate_obj = growth.GrowthPlate(data=plate_data,key=plate_summary)
        
        filepath = "%s/%s/%s.pdf" % (save_dirname,
        							plate_summary.Isolate[0],
                                    plate_summary.Plate_ID[0])
        
        plates.createFolder('%s/%s/' % (save_dirname,plates.parsePlateName(plate_id)[0]))
        
        if not os.path.isfile(filepath):

            fig,ax = plate_obj.plot(savefig=True,filepath=filepath);

