#!/usr/bin/env python

# Firas Said Midani
# Start date: 2019-09-10
# Final date: 2019-09-15

# DESCRIPTION driver script for analyzing microbial growth curves


##################################
# IMPORT OFF-THE-SHELF LIBRARIES
##################################

import os
import re
import sys
import argparse
import numpy as np
import pandas as pd

import time

import matplotlib.pyplot as plt

#########################
# IMPORT USER LIBRARIES
#########################

from libs import plates,growth
from libs.pipeline import readPlateReaderFolder

######################################
# USER INTERFACING AND INPUT PARSING
######################################

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True)
parser.add_argument('-v','--verbose',action='store_true',default=False)

parser.add_argument('-f','--flag',required=False)
parser.add_argument('-s','--subset',required=False)
parser.add_argument('-H','--hypothesis',required=False)
parser.add_argument('-I','--interval',required=False)

parser.add_argument('--plot-plates-only',action='store_true',default=False)

#parser.add_argument('-m','--merge-results',required=False)

# --plot-wells (plot each desired well)
# --plot-only
# --no-flags
# --no-subsets

args = parser.parse_args();

flag = args.flag;
subset = args.subset;
verbose = args.verbose;
interval = args.interval;
hypothesis = args.hypothesis;
plot_plates_only = args.plot_plates_only # need to check that plates are 8x12

def checkArgText(command,sep=','):

    if command is None:
        return None

    lines = command.strip(';').strip(' ').split(';');
    lines_keys = [ii.split(':')[0] for ii in lines];
    lines_values = [re.split(sep,ii.split(':')[1]) for ii in lines];
    lines_dict = {ii:jj for ii,jj in zip(lines_keys,lines_values)};

    return lines_dict

def integerizeDictValues(dict):

    if dict is None:
        return None

    return {key:[int(vv) for vv in value] for key,value in dict.iteritems()}

##################
# PROCESS INPUTS
##################

directory,mapping,files = {},{},{};

def checkdirectoryExists(directory,generic_name='directory',verbose=False,sys_exit=False,initialize=False):

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

    elif initilaize:
        os.makedirs(directory)
        return True

    elif verbose:
        print('WARNING: %s does not exist.\n' % (directory))
        return False

    else:
        return False

def checkdirectoryNotEmpty(directory,generic_name='Data directory',verbose=verbose):

    if len(os.listdir(directory)) == 0:
        sys.exit('USER ERROR: %s %s is empty.\n' % (generic_name,directory))

    elif verbose:
        print('%s %s has %s files:' % (generic_name,directory,len(os.listdir(directory))))
        printItemizedList(directory)
        print('\n')

    return True

def printItemizedList(directory,sort=True):

    if sort: 
        items = sorted(os.listdir(directory));

    for item in items:
        print item

    return None


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

def dropFlaggedWells(df,flags,verbose=False):
    '''
    df is mapping DataFrame
    flags is a dictionary with Plate_IDs as keys and Wells as vlaues
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


##########################################
print 'VERIFYING directory STRUCTURE AND USER INPUT'
directory['PARENT'] = args.input;
directory['DATA'] = '%s/data' % directory['PARENT']
directory['DERIVED'] = '%s/data_derived' % directory['PARENT']
directory['MAPPING'] = '%s/mapping' % directory['PARENT']
directory['SUMMARY'] = '%s/summary' % directory['PARENT']
directory['PARAMETERS'] = '%s/parameters' % directory['PARENT']
directory['FIGURES'] = '%s/figures' % directory['PARENT']

files['META'] = '%s/meta.txt' % directory['MAPPING']
files['FLAG'] = '%s/flag.txt' % directory['PARAMETERS']
files['HYPO'] = '%s/hypothesis.txt' % directory['PARAMETERS']
files['SUBSET'] = '%s/subset.txt' % directory['PARAMETERS']
files['INTERVAL'] = '%s/interval.txt' % directory['PARAMETERS']

checkdirectoryExists(directory['PARENT'],'Input directory',verbose=True,sys_exit=True)
checkdirectoryExists(directory['DATA'],'Data directory',verbose=True)
checkdirectoryNotEmpty(directory['DATA'],'Data directory',verbose)
checkdirectoryExists(directory['DERIVED'],'Derived data directory',verbose=True,initialize=True)
checkdirectoryExists(directory['MAPPING'],'Summary directory',verbose=True,initialize=True)
checkdirectoryExists(directory['SUMMARY'],'Summary directory',verbose=True,initialize=True)
checkdirectoryExists(directory['FIGURES'],'Figures directory',verbose=True,initialize=True)

##############

##########################################
print 'READING & ASSEMBLING DATA'
list_data = sorted(os.listdir(directory['DATA']));

def initializeParameter(arg_in,arg_value,sep=',',integerize=False):

    if arg_in is None:
        arg_dict = checkDictTxt(files[arg_value],verbose=True);
    elif len(arg_in)>0:
        arg_dict = checkArgText(arg_in,sep=sep);
    else:
        arg_dict = {};

    if integerize:
        return integerizeDictValues(arg_dict);
    else:
        return arg_dict


interval_dict = initializeParameter(interval,'INTERVAL',sep=',',integerize=True)
print interval_dict
data = readPlateReaderFolder(folderpath=directory['DATA'],save=True,save_dirname=directory['DERIVED'],interval=600,interval_dict=interval_dict)
##########################################

if plot_plates_only:
    # check that each plate has 96 well IDs (in 8x12) format

    for fname,table in data.iteritems():

        print fname
        print table.head()
        actual_well_ids = set(sorted(table.columns.values[1:])); print actual_well_ids


        # if a plate is Biolog, it seem to automatically have 96 wells internally even if rows are missing
        expected_well_ids = set(sorted(plates.parseWellLayout().index.values)); print expected_well_ids

        print expected_well_ids.difference(actual_well_ids)

    sys.exit()

##########################################
print 'CHECKING FOR META.TXT'
df_meta, df_meta_plates = checkMetaTxt(files['META'],verbose=True)
##########################################

##########################################
print 'READING & ASSEMBLING mapping'
for filename in list_data:
    filebase = os.path.splitext(filename)[0]; 
    mapping_path = '%s/%s.txt' % (directory['MAPPING'],filebase); print mapping_path
    well_ids = data[filebase].columns[1:];
    mapping[filebase] = smartmapping(filebase,mapping_path,well_ids,df_meta,df_meta_plates)
master_mapping = pd.concat(mapping.values(),ignore_index=True,sort=False)
print 

##########################################

##########################################
print 'READING HYPOTHESIS'

hypo_dict = initializeParameter(hypothesis,'HYPO',sep='\+|,',integerize=False)
print hypo_dict 
#master_mapping.to_csv('%s/stitched_mapping.txt' % directory['mapping'],sep='\t',header=True,index=True)
print
##########################################

##########################################
print 'REMOVING FLAGGED WELLS'

flag_dict = initializeParameter(flag,'FLAG',sep=',',integerize=False)
print flag_dict
master_mapping = dropFlaggedWells(master_mapping,flag_dict,verbose=True);
#master_mapping.to_csv('%s/stitched_mapping.txt' % directory['mapping'],sep='\t',header=True,index=True)
print


##########################################

##########################################
print 'SUBSETTING MAPPING & DATA BASED ON USER INPUT'
subset_dict = initializeParameter(subset,'SUBSET',sep=',',integerize=False)
print subset_dict

master_mapping = subsetWells(master_mapping,subset_dict,verbose=True).reset_index(drop=True)#; 
master_mapping.index.name='Sample_ID';
master_mapping = master_mapping.reset_index(drop=False); 
print master_mapping

if master_mapping.shape[0] == 0:
    print('ERROR: no wells selected. System Exitting!')
    sys.exit()

print 'WIDENING'
wide_data_dict,key_dict = subsetCombineData(data,master_mapping)
gplate = packageGrowthPlate(wide_data_dict,key_dict)
print

##########################################

##########################################

gplate.computeFoldChange()
gplate.convertTimeUnits()
gplate.logData()
gplate.subtractBaseline()

#gplate.addRowColVarbs()

visual_check = False
if visual_check:

    gplate.addRowColVarbs()
    filepath = '%s/plot_visual_check.pdf' % directory['FIGURES']
    fig,axes = gplate.plot(savefig=True,title="",filepath=filepath,modified=True)
    
    sys.exit('DONE')

if len(hypo_dict) > 0 :
    print hypo_dict

    start = time.time()
    gplate.runTestGP(hypothesis=hypo_dict,thinning=9)
    print
    print('ELAPSED TIME:\t%s' % (time.time() - start))

    sys.exit('~~~DONE~~~')

sys.exit()


sample_output_list = [];
sample_mapping_list = [];

for sample_id in sorted(master_mapping.Sample_ID.unique()):

    #sample_id = 32

    sample_curve = gplate.extractGrowthData({'Sample_ID':[sample_id]});

    sample_metrics = growth.GrowthMetrics(sample_curve);

    # check if basicSummary already exists

    sample_metrics.basicSummary(unmodified=True);
    sample_metrics.fitGP();
    sample_metrics.inferGPDynamics();
    sample_metrics.inferDoublingTime(mtype='GP');
    sample_metrics.predictGP()

    # input_time is in seconds, time used is in hours for analysis but not output
    sample_output = pd.concat([sample_metrics.time,
                               sample_metrics.input_data,
                               sample_metrics.data,
                               sample_metrics.pred],axis=1)

    sample_output.columns = ['Time','OD_input','OD_cleaned','OD_predicted']
    sample_id_df = pd.DataFrame([sample_id]*sample_output.shape[0],columns=['Sample_ID']);
    sample_output = sample_output.join(sample_id_df)

    sample_output_list.append(sample_output)
    sample_mapping_list.append(sample_metrics.key);

#    if plot_results:

    # if sample_id == 32:

    #     filepath = '%s/plot_results.pdf' % directory['FIGURES']

    #     #print sample_metrics.key.T
    #     #print type(sample_metrics.key.Row)
    #     #print sample_metrics.key.Row.values[0]
    #     #print sample_metrics.key.Row.values[0]-1

    #     r = sample_metrics.key.Row.values[0]-1
    #     c = sample_metrics.key.Column.values[0]-1

    #     ax = axes[r,c]
    #     x_data = np.ravel(sample_metrics.time.astype(float).values)
    #     y_data = np.ravel(sample_metrics.pred.astype(float).values)
    #     y_data = [1]*len(x_data)
    #     print ax.get_xlim(),ax.get_yowlim()
    #     print x_data,y_data
    #     filepath = '%s/plot_result.pdf' % directory['FIGURES']

    #     ax.plot(x=[0,3],y=[0,3],color=(1,0,0,0.8),lw=10)
    #     axes[0,0].plot(x=x_data,y=    y_data,color=(1,1,0,0.8),lw=10)

    #     fig.savefig(filepath,filetype='pdf')
        
    #     sys.exit('~~~DONE~~~~')


sample_output_df = pd.concat(sample_output_list,axis=0)
sample_mapping_list = pd.concat(sample_mapping_list,axis=0);

OD_input = sample_output_df.pivot(index='Time',columns='Sample_ID',values='OD_input')
OD_cleaned = sample_output_df.pivot(index='Time',columns='Sample_ID',values='OD_cleaned')
OD_predicted = sample_output_df.pivot(index='Time',columns='Sample_ID',values='OD_predicted')

sample_output_df.to_csv('%s/stitched_od_output.txt' % directory['SUMMARY'],sep='\t',header=True,index=False)
sample_mapping_list.to_csv('%s/stitched_key_output.txt' % directory['SUMMARY'],sep='\t',header=True,index=False)


##########################################

print('\n')

sys.exit('~~~DONE~~~')


