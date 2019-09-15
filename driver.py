#!/usr/bin/env python

# Firas Said Midani
# Start date: 2019-09-10
# Final date: 2019-09-12

# DESCRIPTION driver script for analyzing microbial growth curves


##################################
# IMPORT OFF-THE-SHELF LIBRARIES
##################################

import os
import sys
import argparse
import pandas as pd

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

#parser.add_argument('-f','--flag',required=False)
#parser.add_argument('-s','--subset',required=False)
#parser.add_argument('-H','--hypothesis',required=False)
#parser.add_argument('-m','--merge-results',required=False)

# --plot-plates (only plot plate to visualize for any odd things)
# --plot-wells (plot each desired well)
# --plot-only

# --no-flags
# --no-subsets



args = parser.parse_args();

verbose = args.verbose

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

def smartmapping(filebase,mapping_path,DATA,df_meta,df_meta_plates):

    if os.path.exists(mapping_path):

        print '%s: Reading %s' % (filebase,mapping_path)
        df_mapping = pd.read_csv(mapping_path,sep='\t',header=0,index_col=0);

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
##############

##########################################
print 'READING & ASSEMBLING DATA'
list_data = sorted(os.listdir(directory['DATA']));
interval_dict = checkDictTxt(files['INTERVAL'],verbose=True); print interval_dict
data = readPlateReaderFolder(folderpath=directory['DATA'],save=True,save_dirname=directory['DERIVED'],interval=600,interval_dict=interval_dict)
##########################################

##########################################
print 'CHECKING FOR META.TXT'
df_meta, df_meta_plates = checkMetaTxt(files['META'],verbose=True)
##########################################

##########################################
print 'READING & ASSEMBLING mapping'
for filename in list_data:
    filebase = os.path.splitext(filename)[0];
    mapping_path = '%s/%s.txt' % (directory['MAPPING'],filebase);
    well_ids = data[filebase].columns[1:];
    mapping[filebase] = smartmapping(filebase,mapping_path,well_ids,df_meta,df_meta_plates)
master_mapping = pd.concat(mapping.values(),ignore_index=True,sort=False)
print 

##########################################

# plot_plates = True
# if plot_plates:
#     for pid in data.keys():
#         print data[pid].head()
#         filepath = '%s/%s.pdf' % (directory['FIGURES'],pid)
#         data[pid].Time = data[pid].
#         #growth.GrowthPlate(data[pid]).plot(title="",savefig=False,filepath="")
#         print mapping[pid]

# sys.exit('DONE')

##########################################
print 'REMOVING FLAGGED WELLS'
flag_dict = checkDictTxt(files['FLAG'],verbose=True);
master_mapping = dropFlaggedWells(master_mapping,flag_dict,verbose=True)
#master_mapping.to_csv('%s/stitched_mapping.txt' % directory['mapping'],sep='\t',header=True,index=True)
print

##########################################

##########################################
print 'SUBSETTING MAPPING & DATA BASED ON USER INPUT'
subset_dict = checkDictTxt(files['SUBSET'],verbose=True);
master_mapping = subsetWells(master_mapping,subset_dict,verbose=True).reset_index(drop=True)#; 
master_mapping.index.name='Sample_ID';
master_mapping = master_mapping.reset_index(drop=False);
print
print master_mapping

sub_data,sub_key = {},{}
sub_gdata,sub_gkey = {},{}

for pid in data.keys():

    # sub_mapping is wells by meta-data variables (including WEll, Plate_ID)
    # can be used for GrowthData as key
    
    sub_mapping = master_mapping[master_mapping.Plate_ID == pid]; 
    if sub_mapping.shape[0] > 0:
        
        wells = list(sub_mapping.Well.values)
        
        df_pid = data[pid].loc[:,['Time'] + wells]; # same format; just removing wells/rows
        df_pid_melt = meltData(df_pid,pid); # 

        sub_gdata[pid] = df_pid;
        sub_data[pid] = df_pid_melt;

        sub_key[pid] = sub_mapping;

        #gdata = growth.GrowthPlate(data=df_pid,key=sub_mapping);

print 
master_data = pd.concat(sub_data.values(),sort=False).reset_index()
master_data.to_csv('%s/stitched_data_input.txt' % directory['MAPPING'],sep='\t',header=True,index=False)

def packageGrowthPlate(sub_gdata,sub_key):

    ## ll is short for left and rr is short for right
    gplate_data = reduce(lambda ll,rr: pd.merge(ll,rr,on='Time',how='outer'),sub_gdata.values());
    gplate_data = gplate_data.sort_values(['Time']).reset_index(drop=True);
    gplate_data.columns = ['Time'] + range(gplate_data.shape[1]-1);
    gplate_key = pd.concat(sub_key.values()).reset_index(drop=True);

    gplate = growth.GrowthPlate(data=gplate_data,key=gplate_key);

    return gplate
    
##########################################

##########################################
print 'READING HYPOTHESIS'
hypo_dict = checkDictTxt(files['HYPO'],verbose=True,spliton='+'); print hypo_dict
#master_mapping.to_csv('%s/stitched_mapping.txt' % directory['mapping'],sep='\t',header=True,index=True)
print
##########################################


if len(hypo_dict) > 0 :
    #gplate = growth.GrowthPlate(data=gdata_input,key=gdata_key)
    gplate = packageGrowthPlate(sub_gdata,sub_key)
    gplate.convertTimeUnits()
    gplate.logData()
    gplate.subtractBaseline()

    gplate.runTestGP(hypothesis=hypo_dict)

    sys.exit('~~~DONE~~~')

#sys.exit('~~~DONE~~~')

##########################################

key_list = [];
od_list = [];

for _,sample in master_mapping.iterrows():

    Well_ID, Plate_ID = sample[['Well','Plate_ID']];
    
    # initialize data and key
    sample_data = grabCurve(master_data,Well_ID,Plate_ID);
    sample_key = pd.DataFrame([Well_ID,Plate_ID],index=['Well','Plate_ID']).T

    # package into GrwothData object
    curve = packGrowthData(sample_data,sample_key)

    # basic summary: min, max, and T0
    curve.basicSummary()

    # convert time from seconds to hours, base-10 log-transform OD, and subtract T0
    curve.convertTimeUnits()
    curve.logData()
    curve.subtractBaseline()

    # apply GP model fitting, infer growth parameters, and predict OD
    curve_fit = growth.GrowthMetrics(curve);
    curve_fit.fitGP();
    curve_fit.inferGPDynamics();
    curve_fit.inferDoublingTime(mtype='GP');
    curve_fit.predictGP();

    print curve_fit.key.head()

    curve_output = pd.concat([curve_fit.time,curve.input_data,curve_fit.data,curve_fit.pred],axis=1)
    curve_output.columns = ['Time','OD_input','OD_cleaned','OD_predicted']
    
    well_id = curve_fit.key.loc[:,['Well','Plate_ID']].values
    well_id_df = pd.DataFrame(index=curve_output.index,columns=['Well','Plate_ID']);
    well_id_df.loc[:,['Well','Plate_ID']] = [well_id]*curve_output.shape[0]

    curve_output = well_id_df.join(curve_output)

    key_list.append(curve_fit.key)
    od_list.append(curve_output)
    #pred_list.append(curve_fit.pred)

#     fig,ax = curve_fit.plot()
#     plt.savefig('/Users/firasmidani/Downloads/20190911/%s_%s.pdf' % (Plate_ID,Well_ID),filetype='pdf')

sys.exit('~~~DONE~~~')


key_list = pd.concat(key_list)
key_list.to_csv('%s/stitched_key_output.txt' % directory['SUMMARY'],sep='\t',header=True,index=False)

od_list = pd.concat(od_list)
od_list.to_csv('%s/stitched_od_output.txt' % directory['SUMMARY'],sep='\t',header=True,index=False)

##########################################

# do i combine all data (melt it too) ? 
# do i run each plate separately ? 
# do i save all data results together ? 

# I should analyze all seperately but have the option to merge all results!




print('\n')



