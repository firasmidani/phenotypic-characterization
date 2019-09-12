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

#########################
# IMPORT USER LIBRARIES
#########################

from libs import plates
from libs.pipeline import readPlateReaderFolder

######################################
# USER INTERFACING AND INPUT PARSING
######################################

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True)
parser.add_argument('-v','--verbose',action='store_true',default=False)


# -f command --> creates flag.txt file
# -s command --> creates substrate.txt file
# -h command --> creates hypothesis.txt file


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

def checkDictTxt(path,verbose=False):

    exists = os.path.exists(path);

    args = {};

    if exists:

        fid = open(path,'r');
        
        for line in fid:

            key,values = line.split(':');
            values = values.strip('\n').split(',');
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

##########################################
print 'VERIFYING directory STRUCTURE AND USER INPUT'
directory['PARENT'] = args.input;
directory['DATA'] = '%s/data' % directory['PARENT']
directory['DERIVED'] = '%s/data_derived' % directory['PARENT']
directory['mapping'] = '%s/mapping' % directory['PARENT']
directory['PARAMETERS'] = '%s/parameters' % directory['PARENT']

files['META'] = '%s/meta.txt' % directory['mapping']
files['FLAG'] = '%s/flag.txt' % directory['PARAMETERS']
files['HYPO'] = '%s/hypothesis.txt' % directory['PARAMETERS']
files['SUBSET'] = '%s/subset.txt' % directory['PARAMETERS']

checkdirectoryExists(directory['PARENT'],'Input directory',verbose=True,sys_exit=True)
checkdirectoryExists(directory['DATA'],'Data directory',verbose=True)
checkdirectoryNotEmpty(directory['DATA'],'Data directory',verbose)
##############

##########################################
print 'READING & ASSEMBLING DATA'
list_data = sorted(os.listdir(directory['DATA']));
data = readPlateReaderFolder(folderpath=directory['DATA'],save=True,save_dirname=directory['DERIVED'])
##########################################

##########################################
print 'CHECKING FOR META.TXT'
df_meta, df_meta_plates = checkMetaTxt(files['META'],verbose=True)
##########################################

##########################################
print 'READING & ASSEMBLING mapping'
for filename in list_data:
    filebase = os.path.splitext(filename)[0];
    mapping_path = '%s/%s.txt' % (directory['mapping'],filebase);
    well_ids = data[filebase].index;
    mapping[filebase] = smartmapping(filebase,mapping_path,well_ids,df_meta,df_meta_plates)
master_mapping = pd.concat(mapping.values(),ignore_index=True,sort=False)
print 
##########################################

##########################################
print 'REMOVING FLAGGED WELLS'
flag_dict = checkDictTxt(files['FLAG'],verbose=True);
master_mapping = dropFlaggedWells(master_mapping,flag_dict,verbose=True)
master_mapping.to_csv('%s/stitched_mapping.txt' % directory['mapping'],sep='\t',header=True,index=True)
print
##########################################

##########################################
print 'SUBSETTING MAPPING & DATA BASED ON USER INPUT'
subset_dict = checkDictTxt(files['SUBSET'],verbose=True);
master_mapping = subsetWells(master_mapping,subset_dict,verbose=True)
master_mapping.to_csv('%s/stitched_mapping.txt' % directory['mapping'],sep='\t',header=True,index=True)
for pid in data.keys():

    wells = master_mapping[master_mapping.Plate_ID == pid].Well
    data[pid] = data[pid].loc[wells,:]
print
##########################################



print('\n')



