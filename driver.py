#!/usr/bin/env python

# Firas Said Midani
# Start date: 2019-09-10
# Final date: 2019-09-10

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

args = parser.parse_args();

verbose = args.verbose

##################
# PROCESS INPUTS
##################


DIRECTORY,MAPPING,FILES = {},{},{};

def checkDirectoryExists(directory,generic_name='Directory',verbose=False,sys_exit=False,initialize=False):

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

def checkDirectoryNotEmpty(directory,generic_name='Data directory',verbose=verbose):

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
        return criteria

    elif len(criteria)==0:
        return criteria

    df = df[df.isin(criteria).sum(1)==len(criteria)];

    if verbose:
        print 'The following criteria was used to subset the data'
        print criteria

    return df

def smartMapping(filebase,mapping_path,DATA,df_meta,df_meta_plates):

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

    elif plates.isBIOLOG(filebase):

        df_mapping = plates.initializeBiologPlateKey(filebase)
        print '%s: Did not find file or meta-data but seems to be a BIOLOG plate' % filebase
       
    else:
        df_mapping = pd.DataFrame(index=well_ids,columns=['Plate_ID'],
                                  data = [filebase]*len(well_ids))
        
        print '%s Did not find file or meta-data and does not seem to be a BIOLOG plate' % filebase

    df_mapping = df_mapping.reset_index(drop=False)

    return df_mapping

print 'VERIFYING DIRECTORY STRUCTURE AND USER INPUT'
DIRECTORY['PARENT'] = args.input;
DIRECTORY['DATA'] = '%s/data' % DIRECTORY['PARENT']
DIRECTORY['DERIVED'] = '%s/data_derived' % DIRECTORY['PARENT']
DIRECTORY['MAPPING'] = '%s/mapping' % DIRECTORY['PARENT']
DIRECTORY['PARAMETERS'] = '%s/parameters' % DIRECTORY['PARENT']

FILES['META'] = '%s/meta.txt' % DIRECTORY['MAPPING']
FILES['FLAG'] = '%s/flag.txt' % DIRECTORY['PARAMETERS']
FILES['HYPO'] = '%s/hypothesis.txt' % DIRECTORY['PARAMETERS']
FILES['SUBSET'] = '%s/subset.txt' % DIRECTORY['PARAMETERS']

checkDirectoryExists(DIRECTORY['PARENT'],'Input directory',verbose=True,sys_exit=True)
checkDirectoryExists(DIRECTORY['DATA'],'Data directory',verbose=True)
checkDirectoryNotEmpty(DIRECTORY['DATA'],'Data directory',verbose)

print 'READING & ASSEMBLING DATA'
list_data = sorted(os.listdir(DIRECTORY['DATA']));
DATA = readPlateReaderFolder(folderpath=DIRECTORY['DATA'],save=True,save_dirname=DIRECTORY['DERIVED'])

print 'CHECKING FOR META.TXT'
df_meta, df_meta_plates = checkMetaTxt(FILES['META'],verbose=True)

print 'READING & ASSEMBLING MAPPING'
for filename in list_data:

    filebase = os.path.splitext(filename)[0];
    mapping_path = '%s/%s.txt' % (DIRECTORY['MAPPING'],filebase);
    well_ids = DATA[filebase].index;

    MAPPING[filebase] = smartMapping(filebase,mapping_path,well_ids,df_meta,df_meta_plates)
print 

MASTER_MAPPING = pd.concat(MAPPING.values(),ignore_index=True,sort=False)

print 'REMOVING FLAGGED WELLS'
flag_dict = checkDictTxt(FILES['FLAG'],verbose=True);
MASTER_MAPPING = dropFlaggedWells(MASTER_MAPPING,flag_dict,verbose=True)
MASTER_MAPPING.to_csv('%s/stitched_mapping.txt' % DIRECTORY['MAPPING'],sep='\t',header=True,index=True)
print

print 'SUBSETTING BASED ON USER INPUT'
subset_dict = checkDictTxt(FILES['SUBSET'],verbose=True);
MASTER_MAPPING = subsetWells(MASTER_MAPPING,subset_dict,verbose=True)
MASTER_MAPPING.to_csv('%s/stitched_mapping.txt' % DIRECTORY['MAPPING'],sep='\t',header=True,index=True)
print 

print('\n')



