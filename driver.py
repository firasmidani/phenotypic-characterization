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

# check if PARENT_DIRECTORY exists + DATA_DIRECOTRY exists & is not empty

DIRECTORY['PARENT'] = args.input;
DIRECTORY['DATA'] = '%s/data' % DIRECTORY['PARENT']
DIRECTORY['DERIVED'] = '%s/data_derived' % DIRECTORY['PARENT']
DIRECTORY['MAPPING'] = '%s/mapping' % DIRECTORY['PARENT']

FILES['META'] = '%s/meta.txt' % DIRECTORY['MAPPING']

print('\n')
checkDirectoryExists(DIRECTORY['PARENT'],'Input directory',verbose=True,sys_exit=True)
checkDirectoryExists(DIRECTORY['DATA'],'Data directory',verbose=True)
checkDirectoryNotEmpty(DIRECTORY['DATA'],'Data directory',verbose)

# initialize a list of all input files
list_data = sorted(os.listdir(DIRECTORY['DATA']));
DATA = readPlateReaderFolder(folderpath=DIRECTORY['DATA'],save=True,save_dirname=DIRECTORY['DERIVED'])

# look for meta.txt
BOOL_meta = os.path.exists(FILES['META'])
#check format of mapping file. Must have PLATE_ID.

# read meta.txt or initialize
if BOOL_meta:
    df_meta = pd.read_csv(FILES['META'],sep='\t',header=0,index_col=None);
    df_meta_plates = df_meta.Plate_ID.values;

    print 'Meta-Data is'
    print df_meta

for filename in list_data:
  
    filebase = os.path.spliext(filename)[0]
    mapping_path = '%s/%s.txt' % (DIRECTORY['MAPPING'],filebase)

    data = DATA[filebase];
    well_ids = data.index;

    if os.path.exists(mapping_path):
   
        print 'Reading %s' % mapping_path
        df_mapping = pd.read_csv(mapping_path,sep='\t',header=0,index_col=0);

    elif filebase in df_meta_plates:

        print 'Found Meta-Data for %s in meta.txt' % filebase

        metadata = df_meta[df_meta.Plate_ID==filebase] # pd.DataFrame

        biolog = plates.isBiologFromMeta(metadata);

        if biolog:
            df_mapping = plates.expandBiologMetaData(metadata)
            print 'and %s seems to be a BIOLOG plate' % filebase
        else:
            df_mapping = plates.initKeyFromMeta(metadata,well_ids)
            print 'and %s does not seem to be a BIOLOG plate' % filebase

    elif plates.isBIOLOG(filebase):

        print filebase
        df_mapping = plates.initializeBiologPlateKey(filebase)

        print 'Did not find file or meta-data but %s seems to be a BIOLOG plate' % filebase

    else:
        df_mapping = pd.DataFrame(index=well_ids,columns=['Plate_ID'],
                                  data = [filebase]*len(well_ids))
        
        print 'Did not find file or meta-data and %s does not seems to be a BIOLOG plate' % filebase

    MAPPING[filebase] = df_mapping.reset_index(drop=False)

df_mapping = pd.concat(MAPPING.values(),ignore_index=True,sort=False)
df_mapping.to_csv('%s/stitched_mapping.txt' % DIRECTORY['MAPPING'],sep='\t',header=True,index=True)

print('\n')



