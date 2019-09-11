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

from libs.plates import breakDownFilePath,initializeBiologPlateKey
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


DIRECTORY,FILES = {},{};

# def checkDirectoryExists(directory,generic_name='Directory',verbose=verbose):

#     if not os.path.exists(directory):
#         sys.exit('USER ERROR: %s %s does not exist' % (generic_name,directory))

#     elif verbose:
#         print '%s is %s' % (generic_name,directory)

#     return True

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
#BOOL_meta = checkDirectoryExists(FILES['MAPPING'],generic_name='Metadata file',verbose=True,initialize=True)

dict_mapping = {};
BOOL_meta = os.path.exists(FILES['META'])
#check format of mapping file. Must have PLATE_ID.

# read meta.txt or initialize
if BOOL_meta:
    df_meta = pd.read_csv(FILES['META'],sep='\t',header=0,index_col=None);
    print df_meta

for filename in list_data:
  
    filebase = os.path.splitext(filename)[0]
    mapping_path = '%s/%s.txt' % (DIRECTORY['MAPPING'],filebase)

    data = DATA[filebase];
    nwells = data.shape[0];

    # if mapping file exists, read it
    if os.path.exists(mapping_path):
   
        print 'Looking for %s' % mapping_path
    
        df_mapping = pd.read_csv(mapping_path,sep='\t',header=0,index_col=0);

    # if PLATE_ID is found in meta.txt, create and populate a mapping file using meta.txt info
    elif filebase in df_meta.Plate_ID.values:

        metadata = df_meta[df_meta.Plate_ID==filebase]
        biolog = metadata.PM in range(1,7)

        if biolog:
            isolate = metadata.isolate[0];
            pm = metadata.PM[0];
            rep = metadata.Replicate[0];

            filebase = '%s_PM%s-%s' % (isolate,pm,rep);
            df_mapping = plates.initializeBiologPlateKey(filebase);;
        else:
            df_mapping = pd.concat([df_meta[df_meta.Plate_ID==filebase]]*nwells,axis=0)
            df_mapping.index = data.index;

    # if filename follows Biolog nomenclature, create and populate a mapping file accoridngly
    elif isBiolog(filebase):

        # iso,pmn,rep = parsePlateName(filebase,simple=False);

        # df_mapping = pd.DataFrame(columns=['Plate_ID','Isolate','PM','Replicate'],
        #                           index=data.index,
        #                           data=[filebase,iso,pmn,rep]*nwells);
        df_mapping = initializeBiologPlateKey(filebase)

    # otherwise, create a bare-bone mapping file with single column of PLATE_ID
    else:
        df_mapping = pd.DataFrame(index=data.index,columns=['Plate_ID'],
                                  data = [filebase]*nwells)

    MAPPING[filebase] = df_mapping

print('\n')



