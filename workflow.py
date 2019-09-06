#!/usr/env/bin python

##################################
# IMPORT OFF-THE-SHELF LIBRARIES
##################################
import os
import sys
import GPy, scipy
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import copy

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

#########################
# IMPORT USER LIBRARIES
#########################

#sys.path.append("/data/davidalb/users/fsm/biolog/phenotypic-characterization/")
from libs import classical,growth,pipeline,plates
#from config import mapping_df

##################
# PROCESS INPUTS
##################

# PROCESS INPUTS
# 	1. directory
#	2. data
# 	3. flags
# 	4. metatdata
#	5. configurations

# receive arguments (ERROR CHECK NO ARGUMENT ADDED)
parentPath = sys.argv[1]; print 'parent path: %s' % parentPath
#biolog = sys.argv[2]; print 'BIOLOG: %s' % biolog
#mappingPath = sys.argv[3]; print 'mapping path %s' % mappingPath

######################
# CHECK/INIT FOLDERS
######################

# initialize directories
req_folders = ['figures','mapping','results','summary','data_processed','data_formatted']
for folder in req_folders:
	plates.createFolder('%s/%s' % (parentPath,folder))

# initialize paths
dataPath = '%s/data' % parentPath
fmt_dirname = '%s/data_formatted' % parentPath
proc_dirname = '%s/data_processed' % parentPath
mappingPath = "%s/mapping" % parentPath

# ERROR CHECK: check if data is available and provided in the correct location
if len(os.listdir(dataPath)) == 0:
	print("WARNING: Data directory is empty")

##############
# COPY DATA
##############

# copy data files 
cmd = "cp %s/* %s" % (dataPath,proc_dirname); os.system(cmd)

####################
# CHECK/INIT FLAGS
####################

#check for flags (REPORT BACK whether flags are empty or not, how many plates flagged, how many wells flagged)
if os.path.exists("%s/mapping/flags.py" % parentPath):
	sys.path.append("%s/mapping" % parentPath)
	from flags import *
	### CHECK FLAGS FOR DATA FILES ONLY

else:
	print("WARNING: No flag.py detected therefore all plates and wells will be processed.")
	flagged_plates = [];
	flagged_wells = {}

##############
# READ DATA
##############

data_dict = pipeline.readPlateReaderFolder(folderpath=dataPath,save=True,save_dirname=fmt_dirname)

#################
# READ MAPPING
#################

# check if

# if os.path.exists(mappingPath):

# 	mapping_df = pipeline.readPlateReaderMapping('%s/summary/mapping_ribotype.txt' % parentPath,filter_plates=True)
# 	plate_list_df = pipeline.expandPlateMapping(data_dict,mapping_df,flagged_plates)

# if biolog:

# 	mapping_df = pipeline.initializePlateMapping(data_dict,flagged_plates)
# 	mapping_df.to_csv(mappingPath,sep='\t',header=True,index=True)

mapping_dict = {}

# STRUCTURE OF  MAPPING
# | WELL | Plate ID | Isolate | 

for filename in data_dict.keys():

	# ## CHECK FOR BIOLOG FILENAME NOMENCLATURE
	# if plates.isBIOLOG(filename) is not None:

	mapping_file_path = "%s/%s.txt" % (mappingPath,filename)

	# if mapping file already exists 
	if os.path.exists(mapping_file_path):
		mapping_dict[filename] = pd.read_csv(mapping_file_path,sep='\t',header=0,index_col=0)

	# if not, but data is BIOLOG, create manifest
	elif plates.isBIOLOG(filename): 

		mapping_df = plates.initializeBiologPlateKey(filename);
		mapping_df.to_csv(mapping_file_path,sep='\t',header=True,index=True);
		mapping_dict[filename] = mapping_df;

	# otherwise, gently yell at user
	else:
		error_msg = '\n'
		error_msg += 'WARNING: \'%s\' does not have a mapping file ' % filename
		error_msg += 'and does not follow BIOLOG PM plate nomenculature. '
		error_msg += 'If data file corresponds to a BIOLOG plate, '
		error_msg += 'please rename file using accepted nomenclature '
		error_msg += '{name}_PM{number}-{replicate} (e.g. ECOLI_PM1-2), '
		error_msg += 'where {name} indicates arbitrary plate ID '
		error_msg += '{number} indicates BIOLOG PM plate number, '
		error_msg += '{replicate} is a digit that indicates the replicate number for the plate. '
		error_msg += 'Otherwise, please create a mapping file and save it in the mapping directory. '
		error_msg += 'To continue, we will initialize a minimalist mapping file for the plate.'
		error_msg += '\n'

		print(error_msg)

		mapping_df = plates.parseWellLayout()
		mapping_df.loc[:,'Plate'] = [filename] * mapping_df.shape[0];
		mapping_df = mapping_df.drop(['Row','Column'],axis=1)
		mapping_df.to_csv(mapping_file_path,sep='\t',header=True,index=True);
		mapping_dict[filename] = mapping_df;



# 1. If NO MAPPING FILE & BIOLOG --> AUTOMATICALLY CREATE A MANIFEST FROM DATA FILENAME
# 2. IF NO MAPPING FILE & NOT BIOLOG --> AUTOMATICALLY CREATE A MINIMALIST MANIFEST
# 2. IF MAPPING FILE & NOT BIOLOG --> IMPORT AND CHECK FOR CORRECT FORMAT

# DISTINGUISH BETWEEN PLATE MAPPING AND ISOLATE MAPPING


# print
# print 'READ MAPPING'
# mapping_df = pipeline.readPlateReaderMapping('%s/summary/mapping_ribotype.txt' % parentPath,filter_plates=True)
# print
# print 'SELECT PLATES'
# plate_list_df = pipeline.expandPlateMapping(data_dict,mapping_df,flagged_plates)
# # plate_list_df = plate_list_df.loc[['PRB1054_PM1-1','PRB1054_PM1-2','PRB1054_PM2-1','PRB1054_PM2-2',
# # 								   'PRB1055_PM1-1','PRB1055_PM1-2','PRB1055_PM2-1','PRB1055_PM2-2',
# # 								   'PRB1056_PM1-1','PRB1056_PM1-2','PRB1056_PM2-1','PRB1056_PM2-2',
# # 								   'PRB1058_PM1-1','PRB1058_PM1-2','PRB1058_PM2-1','PRB1058_PM2-2',
# # 								   'PRB1059_PM1-1','PRB1059_PM1-2','PRB1059_PM2-1','PRB1059_PM2-2'],:]
# print
# print 'SUMMARIZE PLATES'
# summary_dict = pipeline.summarizeMultiplePlates(data_dict,plate_list_df.index)
# print
# print 'PREPARING DATA FOR GP REGRESSION'
# new_data_dict,new_summary_dict = pipeline.preModellingDataFormatting(data_dict,summary_dict,plate_list_df.index)
# print
# print 'VISUAL CHEK FOR ABERRANT WELLS OR PLATES'
# pipeline.visualCheck(new_data_dict,new_summary_dict,plate_list_df,save_dirname='%s/figures' % parentPath)
# print
# print 'PERFORMING GAUSSIAN PROCESS REGRESSING'
# new_summary_dict,pred_data_dict = pipeline.modelMultiplePlates(new_data_dict,new_summary_dict,plate_list_df.index)

# summary_df = pd.concat(new_summary_dict).reset_index(drop=True);
# data_df = pd.concat(new_data_dict,axis=1)
# pred_data_df = pd.concat(pred_data_dict,axis=1)

# timestamp = plates.getFormattedTime()

# summary_df.to_csv('%s/summary/summary_%s.txt' % (parentPath,timestamp),
#     sep='\t',header=True,index=True)
# data_df.to_csv('%s/summary/data_%s.txt' % (parentPath,timestamp),
#     sep='\t',header=True,index=True)
# pred_data_df.to_csv('%s/summary/pred_%s.txt' % (parentPath,timestamp),
#     sep='\t',header=True,index=True)

#CLEAN_UP
#cmd = "rm %s/*" % dataPath
