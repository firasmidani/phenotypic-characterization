#!/usr/bin/env python

# Firas Said Midani
# Start date: 2019-09-10
# Final date: 2019-10-14

# DESCRIPTION driver script for analyzing microbial growth curves


##################################
# IMPORT OFF-THE-SHELF LIBRARIES
##################################

import os
import sys
import argparse
import numpy as np
import pandas as pd

import time

import matplotlib.pyplot as plt

#########################
# IMPORT USER LIBRARIES
#########################

from libs import plates,growth,pipeline
#from
#from libs.pipeline import readPlateReaderFolder

######################################
# USER INTERFACING AND INPUT PARSING
######################################

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True)

parser.add_argument('-f','--flag',required=False)
parser.add_argument('-s','--subset',required=False)
parser.add_argument('-H','--hypothesis',required=False)
parser.add_argument('-I','--interval',required=False)

parser.add_argument('-v','--verbose',action='store_true',default=False)

#parser.add_argument('--plot-plate-only',action='store_true',default=False)

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
#plot_plate_only = args.plot_plate_only # need to check that plates are 8x12

directory,mapping,files = {},{},{};


print '#############################################'
print 'VERIFYING directory STRUCTURE AND USER INPUT'
print 
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

pipeline.checkDirectoryExists(directory['PARENT'],'Input directory',verbose=True,sys_exit=True)
pipeline.checkDirectoryExists(directory['DATA'],'Data directory',verbose=True)
pipeline.checkDirectoryNotEmpty(directory['DATA'],'Data directory',verbose)
pipeline.checkDirectoryExists(directory['DERIVED'],'Derived data directory',verbose=True,initialize=True)
pipeline.checkDirectoryExists(directory['MAPPING'],'Summary directory',verbose=True,initialize=True)
pipeline.checkDirectoryExists(directory['SUMMARY'],'Summary directory',verbose=True,initialize=True)
pipeline.checkDirectoryExists(directory['FIGURES'],'Figures directory',verbose=True,initialize=True)
print

print '#############################################'
print 'READING PARAMETERS FILES OR COMMAND ARGUMENTS'
print 
print '+ READING INTERVALS PARAMETER'
interval_dict = pipeline.initializeParameter(files,interval,'INTERVAL',sep=',',integerize=True)
print ' ',interval_dict
print
print '+ READING SUBSETTING REQUESTS'
subset_dict = pipeline.initializeParameter(files,subset,'SUBSET',sep=',',integerize=False)
print ' ',subset_dict
print
print '+ READING FLAGS'
flag_dict = pipeline.initializeParameter(files,flag,'FLAG',sep=',',integerize=False)
print ' ',flag_dict
print
print '+ READING HYPOTHESIS'
hypo_dict = pipeline.initializeParameter(files,hypothesis,'HYPO',sep='\+|,',integerize=False)
print ' ',hypo_dict 
print
print '+ READING META-DATA'
df_meta, df_meta_plates = pipeline.checkMetaTxt(files['META'],verbose=True)
print '#############################################'
print 'READING & ASSEMBLING DATA'
print
#if filename:
#    list_data = [os.path.basename(filename)]
#else:


data = pipeline.readPlateReaderFolder(folderpath=directory['DATA'],save=True,save_dirname=directory['DERIVED'],interval=600,interval_dict=interval_dict)
print 

# if plot_plate_only:
#     # check that each plate has 96 well IDs (in 8x12) format

#     for fname,table in data.iteritems():

#         print fname
#         print table.head()
#         actual_well_ids = set(sorted(table.columns.values[1:])); print actual_well_ids


#         # if a plate is Biolog, it seem to automatically have 96 wells internally even if rows are missing
#         expected_well_ids = set(sorted(plates.parseWellLayout().index.values)); print expected_well_ids

#         print expected_well_ids.difference(actual_well_ids)

#     sys.exit()

print '#############################################'
print 'READING & ASSEMBLING mapping'
print
list_data = sorted(os.listdir(directory['DATA']));
for filename in list_data:
    filebase = os.path.splitext(filename)[0]; 
    mapping_path = '%s/%s.txt' % (directory['MAPPING'],filebase); print mapping_path
    well_ids = data[filebase].columns[1:];
    mapping[filebase] = pipeline.smartmapping(filebase,mapping_path,well_ids,df_meta,df_meta_plates)
master_mapping = pd.concat(mapping.values(),ignore_index=True)#,sort=False)
print

master_mapping =pipeline. dropFlaggedWells(master_mapping,flag_dict,verbose=True);
#master_mapping.to_csv('%s/stitched_mapping.txt' % directory['mapping'],sep='\t',header=True,index=True)
print

print '#############################################'
print 'SUBSETTING MAPPING & DATA BASED ON USER INPUT'
subset_dict = pipeline.initializeParameter(files,subset,'SUBSET',sep=',',integerize=False)
print
master_mapping = pipeline.subsetWells(master_mapping,subset_dict,verbose=True).reset_index(drop=True)#; 
print 
master_mapping.index.name='Sample_ID';
master_mapping = master_mapping.reset_index(drop=False); 
print master_mapping

if master_mapping.shape[0] == 0:
    print('ERROR: no wells selected. System Exitting!')
    sys.exit()

wide_data_dict,key_dict = pipeline.subsetCombineData(data,master_mapping)
gplate = pipeline.packageGrowthPlate(wide_data_dict,key_dict)
print


print '#############################################'
print 'Running GP Regression'
print
gplate.computeFoldChange()
gplate.convertTimeUnits()
gplate.logData()
gplate.subtractBaseline()

#gplate.addRowColVarbs()

visual_check = False;#False;#plot_plate_only

if visual_check:

    gplate.addRowColVarbs()
    filepath = '%s/figure.pdf' % (directory['FIGURES'])
    fig,axes = gplate.plot(savefig=True,title="",filepath=filepath,modified=True)
    
    sys.exit('DONE')

if len(hypo_dict) > 0 :
    print hypo_dict

    start = time.time()
    gplate.runTestGP(hypothesis=hypo_dict,thinning=11,permute=True)
    print
    print('ELAPSED TIME:\t%s' % (time.time() - start))

    sys.exit('~~~DONE~~~')

#sys.exit()


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


