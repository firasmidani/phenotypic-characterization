#!/usr/bin/env python

# Firas Said Midani
# Start date: 2018-10-08
# Final date: 2018-10-15

# DESCRIPTION Library of functions for processing plate reader data at the Britton Lab

# TABLE OF CONTENTS
#
#
#|-- Data import


# IMPORT NECESSARY LIBRARIES

import os
import csv
import imp
import sys
import time
import codecs
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from codecs import BOM_UTF8, BOM_UTF16_BE, BOM_UTF16_LE, BOM_UTF32_BE, BOM_UTF32_LE


# IMPORT PERSONAL LIBRARIES

foo = imp.load_source('biolog_pm_layout','./biolog_pm_layout.py');
from biolog_pm_layout import *

# SET PARAMETERS & STYLES
#
#|-- Direcotry Parsing
#    |-- breakDownFilePath
#
#|-- Text Parsing
#    |-- BOM_to_CSV
#    |-- check_bom
#    |-- determineLineSkips
#    |-- isASCII
#
#|-- Data Processing
#    |-- readPlateReaderData
#
#|-- Plotting
#    |-- plotPlateGrowth
#    |-- plotPositivePlateGrowth
#    |-- subPlotSplit
#
#|-- DataFrame Initializing
#    |-- listTimePoints
#    |-- parseBiologLayout
#    |-- parseWellLayout
#
#|-- Data Summarizing
#    |-- summarizeGrowthData
#
#|-- Auxiliary
#    |-- getFormattedtime


sns.set_style('whitegrid');

def determineLineSkips(filepath):

	fid = open(filepath,'r');

	header_found = 0;
	count = 0;
	for line in fid.readlines():
	    if line.startswith('A1'):
	    	header_found=1;
	        break
	    count+=1;

	if header_found ==0:
		count=0;
	    
	fid.close()

	return count

def getFormattedTime():
	'''
	returns time stamp formatted as Year-Month-Day-Hour_Minute_Second

	e.g. '2018-10-10-13-51-12'
	'''

	ts = time.localtime()
	ts = time.strftime("%Y-%m-%d-%H-%M-%S",ts)

	return ts

def parseBiologLayout():

	biolog_layout = pd.DataFrame([Carbon1,Carbon2,PhosphorusAndSulfur,PeptideNitrogen1,
                                 PeptideNitrogen2,PeptideNitrogen3],
                                 index=['PM1','PM2','PM3','PM4','PM5','PM^'],
                                 columns=parseWellLayout().index).T

	return biolog_layout

def parseWellLayout():

	legend_row = {'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8}

	rows_list = np.ravel([[ii]*12 for ii in sorted(legend_row.keys())])
	cols_list = range(1,13)*8;

	df = pd.DataFrame([(xx,yy,'%s%s' % (xx,yy)) for xx,yy in zip(rows_list,cols_list)],
	                       columns=['Row','Col','Well'])
	df = df.set_index('Well')
	
	return df

def plotPlateGrowth(df,summary,threshold=1.5,title="",savefig=0,filepath=""):

	fig,axes = plt.subplots(8,12,figsize=[12,8])

	# round up window limits to integers
	ymax = np.ceil(df.max().max())
	xmax = float(df.columns[-1])
	xmax_h = int(np.ceil(float(df.columns[-1])/60/60))

	for idx in df.index:
	    
	    r,c = summary.loc[idx,['Row','Col']].values-1;

	    ax = axes[r,c]
	    
	    # green if above threshoold, gray if below
	    if summary.loc[idx,'Growth Fold']>threshold:
	        color_l = (0.0,0.40,0.0,1.00)
	        color_f = (0.0,0.40,0.0,0.35)
	    else:
	        color_l = (0.,0.,0.,1.00)
	        color_f = (0.,0.,0.,0.15)
	    
	    ax.plot(df.columns,df.loc[idx,:],color=color_l,lw=1.5)
	    
	    ax.fill_between(x=df.columns,y1=[0]*df.shape[1],y2=df.loc[idx,:],color=color_f)
	    
	    ax.set_ylim([0,ymax])
	    ax.set_xlim([0,xmax])
	    
	    # show tick labels for bottom left subplot only
	    if (r==7 and c==0):
	        plt.setp(ax,yticks=[0,ymax])
	        plt.setp(ax,xticks=[0,xmax],xticklabels=[0,xmax_h])
	    else:
	        plt.setp(ax,yticks=[0,ymax],yticklabels=[])
	        plt.setp(ax,xticks=[0,xmax],xticklabels=[])
	    
	    # add well identifier on top left of each subplot
	    ax.text(0., 1., idx, color=(0,0,1,0.5),
	            horizontalalignment='left', verticalalignment='top', 
	            transform=ax.transAxes)
	    
	fig.text(0.515, 0.07, 'Time (hours)', fontsize=15, 
	         ha='center', va='bottom', 
	         transform=ax.transAxes)
	fig.text(0.1, 0.5, 'Optical Density (620 nm)', fontsize=15, 
	         va='center', ha='right', rotation='vertical',
	         transform=ax.transAxes)

	fig.suptitle(title,fontsize=15)

	if savefig:

		if filepath=="":

			filepath = "/Users/firasmidani/Downloads/plotPlateGrowth-%s.pdf" % getFormattedTime()
	
		plt.savefig(filepath,filetype='pdf')

def subPlotSplit(df,nCols=4):
    ''' 
    With a limit of four columns in a sub-plot grid
    '''
    
    number = df.shape[0];
    
    numCols = nCols # number of columns in grid
    
    numRows = (number / numCols) + (number % numCols); # number of needed rows in grid 
    
    df_subplot = pd.DataFrame(index=df.index,columns=['PlotRow','PlotCol'])
    
    for ii in range(df.shape[0]):

        df_subplot.iloc[ii,:] = [ii/nCols,ii%nCols]
        
    df = df.join(df_subplot)
    
    return numRows,numCols, df

def plotPositivePlateGrowth(df_od,df_sugars,nCols=4,title="",savefig=0,filepath=""):
	'''



	'''

	# determine layout of grid
	nR,nC,df_sugars = subPlotSplit(df_sugars,nCols)

	print df_sugars, [2*nR,nC]

	fig,axes = plt.subplots(nR,nC,figsize=[nC+3,1.5*nR])

	df = df_od.loc[df_sugars.index]

	# round up window limits to integers
	ymax = np.ceil(df.max().max())
	xmax = float(df.columns[-1])
	xmax_h = int(np.ceil(float(df.columns[-1])/60/60))

	count = 1;

	for idx,row in df_sugars.iterrows():
	    
	    rr,cc = row.loc[['PlotRow','PlotCol']].values; print rr,cc

	    if nC==1:
	    	ax = axes[rr]
	    else:
	    	ax = axes[rr,cc]
	    print ax

	    color_l = (0.0,0.40,0.0,1.00)
	    color_f = (0.0,0.40,0.0,0.35)

	    ax.plot(df.columns,df.loc[idx,:],color=color_l,lw=1.5)
	    
	    ax.fill_between(x=df.columns,y1=[0]*df.shape[1],y2=df.loc[idx,:],color=color_f)

	    ax.set_ylim([0,ymax])
	    ax.set_xlim([0,xmax])

	    # show tick labels for bottom left subplot only
	    if (rr==nR-1 and cc==0):
	        plt.setp(ax,yticks=[0,ymax])
	        plt.setp(ax,xticks=[0,xmax],xticklabels=[0,xmax_h])
	        
	        ax.set_xlabel('Time (hours)')
	        ax.set_ylabel('OD (620 nm)')
	    else:
	        plt.setp(ax,yticks=[0,ymax],yticklabels=[])
	        plt.setp(ax,xticks=[0,xmax],xticklabels=[])
	        
	    sub_title = '%02i. %s' % (count,df_sugars.loc[idx,'PM1'])
	    
	    if nC==1:
	    	transform_ax = axes[rr];
	    else:
	    	transform_ax = axes[rr,-1];

	    plt.text(1.5, 1-float(cc)/nCols, sub_title, fontsize=13,
	             va='top',ha='left',transform=transform_ax.transAxes)

	    # add well identifier on top left of each subplot
	    ax.text(1., 1.,'%02i' % count, color=(0,0,0,1.),
	            horizontalalignment='right', verticalalignment='top', 
	            transform=ax.transAxes)

	    count+=1;
	    
	# turn off axes for remaining unused subplots
	[axes[rr,col].axis('off') for col in range(cc+1,nCols)];

	if nC==1:
		transform_ax = axes[0];
	else:
		transform_ax = axes[0,0];

	ax.text(0.,1.4,title,fontsize=15,transform=transform_ax.transAxes)

	if savefig:

		if filepath=="":

			filepath = "/Users/firasmidani/Downloads/plotPositivePlateGrowth-%s.pdf" % getFormattedTime()
	
		plt.subplots_adjust(right=0.5)
		plt.savefig(filepath,filetype='pdf')

def listTimePoints(interval,numTimePoints):

	return np.arange(start=0,stop=interval*numTimePoints,step=interval)

def isASCII(data):
    try:
        data.decode('ASCII')
    except UnicodeDecodeError:
        return False
    else:
        return True

def check_bom(data):

	BOMS = (
    	(BOM_UTF8, "UTF-8"),
    	(BOM_UTF32_BE, "UTF-32-BE"),
    	(BOM_UTF32_LE, "UTF-32-LE"),
    	(BOM_UTF16_BE, "UTF-16-BE"),
    	(BOM_UTF16_LE, "UTF-16-LE"),
	)

	return [encoding for bom, encoding in BOMS if data.startswith(bom)]

def BOM_to_CSV(filepath,newfile,encoding):

	csvReader = csv.reader(codecs.open(filepath, 'rU', encoding))

	fid = open(newfile,'w')
	for line in csvReader:
		fid.write('%s\n' % line[0].strip('\t'))  
	fid.close()

	return newfile

def readPlateReaderData(filepath,interval=6000):

	filename, filebase, newfile = breakDownFilePath(filepath)

	content = open(filepath).readlines();
	sneakPeak = content[0];

	if isASCII(sneakPeak):
		
		print '%s is encoded with ASCII' % filename

	elif check_bom(sneakPeak):
		
		encoding = check_bom(content[0])[0]

		filepath = BOM_to_CSV(filepath,newfile,encoding[0:6])

		print '%s is encoded with %s ' % (filename,encoding)

	else:

		print 'Parsing Error: Encoding for %s is unknown.'

		return

	skiprows = determineLineSkips(filepath); #print skiprows

	df = pd.read_csv(filepath,sep='\t',header=None,index_col=0,skiprows=skiprows);

	df.columns = listTimePoints(interval,df.shape[1])

	df.index.name = 'Well'
	df.T.index.name = 'Time'

	df.to_csv(newfile, sep='\t')

	return df

def breakDownFilePath(filepath):

	filename = os.path.basename(filepath);
	filebase = "".join(filename.split('.')[:-1]);
	dirname = os.path.dirname(filepath);
	newfile = '%s/%s.tsv' % (dirname,filebase);

	return filename, filebase, newfile

def summarizeGrowthData(df):
    '''
    summarizes the locationa and growth statistics for each well in a plate 

    Keyword arguments:
    df -- pandas.DataFrame (well x time)

    Returns pandas.DataFrame
    '''

    # initialize dataframe
    legend = pd.DataFrame(index=df.index,columns=['Row','Col','Let'])
    
    # map row letters to numbers
    legend_row = {'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8}

    # subtract T0 from other time points in each well
    df = df.apply(lambda col: col - df.loc[:,0], axis=0)

    # caculate fold growth as maximum OD in each well relative to negative control
    df_max = df.max(1)
    df_fold = df_max/df_max.loc['A1']

    	# add metadata to each well: well row as number, well col as number, well identifier
    for idx in df.index:
   	    legend.loc[idx,:] = [legend_row[idx[0]],int(idx[1:]),idx[0]]

    # add fold growth data to summary_table
    summary = legend.join(pd.DataFrame(df_max,columns=['Max OD']))
    summary = summary.join(pd.DataFrame(df_fold,columns=['Growth Fold']))

    return summary

