#!/usr/bin/env python

# Firas Said Midani
# Start date: 2018-10-08
# Final date: 2018-10-23	

# DESCRIPTION Library of functions for processing plate reader data at the Britton Lab

# TABLE OF CONTENTS
#
#
#|-- Direcotry Parsing
#    |-- breakDownFilePath
#    |-- createFolder
#
#|-- Text Parsing
#    |-- BOM_to_CSV
#    |-- check_BOM
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
#    |-- summarizeSugarData
#
#|-- DataFrame Initializing
#    |-- findSugarBiolog
#    |-- listTimePoints
#    |-- parseBiologLayout
#    |-- parseWellLayout
#
#|-- Data Summarizing
#    |-- summarizeGrowthData
#
#|-- Auxiliary
#    |-- getFormattedtime
#    |-- nRGB


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

sns.set_style('whitegrid');


# BEGIN FUNCTIONS

def createFolder(directory):
	'''
	Creates a folder only if it does not exist

    Keyword arguments:
    directory -- string

    Returns None
    '''

	if not os.path.exists(directory):
   		os.makedirs(directory)

   	return None

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

def findSugarBiolog(sugar):
	'''
	Identifies location of well with a particular sugar in Biolog PM plates (1-7)

	Keyword arguments:
	sugar -- string (must be an exact match with a sugar in biolog_pm_layout.py) 

	Returns an numpy.array of tuples
	'''

    bm = parseBiologLayout();
    hits = bm.where(bm == sugar).dropna(how='all')
    hits = hits.dropna(axis=1,how='all').T

    hits_list = []
    for idx, row in hits.iterrows():
        hits_list.append((idx[-1],(row==sugar).idxmax()))
        
    return hits_list

def getFormattedTime():
	'''
	Constructs time stamp formatted as Year-Month-Day-Hour_Minute_Second

	e.g. '2018-10-10-13-51-12'
	'''

	ts = time.localtime()
	ts = time.strftime("%Y-%m-%d-%H-%M-%S",ts)

	return ts

def nRGB(tup):
    '''
    Normalize RGB coordinates to values between 0 and 1

    Keyword arguments:
    tup -- tup with three values, where each value ranges between (and including) 0 and 255

    Returns tuple with three values, where each value ranges between (and including) 0 and 1
    '''
    
    return tuple([float(ii)/255 for ii in tup])
    
def parseBiologLayout():
	'''
	Initializes a pandas.DataFrame that maps the location (well and plate number) of each sugar in Biolog PM plates

	Returns pandas.DataFrame
	'''

	biolog_layout = pd.DataFrame([Carbon1,Carbon2,PhosphorusAndSulfur,PeptideNitrogen1,
                                 PeptideNitrogen2,PeptideNitrogen3],
                                 index=['PM1','PM2','PM3','PM4','PM5','PM^'],
                                 columns=parseWellLayout().index).T

	return biolog_layout

def parseWellLayout():
	'''
	Initializes a pandas.DataFrame where indices are well identifiers (e.g. C8)
	and columns indicate row number and column number

	Returns pandas.DataFrame
	'''

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
	ymax = np.ceil(df.max().max()); 
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
	    
	    ax.set_ylim([0,ymax])
	    ax.set_xlim([0,xmax])

	    ax.plot(df.columns,df.loc[idx,:],color=color_l,lw=1.5)
	    
	    ax.fill_between(x=df.columns,y1=[0]*df.shape[1],y2=df.loc[idx,:],color=color_f)

	    
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

	    ax.text(1., 1., "%0.2f" % summary.loc[idx,'Max OD'], color='black',
	            horizontalalignment='right', verticalalignment='top', 
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

	plt.close()

	return fig,axes

def subPlotSplit(df,nCols=4):
    ''' 
    With a limit of four columns in a sub-plot grid
    '''
    
    number = df.shape[0];
    
    numCols = nCols # number of columns in grid
    
    numRows = (number / numCols) + [1 if (number % numCols) else 0][0]; # number of needed rows in grid 
    
    df_subplot = pd.DataFrame(index=df.sort_index().index,columns=['PlotRow','PlotCol'])
    
    for ii in range(df.shape[0]):

        df_subplot.iloc[ii,:] = [ii/nCols,ii%nCols]
        
    df = df.join(df_subplot)
    
    return numRows,numCols, df

def plotPositivePlateGrowth(df_od,df_sugars,nCols=4,title="",savefig=0,filepath=""):
	'''



	'''

	# determine layout of grid
	nR,nC,df_sugars = subPlotSplit(df_sugars,nCols)

	fig,axes = plt.subplots(nR,nC,figsize=[nC+3,1.5*nR])

	df = df_od.loc[df_sugars.index]

	# round up window limits to integers
	ymax = np.ceil(df.max().max())
	xmax = float(df.columns[-1])
	xmax_h = int(np.ceil(float(df.columns[-1])/60/60))

	count = 1;

	for idx,row in df_sugars.iterrows():
	    
	    rr,cc = row.loc[['PlotRow','PlotCol']].values;

	    if nC==1:
	    	ax = axes[rr]
	    else:
	    	ax = axes[rr,cc]

	    color_l = (0.0,0.40,0.0,1.00)
	    color_f = (0.0,0.40,0.0,0.35)

	    ax.set_ylim([0,ymax])
	    ax.set_xlim([0,xmax])

	    ax.plot(df.columns,df.loc[idx,:],color=color_l,lw=1.5)
	    
	    ax.fill_between(x=df.columns,y1=[0]*df.shape[1],y2=df.loc[idx,:],color=color_f)

	    # show tick labels for bottom left subplot only
	    if (rr==nR-1 and cc==0):
	        plt.setp(ax,yticks=[0,ymax])
	        plt.setp(ax,xticks=[0,xmax],xticklabels=[0,xmax_h])
	        
	        ax.set_xlabel('Time (hours)')
	        ax.set_ylabel('OD (620 nm)')
	    else:
	        plt.setp(ax,yticks=[0,ymax],yticklabels=[])
	        plt.setp(ax,xticks=[0,xmax],xticklabels=[])
	        
	    
	    sub_title = '%02i. %s' % (count,df_sugars.loc[idx,'PM'])
	    
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

	return fig,axes

def listTimePoints(interval,numTimePoints):
	'''
	Constructs a numpy.array of a time series based on time interval length and number of time points

	Keyword arguments:
	interval -- int or float (latter perferred)
	numTimePoints -- int

	Returns numpy.array
	'''

	return np.arange(start=0,stop=interval*numTimePoints,step=interval)

def isASCII(data):
	'''
	Reference https://unicodebook.readthedocs.io/guess_encoding.html
	'''

    try:
        data.decode('ASCII')
    except UnicodeDecodeError:
        return False
    else:
        return True

def check_BOM(data):
	'''
	If a string starts with a BOM marker, return marker if it corresponds to one of several UTF encoding types. 

	Returns list of strings

	Reference https://unicodebook.readthedocs.io/guess_encoding.html
	'''

	BOMS = (
    	(BOM_UTF8, "UTF-8"),
    	(BOM_UTF32_BE, "UTF-32-BE"),
    	(BOM_UTF32_LE, "UTF-32-LE"),
    	(BOM_UTF16_BE, "UTF-16-BE"),
    	(BOM_UTF16_LE, "UTF-16-LE"),
	)

	return [encoding for bom, encoding in BOMS if data.startswith(bom)]

def BOM_to_CSV(filepath,newfile,encoding):
	'''
	For text files marked with BOM, convert to CSV format.

    Keyword arguments:
    filepath -- string
    newfile -- new filename (string)
    encoding -- string of encoding (e.g. 'utf-8')

    Returns newfile (stirng)
    '''

	csvReader = csv.reader(codecs.open(filepath, 'rU', encoding))

	fid = open(newfile,'w')
	for line in csvReader:
		fid.write('%s\n' % line[0].strip('\t'))  
	fid.close()

	return newfile

def readPlateReaderData(filepath,interval=600):

	filename, filebase, newfile = breakDownFilePath(filepath)

	content = open(filepath).readlines();
	sneakPeak = content[0];

	if isASCII(sneakPeak):
		
		print '%s is encoded with ASCII' % filename

	elif check_BOM(sneakPeak):
		
		encoding = check_BOM(content[0])[0]

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
	'''
	Breaks down a file path into several components (e.g. /Users/firasmidani/RandomFileName.asc)
	* filename (basename without path) --> RandomFileName.asc 
	* filebase (basename without path and extension) --> RandomFileName
	* newfile (filename with extension replaced to .tsv) --> /Users/firasmidani/RandomFileName.tsv)

	Returns list of three strings
	'''


	filename = os.path.basename(filepath);
	filebase = "".join(filename.split('.')[:-1]);
	dirname = os.path.dirname(filepath);
	newfile = '%s/%s.tsv' % (dirname,filebase);

	return filename, filebase, newfile

def summarizeGrowthData(df):
    '''
    summarizes the location and growth statistics for each well in a plate 

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

def summarizeSugarData(df_dict,sub_plate_list,nCols=4,title=sugar,savefig=0,filepath=""):

    nR, nC, df_plot = subPlotSplit(sub_plate_list,nCols=nCols);
    print nR, nC, df_plot

    fig,axes = plt.subplots(nR,nC,figsize=[nC+3,2.5*nR],sharey=True)   

    g_xmax,g_ymax = 0,0;

    for idx,row in df_plot.iterrows():

        df_sub_plot = row;
        rr,cc = df_sub_plot.loc[['PlotRow','PlotCol']]; #print rr, cc
        #plates = row.drop(['PlotRow','PlotCol']).values; #print plates
        plates = df_sub_plot.loc['Plates'];
        
        max_od = 0;

        for plate in plates:

            ax = axes[rr,cc]; #print plate

            od_ctrl = summarizeGrowthData(df_dict[plate]).loc['A1','Max OD']
            od_case = summarizeGrowthData(df_dict[plate]).loc[well,'Max OD']
            od_ratio = float(od_case)/od_ctrl        
            max_od += od_ratio

            x = df_dict[plate].columns; #print len(x)
            y_control = df_dict[plate].loc['A1'].values;
            y_case = df_dict[plate].loc[well].values;

            y_control = y_control - y_control[0]
            y_case = y_case - y_case[0]

            c_ctrl = nRGB((0,0,0));

            c_case =  [(0,0,1) if od_ratio>1.5 else (1,0,0)][0]

            ax.plot(x,y_control,color=c_ctrl,alpha=0.5,lw=3.5)
            ax.plot(x,y_case,color=c_case,alpha=0.5,lw=5)

            ax.fill_between(x=x,y1=[0]*len(y_control),y2=y_control,color=(0,0,0,0.1))

            ymax = np.ceil(np.max([y_control[:-1].max(),y_case[:-1].max()]))
            xmax = float(df_dict[plate].columns[-1]); #print xmax

            g_ymax = max(g_ymax,ymax)
            g_xmax = max(g_xmax,xmax)

            if (rr==nR-1 and cc==0):
                ax.set_xlabel('Time (hours)',fontsize=12,fontweight='bold')
                ax.set_ylabel('OD (620 nm)',fontsize=12,fontweight='bold')

            ax.set_title(row.name,fontsize=15,fontweight='bold')

            [ii.set(fontsize=12,fontweight='bold') for ii in ax.get_xticklabels()+ax.get_yticklabels()]

            plt.setp(ax,xticks=ax.get_xlim())


        # add Max OD to top right of each subplot
        max_od = max_od/len(plates);
        fontcolor = [(0,0,1) if max_od>1.5 else (1,0,0)][0]

        ax.text(1., 1., "%0.2f" % max_od, color=fontcolor,
            horizontalalignment='right', verticalalignment='top', 
            transform=ax.transAxes,fontsize=15)

    for ax in np.ravel(axes):
        
        xmax_h = int(np.ceil(float(g_xmax)/60/60))

        ax.set_ylim([0,g_ymax])
        ax.set_xlim([0,g_xmax])

        plt.setp(ax,yticks=[0,g_ymax])
        plt.setp(ax,xticks=[0,g_xmax],xticklabels=[0,xmax_h])

    l_row,l_col = df_plot.loc[idx,['PlotRow','PlotCol']].values
    [axes[l_row,col].axis('off') for col in range(l_col+1,nCols)];

    plt.suptitle(title,fontsize=20)

    plt.subplots_adjust(hspace=0.5,wspace=0.2)

    if savefig:
        if filepath=="":
            filepath = "/Users/firasmidani/Downloads/summarizeSugarData-%s.pdf" % getFormattedTime()
        plt.savefig(filepath,filetype='pdf')
            
    plt.close()

    return fig,axes