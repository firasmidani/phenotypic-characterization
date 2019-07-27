#!/usr/bin/env python

# Firas Said Midani
# Start date: 2018-10-08
# Final date: 2019-05-06

# DESCRIPTION Library of functions for processing plate reader data at the Britton Lab

# TABLE OF CONTENTS (functions are scripted in alphabetical order)
#
#
#|-- Direcotry Parsing
#    |-- breakDownFilePath
#    |-- findPlateReaderFiles
#    |-- createFolder
#
#|-- Text Parsing
#    |-- BOM_to_CSV
#    |-- check_BOM
#    |-- findFirstRow
#    |-- isASCII
#    |-- parsePlateName
#
#|-- Data Processing
#    |-- readPlateReaderData
#    |-- smoothGrowthCurves
#
#|-- Plotting
#    |-- definePlotLayout
#    |-- modifyTickLabels
#    |-- plotGroupedGrowth
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
#    |-- populatePlateKey
#
#|-- Data Summarizing
#    |-- summarizeGrowthData
#
#|-- Auxiliary
#    |-- getFormattedtime
#    |-- nRGB

# TO DO

# 1. resolve differences between summarizeGrowthData and summarizeGrowthDataModified
# 2. stop at H12 for readPLateReader
# IMPORT NECESSARY LIBRARIES

import os
import csv
import imp
import sys
import time
import codecs
import string
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from codecs import BOM_UTF8, BOM_UTF16_BE, BOM_UTF16_LE, BOM_UTF32_BE, BOM_UTF32_LE

from scipy.signal import savgol_filter

# IMPORT PERSONAL LIBRARIES

sys.path.append('..')

from config import biolog_pm_layout as bpl 

# SET PARAMETERS & STYLES

sns.set_style('whitegrid');

# BEGIN FUNCTIONS

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

def breakDownFilePath(filepath,save_dirname=None):
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

    if save_dirname:
        newfilepath = '%s/%s.tsv' % (save_dirname,filebase);
    else:
        newfilepath = '%s/%s.tsv' % (dirname,filebase);

    return filename, filebase, newfilepath

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

def definePlotLayout(ISOLATE_MAPPING,ROW_MAPPING_DICT,ROW_COLOR_DICT):
    '''
    definePLotLayout determines the rows on which each isolate will be plotted (ROW MAPPING). 
    Similar to subPlotSplit but better!

    Keyword Arguments:


    ISOLATE_MAPPING -- pandas.DataFrmae where each row is an isolate, 
                       column 'Plates' has values that are LIST of biolog plate IDS.
    ROW_MAPPING_DICT -- dictionary where keys are INT indicators for rows and values are lists of isolate names
    ROW_COLOR_DICT -- dictionary where keys are INT indicators for rows and values are STRING for colors

    Returns pandas.DataFrame

    ISOLATE_MAPPING -- ISOALTE_MAPPING but with three additional columns: 
                       PlotRow, PlotCol columns store INT indicators of sub-plot location, 
                       GroupColor column storess STRING indicator of title (?) color
    '''
    
    ROWS_MAPPING_LIST = [];

    ## map each key and value in dict into p x 2 list
    for key,value_list in ROW_MAPPING_DICT.iteritems():
        for value in value_list:
            ROWS_MAPPING_LIST.append([value,key])
            
    ## convert ROW MAPPING to a pandas.DataFrmae
    PLOT_MAPPING_DF = pd.DataFrame(ROWS_MAPPING_LIST,columns=['Isolate','PlotRow'])
    PLOT_MAPPING_DF.set_index(['Isolate'],inplace=True)
    PLOT_MAPPING_DF

    ## add COL MAPPING based on alphabetic ordering of isolate ID
    for row,isolates in PLOT_MAPPING_DF.groupby('PlotRow').groups.iteritems():
        PLOT_MAPPING_DF.loc[sorted(isolates),'PlotCol'] = [int(ii) for ii in range(len(isolates))]
        
    ## add ROW MAPPING to ISOLATE MAPPING
    ISOLATE_MAPPING = ISOLATE_MAPPING.join(PLOT_MAPPING_DF.astype(int))
    
    for row,isolates in ISOLATE_MAPPING.groupby('PlotRow').groups.iteritems():
        ISOLATE_MAPPING.loc[isolates,'GroupColor'] = ROW_COLOR_DICT[row]
        
    return ISOLATE_MAPPING

def findFirstRow(filepath):
    '''
    searches for the line beginning with well ID (A1), 
    determines the number of rows to skip for reading this line,
    and indicates in return if index column was not found (no A1 found)
    '''

    fid = open(filepath,'r');

    count = 0;
    for line in fid.readlines():
        if line.startswith('A1'):
            fid.close()
            index_column = 0
            return count,index_column
        count += 1;
    else:
        fid.close()
        count = 0
        index_column = None
        return count, index_column

def findPlateReaderFiles(directory):
    '''Given a path to a directory, find all files nested in the folder
     that correspond to plate reader data in .TXT or .asc ofrmats'''

    ls_files = [];

    for (dirpath, dirnames, filenames) in os.walk(directory):
        for filename in filenames:
            if filename.endswith(".TXT") or  filename.endswith(".asc"):
                ls_files.append('%s/%s' % (dirpath,filename))
    
    return ls_files

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

def listTimePoints(interval,numTimePoints):
    '''
    Constructs a numpy.array of a time series based on time interval length and number of time points

    Keyword arguments:
    interval -- int or float (latter perferred)
    numTimePoints -- int

    Returns numpy.array
    '''

    return np.arange(start=0,stop=interval*numTimePoints,step=interval)

def modifyTickLabels(axis,fontsize=15,fontweight='normal',which='both'):

    if which=='both':

        [ii.set(fontsize=fontsize,fontweight=fontweight) for ii in ax.get_xticklabels()+ax.get_yticklabels()];

    elif which=='x':

        [ii.set(fontsize=fontsize,fontweight=fontweight) for ii in ax.get_xticklabels()];

    elif which=='y':

        [ii.set(fontsize=fontsize,fontweight=fontweight) for ii in ax.get_yticklabels()];

    return axis

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

    biolog_layout = pd.DataFrame([bpl.Carbon1,bpl.Carbon2,bpl.PhosphorusAndSulfur,
                                  bpl.PeptideNitrogen1,bpl.PeptideNitrogen2,bpl.PeptideNitrogen3],
                                 index=['PM1','PM2','PM3','PM4','PM5','PM^'],
                                 columns=parseWellLayout(order_axis=0).index).T

    return biolog_layout

def parsePlateName(plate):
    
    isolate = str(plate.split('_')[0]);
    pm = str(plate.split('PM')[1][0]);
    
    return isolate,pm

def parseWellLayout(order_axis=1):
    '''
    Initializes a pandas.DataFrame where indices are well identifiers (e.g. C8)
    and variables indicate row letter and column number

    order_axis = 0 indicates order IDs by column (i.e. A1,B1,C1, ... , A2,B2,C2, ...)
    order_axis = 1 indicates order IDs by column (i.e. A1,A2,A3, ... , B1,B2,B3, ...)

    Returns pandas.DataFrame
    '''

    rows = list(string.ascii_uppercase[0:8]);
    cols = [str(ii) for ii in range(1,13)];

    list_wells = [];
    list_cols = [];
    list_rows = [];

    if order_axis==1:
        for col in cols:
            for row in rows:
                list_wells.append('{}{}'.format(row,col))
                list_rows.append(row)
                list_cols.append(col)
    else:
        for row in rows:
            for col in cols:
                list_wells.append(('{}{}'.format(row,col)))
                list_rows.append(row)
                list_cols.append(col)

    df = pd.DataFrame([list_wells,list_rows,list_cols],index=['Well','Row','Column'])
    df = df.T
    df = df.set_index('Well')

    return df

def plotGroupedGrowth(subtrate,df_growth_means,df_plate,df_dict):
    '''
    plotGroupedGrowth returns a figure with multiple panels. Each panel is a 
    ribotype with multiple lines corresponding to growth of distinct isolates.

    Keyword arguments:

    substrate -- list of strings, e.g. ['D-Trehalose]
    df_growth_means -- pandas.DataFrame with substrates as indices and isolates as columns;
                       needed to extract fold change vlaues
    df_plate -- pandas.DataFrame where each row corresponds to a unique Biolog run.
                Columns include Isoalte, PM (biolog number), Rep (replicate number), and Ribotype
    df_dict -- dictionary of pandas.DataFrames, where each value is the formatted Biolog output (well x time)

    Returns:

    fig --- matplotlib.figure.Figure
    ax -- matplotlib.axes._subplots.AxesSubplot
    '''

    pid, wid = findSugarBiolog(subtrate[0])[0];

    color_list = ['purple','blue','red','orange','gold','green','brown']
    #color_list = ['purple','green','brown']

    # isolate as index, ribotype as column varibale
    map_ribotype = df_plate.reset_index().loc[:,['Isolate','Ribotype']]
    map_ribotype = map_ribotype.drop_duplicates().set_index(['Isolate']);

    # isolate as index, fold-change and ribotype as column variables
    df_carbon = df_growth_means.loc[subtrate[0],:].T

    #print pd.DataFrame(df_carbon)

    df_carbon = pd.DataFrame(df_carbon).join(map_ribotype).dropna()

    # group isoaltes by ribotype
    df_carbon_g = df_carbon.groupby(['Ribotype']).groups

    # ribotypes as index, isolates, rows, and columns as column variables
    df_carbon_plot = pd.DataFrame(columns=['Isolates']);
    for key,value in df_carbon_g.iteritems():
        df_carbon_plot.loc[key,'Isolates'] = value.values 

    # initialize parameters for sub-plotting (i.e. ribotyes mapping to row and col)
    nR,nC,df_carbon_plot = subPlotSplit(df_carbon_plot,nCols=3)

    fig,axes = plt.subplots(nR,nC,figsize=[nC+7.5,5*nR],sharex=False,sharey=True)

    # manually define xticks and xticklabels
    xticks  = [float(ii)*3600 for ii in [0,5,10,15,17]]
    xticklabels = [0,'','','',17];

    for group, row in df_carbon_plot.sort_index().iterrows():
        
        rr,cc = row.loc[['PlotRow','PlotCol']].values;
        
        if nR>1:
            ax = axes[rr,cc];
        else:
            ax = axes[cc]

        members = row.loc['Isolates'];

        #print members
        
        for cc,member in enumerate(members):

            plates = df_plate[df_plate.isin({'Isolate':[member],'PM':[int(pid)]}).sum(1)==2].index

            # parse through techical replicates
            for plate in plates:

                # subtract no-carbon control growth
                x = df_dict[plate].columns;
                y_neg = df_dict[plate].loc['A1',:];
                y = df_dict[plate].loc[wid,:];
                y = y - y_neg;

                ax.plot(x,y,lw=5,color=color_list[cc],alpha=0.65,label=member)
                
        # subplot aesthetics
        [ii.set(fontsize=20) for ii in ax.get_xticklabels()+ax.get_yticklabels()]
        plt.setp(ax,yticks=[0,0.25,0.50,0.75,1.0],yticklabels=[0,'','','',1.0])
        plt.setp(ax,xticks=xticks,xticklabels=[])

        ax.set_ylim([0,1.0])
        ax.set_xlim([0,xticks[-1]])

        ax.legend(fontsize=15)

        ax.set_title('RT%s' % ''.join(group.split('RT')[1:]),fontsize=15)

    # figure asesthetics
    if nR>1:
        axes[-1,0].set_xlabel('Time (hours)',fontsize=15);
        axes[-1,0].set_ylabel('OD (620 nm)',fontsize=15);
    else:
        axes[0].set_xlabel('Time (hours)',fontsize=20);
        axes[0].set_ylabel('OD (620 nm)',fontsize=20);
        
    # remove spines for unused subplots
    to_plot = list(zip(df_carbon_plot.PlotRow, df_carbon_plot.PlotCol));
    to_kill = list(set(itertools.product(range(nR),range(nC))).difference(to_plot));

    if nR>1:
        [axes[kill_row,kill_col].axis('off') for kill_row,kill_col in to_kill];
    else:
        [axes[kill_col].axis('off') for kill_row,kill_col in to_kill];

    # only include xticklabels for subplots on the bottom boundary

    ultimate_row = df_carbon_plot.PlotRow.max();
    penultimate_row = ultimate_row-1;
    penultimate_col = df_carbon_plot[df_carbon_plot.isin({'PlotRow':[ultimate_row]}).any(1)];
    penultimate_col = penultimate_col.PlotCol.max();

    for col in range(penultimate_col+1,nC):

        if nR>1:
            ax = axes[penultimate_row,col];
        else:
            ax = axes[col];
            
        plt.setp(ax,xticks=xticks,xticklabels=xticklabels);

    for col in range(0,penultimate_col+1):

        if nR>1:
            ax = axes[ultimate_row,col];
        else:
            ax = axes[col];

        plt.setp(ax,xticks=xticks,xticklabels=xticklabels);

    plt.subplots_adjust(hspace=0.3,wspace=0.3)

    fig.suptitle(subtrate[0],fontsize=30);

    return fig,ax

def plotPlateGrowth(df,summary,threshold=1.5,title="",savefig=0,filepath="",logged=False):

    fig,axes = plt.subplots(8,12,figsize=[12,8])

    # subtract T0 from other time points in each well
    #df = df.apply(lambda col: col - df.loc[:,0], axis=0)

    # round up window limits to integers
    ymax = np.ceil(df.max().max()); 
    ymin = np.floor(df.min().min());
    xmax = float(df.columns[-1])
    xmax_h = int(np.ceil(float(df.columns[-1])/60/60))

    for idx in df.index:
        
        r,c = summary.loc[idx,['Row','Column']].values-1;

        ax = axes[r,c]
        
        # green if above threshoold, gray if below
        if summary.loc[idx,'Fold Change']>threshold:
            #color_l = (0.0,0.40,0.0,1.00) # green
            #color_f = (0.0,0.40,0.0,0.35)
            color_l = (0.0,0.00,1.0,1.00) # blue
            color_f = (0.0,0.00,1.0,0.15)
        elif summary.loc[idx,'Fold Change']<0.50:
            #color_l = (1.0,0.6,0.0,1.00) # orange
            #color_f = (1.0,0.6,0.0,0.35)
            color_l = (1.0,0.0,0.0,1.00) # red
            color_f = (1.0,0.0,0.0,0.15)
        else:
            color_l = (0.,0.,0.,1.00) # black
            color_f = (0.,0.,0.,0.15)

        ax.set_ylim([0,ymax]);
        ax.set_xlim([0,xmax]);

        if logged:
            ax.set_ylim([ymin,ymax]);

        x = df.columns;
        y = df.loc[idx,:];

        ax.plot(x,y,color=color_l,lw=1.5);
        

        if logged:
            ax.fill_between(x=x,y1=[ax.get_ylim()[0]]*df.shape[1],y2=y,color=color_f);
        else:
            ax.fill_between(x=x,y1=[0]*df.shape[1],y2=y,color=color_f);

        # show tick labels for bottom left subplot only
        if (r==7 and c==0) and (logged):
            plt.setp(ax,yticks=[ymin,ymax]);
            plt.setp(ax,xticks=[0,xmax],xticklabels=[0,xmax_h]);
        elif (r==7 and c==0):
            plt.setp(ax,yticks=[0,ymax]);
            plt.setp(ax,xticks=[0,xmax],xticklabels=[0,xmax_h]);
        elif logged:
            plt.setp(ax,yticks=[ymin,ymax],yticklabels=[]);
            plt.setp(ax,xticks=[0,xmax],xticklabels=[]);
        else:
            plt.setp(ax,yticks=[0,ymax],yticklabels=[]);
            plt.setp(ax,xticks=[0,xmax],xticklabels=[]);

        # add well identifier on top left of each subplot
        well_color = (0.65,0.165,0.165,0.8);#(0,0,1,0.5)
        ax.text(0., 1., idx, color=well_color,
                horizontalalignment='left', verticalalignment='top', 
                transform=ax.transAxes);

        ax.text(1., 1., "%0.2f" % summary.loc[idx,'Max OD'], color='black',
                horizontalalignment='right', verticalalignment='top', 
                transform=ax.transAxes);
    
    # if logged:
    #   ax.set_yscale('log')

    fig.text(0.515, 0.07, 'Time (hours)', fontsize=15, 
             ha='center', va='bottom', 
             transform=ax.transAxes);

    
    if logged:
        ylabel_text = 'log(Optical Density) (620 nm)';
    else:
        ylabel_text = 'Optical Density (620 nm)';
            
    fig.text(0.1, 0.5, ylabel_text, fontsize=15, 
             va='center', ha='right', rotation='vertical',
             transform=ax.transAxes);

    fig.suptitle(title,fontsize=15);

    if savefig:

        if filepath=="":

            filepath = "/Users/firasmidani/Downloads/plotPlateGrowth-%s.pdf" % getFormattedTime()
    
        plt.savefig(filepath,filetype='pdf')

    plt.close()

    return fig,axes

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

def initializePlateKey(plate):
    
    isolate,pm = parsePlateName(plate);
    
    pm = 'PM%s' % pm
    
    biolog = parseBiologLayout().loc[:,pm];
    
    wells = biolog.index;
    
    substrate = biolog.values;
    
    isolate = [isolate]*len(substrate);
    
    plate = [plate]*len(substrate);
    
    key = pd.DataFrame([wells,plate,isolate,substrate],
                      index=['Well','Plate','Isolate','Substrate']);
    
    key = key.T

    key = key.set_index(['Well'])
    
    return key

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

def readPlateReaderData(filepath,interval=600,save=False,save_dirname=None):

    '''
    if save_dirname is None, the new files will be saved inplace (i.e. in same directory as original) 
    
    assumptions:

    * tab-delimited files

    '''

    filename, filebase, newfile = breakDownFilePath(filepath,save_dirname)

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

    skiprows,index_col = findFirstRow(filepath); #print skiprows

    df = pd.read_csv(filepath,sep='\t',header=None,index_col=index_col,skiprows=skiprows);

    df.columns = listTimePoints(interval,df.shape[1])

    if not index_col:
        df.index = parseWellLayout(order_axis=0).index.values
#    else: 
#        df = df.loc[parseWellLayout(order_axis=1).index.values,:]

    df.index.name = 'Well'
    df.T.index.name = 'Time'

    # this makes sure to grab only rows that begin with well ids
    df = df.loc[parseWellLayout(order_axis=1).index]

    if save:
        df.to_csv(newfile, sep='\t')

    return df

def smoothGrowthCurves(data,window,polyorder):

    if isinstance(data,pd.DataFrame):

        for idx,row in data.iterrows():

            data.loc[idx,:] = savgol_filter(row.values,window,polyorder)

        return data

    elif isinstance(data,list):

        return savgol_filter(data,window,polyorder)

    else:

        raise ValueError('data should be either a pandas.DataFrame or list')

def summarizeGrowthDataModified(df,subtract=1,smooth=1,smooth_args=(19,3)):
    '''
    summarizes the location and growth statistics for each well in a plate 

    Keyword arguments:
    df -- pandas.DataFrame (well x time)

    Returns pandas.DataFrame
    '''

    # initialize dataframe
    legend = pd.DataFrame(index=df.index,columns=['Well','Row','Column','Letter'])
    
    # map row letters to numbers
    legend_row = {'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8}

    df_baseline = df.loc[:,0]; 

    # subtract T0 from other time points in each well
    if subtract:
        df = df.apply(lambda col: col - df.loc[:,0], axis=0)

    if smooth:
        df = smoothGrowthCurves(df,smooth_args[0],smooth_args[1])

    # caculate fold growth as maximum OD in each well relative to negative control
    df_max = df.max(1); 
    df_fold = (df_max/df_max.loc['A1']); 

    # add metadata to each well: well row as number, well col as number, well identifier
    for idx in df.index:
        legend.loc[idx,:] = [idx,legend_row[idx[0]],int(idx[1:]),idx[0]]

    summary = pd.DataFrame(index=['Max OD','Growth Fold','Baseline'],
                           columns=df.index,
                                   data=[df_max.values,df_fold.values,df_baseline.values]).T

    summary = legend.join(summary)


    return summary

def summarizeGrowthData(df,subtract=1,smooth=1,smooth_args=(19,3),expand_well_id=True):
    '''
    summarizes the location and growth statistics for each well in a plate 

    Keyword arguments:
    df -- pandas.DataFrame (well x time)

    Returns pandas.DataFrame
    '''

    # initialize dataframe
    #legend = pd.DataFrame(index=df.index,columns=['Well','Row','Column','Letter'])
    legend = pd.DataFrame(index=df.index,columns=['Row','Column','Letter','Well_ID'])
    
    # map row letters to numbers
    legend_row = {'A':1,'B':2,'C':3,'D':4,'E':5,'F':6,'G':7,'H':8}

    df_baseline = df.loc[:,0]; 

    # subtract T0 from other time points in each well
    if subtract:
        df = df.apply(lambda col: col - df.loc[:,0], axis=0)

    if smooth:
        df = smoothGrowthCurves(df,smooth_args[0],smooth_args[1])

    # caculate fold growth as maximum OD in each well relative to negative control
    df_max = df.max(1);     
    df_fold = df_max/df_max.loc['A1']; 

    summary = pd.DataFrame(index=['Baseline OD','Max OD','Fold Change'],
                           columns=df.index,
                           data=[df_baseline.values,df_max.values,df_fold.values]).T
    
    if expand_well_id:
        
        # add metadata to each well: well row as number, well col as number, well identifier
        for idx in df.index:
            row = legend_row[idx[0]];
            col = int(idx[1:]);
            letter = idx[0];
            well_id = '%s%s' % (letter,col);

            legend.loc[idx,:] = [row,col,letter,well_id];

        summary = legend.join(summary)

    return summary

def summarizeSugarData(df_dict,sub_plate_df,sugar,subtract=0,nCols=6,title="",savefig=0,filepath=""):

    # find well location of sugar
    plate, well = findSugarBiolog(sugar)[0]; 

    # initialize parameters for figure panels
    #nR, nC, df_plot = subPlotSplit(sub_plate_df,nCols=nCols);
    nR = sub_plate_df.PlotRow.max()+1; 
    nC = sub_plate_df.PlotCol.max()+1; 
    df_plot = sub_plate_df;

    # columns for df_Plot are Plates, PlotRow, PlotCol, index is isolate ID
    #print nR, nC, df_plot

    # initialize multi-panel figure 
    #print nR,nC
    fig,axes = plt.subplots(int(nR),int(nC),figsize=[nC+3,2.5*nR],sharey=True)   

    # initialize global x-axis and y-axis window limits
    g_xmax,g_ymax = 0,0;

    # loop through isolates
    for idx,row in df_plot.iterrows():

        # extract info on location and data locations
        df_sub_plot = row; 
        rr,cc = df_sub_plot.loc[['PlotRow','PlotCol']]; #print rr, cc
        plates = df_sub_plot.loc['Plates'];
        
        # initialize od ratio (OR)
        mean_or = 0;

        # loop through each replicate/plate
        for plate in plates:

            ax = axes[rr,cc]; #print plate

            # 
            od_ctrl = summarizeGrowthData(df_dict[plate]).loc['A1','Max OD']
            od_case = summarizeGrowthData(df_dict[plate]).loc[well,'Max OD']
            od_ratio = float(od_case)/od_ctrl        
            mean_or += od_ratio

            x = df_dict[plate].columns; #print len(x)
            y_control = df_dict[plate].loc['A1'].values;
            y_case = df_dict[plate].loc[well].values;

            y_control = y_control - y_control[0]
            y_case = y_case - y_case[0]

            c_ctrl = nRGB((0,0,0));

            #c_case = [(0,0,1) if od_ratio>1.5 else (1,0,0)][0]
            #c_case = [(0,0.4,0) if od_ratio>1.5 else (1,0.6,0)][0]

            if od_ratio>1.5:

                c_case = (0,0,1);

            elif od_ratio<0.5:

                c_case = (1,0,0);

            else:
                c_case = (1,0.6,0);
                c_case = (1,0,0);

            if subtract==1:

                ax.plot(x,y_case-y_control,color=c_case,alpha=0.5,lw=5);

            elif subtract==0:
                
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

            ax.set_title(row.name,fontsize=15,fontweight='bold',color=df_sub_plot.loc['GroupColor'])

            [ii.set(fontsize=12,fontweight='bold') for ii in ax.get_xticklabels()+ax.get_yticklabels()]

            plt.setp(ax,xticks=ax.get_xlim())


        # add Mean OD to top right of each subplot
        mean_or = float(mean_or)/len(plates);
        fontcolor = [(0,0,1) if mean_or>1.5 else (1,0,0)][0]
        fontcolor = [(0,0.4,1) if mean_or>1.5 else (1,0.6,0)][0]

        if mean_or>1.5:

            fontcolor = (0,0,1);

        elif mean_or<0.5:

            fontcolor = (1,0,0);

        else:

            fontcolor = (1,0.6,0);
            fontcolor = (1,0,0);


        ax.text(1., 1., "%0.2f" % mean_or, color=fontcolor,
            horizontalalignment='right', verticalalignment='top', 
            transform=ax.transAxes,fontsize=15)

    for ax in np.ravel(axes):
        
        xmax_h = int(np.ceil(float(g_xmax)/60/60))

        ax.set_ylim([0,g_ymax])
        ax.set_xlim([0,g_xmax])

        plt.setp(ax,yticks=[0,g_ymax])
        plt.setp(ax,xticks=[0,g_xmax],xticklabels=[0,xmax_h])

    #l_row,l_col = df_plot.loc[idx,['PlotRow','PlotCol']].values
    #[axes[l_row,col].axis('off') for col in range(l_col+1,nCols)];

    to_plot = list(zip(df_plot.PlotRow, df_plot.PlotCol));
    to_kill = list(set(itertools.product(range(nR),range(nC))).difference(to_plot));

    [axes[kill_row,kill_col].axis('off') for kill_row,kill_col in to_kill]

    if title=="":
        plt.suptitle(sugar,fontsize=20)
    else:
        plt.suptitle(title,fontsize=20)

    plt.subplots_adjust(hspace=0.5,wspace=0.2)

    if savefig:
        if filepath=="":
            filepath = "/Users/firasmidani/Downloads/summarizeSugarData-%s.pdf" % getFormattedTime()
        plt.savefig(filepath,filetype='pdf')
            
    plt.close()

    return fig,axes
