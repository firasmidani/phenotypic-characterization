#!/usr/bin/env python

# Firas Said Midani
# Start date: 2018-10-08
# Final date: 2018-10-08

# DESCRIPTION Library of functions for processing plate reader data at the Britton Lab

# TABLE OF CONTENTS
#
#
#|-- Data import


# IMPORT NECESSARY LIBRARIES

import os
import csv
import sys
import codecs
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# SET PARAMETERS & STYLES

sns.set_style('whitegrid');

def DetermineLineSkips(filepath):

	fid = open(filepath,'r');

	count = 0;
	for line in fid.readlines():
	    if line.startswith('Time'):
	        break
	    count+=1;
	    
	fid.close()

	return count


def readPlateReaderData(filepath,machine):

	filename = "".join(os.path.basename(filepath).split('.')[:-1])
	dirname = os.path.dirname(filepath)

	newfile = '%s/%s.tsv' % (dirname,filename)

	if machine=="Ackuskan":

		skiprows = DetermineLineSkips(filepath)

		df = pd.read_csv(filepath,sep='\t',header=0,index_col=0,skiprows=skiprows);

		# convert column headers to int
		df.columns = [float(ii) if not isinstance(ii,float) else ii for ii in df.columns ]

		df.index.name = 'Well'
		df.T.index.name = 'Time'

		df.to_csv(newfile, sep='\t')

		return df

	elif machine=="Magellan":


		csvReader = csv.reader(codecs.open(filepath, 'rU', 'utf-16'))

		fid = open(newfile,'w')

		for line in csvReader:
			fid.write('%s\n' % line[0].strip('\t'))

		fid.close()

		df = pd.read_csv(newfile,sep='\t',skiprows=1)

		# convert column headers to int
		df.columns = [float(ii[:-1]) if not isinstance(ii,float) else ii for ii in df.columns]

		df.index.name = 'Well'
		df.T.index.name = 'Time'

		df.to_csv(newfile, sep='\t')

		return df 

