#!/usr/env/bin python

# Firas Said Midani
# Start date: 2018-06-12
# Final date: 2019-06-12

# DESCRIPTION Library of functions for input/output of data

# TABLE OF CONTENTS (functions are scripted in alphabetical order)

#|-- Direcotry Parsing
#    |-- breakDownFilePath
#    |-- createFolder
#
#|-- Text Parsing
#    |-- BOM2CSV
#    |-- determineLineSkips
#    |-- isASCII
#    |-- parsePlateName
#    |-- whichBOM
#
#|-- Data Processing
#    |-- readPlateReaderData

# TO DO

# IMPORT NECESSARY LIBRARIES

def BOM2CSV(filepath,newfile,encoding):
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
	'''
	determineLineSkips searches text file for row starting with 'A1' value and returns row number

	Returns integer
	'''

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


def isASCII(data):
	'''
	Checks if a file is encoded with ASCII.

	Reference https://unicodebook.readthedocs.io/guess_encoding.html

	Returns a boolean
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

def readPlateReaderData(filepath,interval=600):
	'''
	Read plate reader output in text files encoded in ASCII or BOM (e.g. UTF-8),
	ignore meta-data stored as headers, and save only OD data.

	Args: 

		filepath (string): file with full or relatie path
		interval (int): time interval between consecutive measurements in seconds (britton lab default is 600)

	Returns pandas.DataFrame where each row is a well and each column is a time-specific measurement
	'''
	filename, filebase, newfile = breakDownFilePath(filepath)

	content = open(filepath).readlines();
	sneakPeak = content[0];

	if isASCII(sneakPeak):
		
		print '%s is encoded with ASCII' % filename

	elif whichBOM(sneakPeak):
		
		encoding = whichBOM(content[0])[0]

		filepath = BOM2CSV(filepath,newfile,encoding[0:6])

		print '%s is encoded with %s ' % (filename,encoding)

	else:

		print 'Parsing Error: Encoding for %s is unknown.'

		return

	# plate reader files typically keeps meta-data in the first few rows; skip them
	skiprows = determineLineSkips(filepath); 

	df = pd.read_csv(filepath,sep='\t',header=None,index_col=0,skiprows=skiprows);

	# interval is set at 600 by default; can be changed as input argument
	df.columns = listTimePoints(interval,df.shape[1])

	df.index.name = 'Well'
	df.T.index.name = 'Time'

	# newfile ends with .tsv extension by default, 
	# if plate reader output is also .tsv, this would over-write it
	df.to_csv(newfile, sep='\t')

	return df

def parsePlateName(plate):
    '''
    extracts isolate name and Biolog PM plate number from plate name 

    Args:

    	plate (string) - plate name

    Returns

    	isolate (string) - isolate name
    	pm (string) -- number, most likely 1 or 2 
    '''
    isolate = str(plate.split('_')[0]);
    pm = str(plate.split('PM')[1][0]);
    
    return isolate,pm

def whichBOM(data):
	'''
	If a string starts with a BOM marker, return marker if it corresponds to one of several UTF encoding types. 

	Reference https://unicodebook.readthedocs.io/guess_encoding.html

	Returns list of strings
	'''

	BOMS = (
    	(BOM_UTF8, "UTF-8"),
    	(BOM_UTF32_BE, "UTF-32-BE"),
    	(BOM_UTF32_LE, "UTF-32-LE"),
    	(BOM_UTF16_BE, "UTF-16-BE"),
    	(BOM_UTF16_LE, "UTF-16-LE"),
	)

	return [encoding for bom, encoding in BOMS if data.startswith(bom)]

