# phenotypic-characterization

repository for wrangling and analysing data from biolog-based phenotypic characterization


## Download Repository (code) [required]**

```git clone https://github.com/firasmidani/phenotypic-characterization.git```

## Set-up a local python environment [optional]**

**Make sure your computer has virtual environments for Python (see <a href="https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/">here</a>)**
Virtual environments allow you to create a virtual copy of your machine’s python without affecting the set-up the native python. This way you can download modules/packages without affecting the dependencies for other applications that require python.

on macOs and Linux: 

```python -m pip install —user virtualenv```

on Windows: 

```py -m pip install —user virtualenv```

**Setup the environment**

```virtualenv .```


```source bin/activate``` 

**Install requirements**

```pip install -r requirements.txt```

## Package dependencies

If you have matplotlib, seaborn, pandas, numpy, scipy, GPy, you should be able to test it right away. The other requirements are dependencies for the before-mentioned packages. 

**Requirements**
make sure that your python environment has the following requirements:
```
backports.functools-lru-cache==1.5
cycler==0.10.0
decorator==4.4.0
GPy==1.9.8
kiwisolver==1.1.0
matplotlib==2.2.4
numpy==1.16.5
pandas==0.24.2
paramz==0.9.5
pyparsing==2.4.2
python-dateutil==2.8.0
pytz==2019.3
scipy==1.2.2
seaborn==0.9.0
six==1.12.0
subprocess32==3.5.4
```



## How to run AMiGA and pass arguments via text files

simply call ```amiga.py``` with python and provide the only required argument  that points to the working directory

```python
python amiga.py 
	-i /Users/firasmidani/Downloads/20191011/ 
```

See readme_metadata_and_parameters.pdf for information on how to pass arguments via text file. 


## How to run AMiGA and pass arguments via the command line

**Basic Usage**

simply call ```amiga.py``` with python and provide the only required argument  that points to the working directory

```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
```

Let's say you have 20 plates in your data directory, but you want to only analyze a specifc subset. You can use the *subset* (```-s``` or ```-subset```) argument to specify these conditions. If you are using Biolog plates, you can for example only include certain isolates and certain substrates in your analysis.

```python

python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53;Substrate:Negative Control;D-Trehalose'
```

Maybe some of the wells in your data were noisy, you can flag those wells with the *flag* argument (```-f``` or ```-flag```) as follows. 

```python

python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-f 'PRB953_PM1-1:G10;PRB952_PM1-1:C3'
```

Of course, you can also both subset and flag.
```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53'
	-f 'PRB953_PM1-1:G10;PRB952_PM1-1:C3'
```

If you want to run a specific hypothesis, you can call it as follows with the *hypothesis* argument (```-h``` or ```-hypothesis```)
```python

python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53;Substrate:Negative Control;D-Trehalose'
	-h 'H0:Time;H0:Time+Substrate'
```


