# phenotypic-characterization

AMiGA is a python-based program that facilitates the high-throughput analysis of microbial growth data. It models growth curves with Gaussian Processes (GP) to infer microbial growth parameters such as maximum specific growth rate, doubling time, lag phase, and carrying capacity. It is especially useful for the analysis of Biolog Phenotypic Microarray (PM) data. The flexibility and utility of GP regression enables:
1. the analysis of microbial growth data that does not follow standard logistic or sigmoidal growth,
2. inference of non-standard microbial dynamics such as diauxic shifts, and
3. hypothesis-driven statistical testing of differences in microbial growth under different environmental conditions. 

AMiGA has been designed to be a minimalist, modular, and user-friendly program that allows for the analysis of single or multiple files in a single batch. It requires a single command line in the terminal. User arguments can be passed via the terminal or simply using the text-based parameter files described below.

## [Required] Download repository/code 

```git clone https://github.com/firasmidani/phenotypic-characterization.git```

or simply download as zip folder and extract. See green button on top right. 

## [Required] Python

If you have not previously worked with python, I would recommending a python distribution such as <a href="http://docs.continuum.io/anaconda/">Anaconda</a>. See this useful <a href="https://fangohr.github.io/blog/installation-of-python-spyder-numpy-sympy-scipy-pytest-matplotlib-via-anaconda.html">guide</a> on installation of Python. If you run into difficulties with installing ```GPy``` using Anaconda, please try to set-up a local python environment. 

## [Optional] Set-up a local python environment 

Virtual environments allow you to create a virtual copy of your machine’s Python without affecting the set-up the native python. This way you can download modules/packages without affecting the dependencies for other applications that require python. For more info, see (see <a href="https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/">here</a>). __I recommend that you follow the instructions in the link instead of mine below, becasue the linked instructions are more thorough, up-to-date, and delineate differences between installations for Python 2 vs Python 3 and Windows vs Unix.__

**To install virtualenv**

on macOs and Linux: 

```python -m pip install —-user virtualenv```

on Windows: 

```py -m pip install —-user virtualenv```

**Setup the environment**

On macOS and Linux

```virtualenv .```
or
```python -m virtualenv env```

On Windows

```py -m venv env``` # this is will use "env" as the name of the virtual environment

**Activate the environment (you will need to do this everytime you run AMiGA**

On macOS and Linux

```source env/bin/activate``` 

On Windows

```.\env\Scripts\activate```

**Install requirements**

```pip install -r requirements.txt```

## [Required] Package dependencies or requirements

See `requirements.txt` for full list of dependencies. 

If you have `matplotlib`, `seaborn`, `pandas`, `numpy`, `scipy`, `GPy`, you should be able to test AMiGA right away. The other packages in requirements.txt are dependencies for these main ones. Anaconda distributions typically have all of these except for `GPy`. You can install `GPy` in Anaconda with conda as follows (see <a href="https://docs.anaconda.com/anaconda/user-guide/tasks/install-packages/">documentation</a>):

```python
conda install -c conda-forge gpy
```

## How to set-up your working directory and format data?

See `instructions.pdf`. At the very bare minimum, you need a ```data``` folder and your data files should be saved inside it. Each data file should be structured as wells x time. The first column must be your Well ID (i.e. A1, B1, ... H11, H12).

## How to run AMiGA and pass arguments via text files?

See `instructions.pdf` for information on how to format your input data and pass arguments via text file. 

Call ```amiga.py``` with python and provide the only required argument that points to the working directory or individual filename in the working directory.

```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
```
or
```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/data/od_bacteria.asc
```

## How to run AMiGA and pass arguments via the command line?

See ```instructions.pdf``` for information on how to format your input data and more details on the different parameters that AMiGA accepts.

Call ```amiga.py``` with python and provide the only required argument of input (```-i``` or ```--input```) that points to the working directory

```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
```
or
```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/data/od_bacteria.asc
````

Let's say you have many plates in your data directory, but you only want to  analyze a specifc subset of your data set. You can use the *subset* (```-s``` or ```-subset```) argument to specify the desired conditions. For example, if you are using Biolog plates, you can restrict analysis to speific set of isolates and substrates.

```python

python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53;Substrate:Negative Control;alpha-D-Glucose'
```

Maybe some of the wells in your data were noisy, you can flag those wells with the *flag* argument (```-f``` or ```-flag```) as follows. 

```python

python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-f 'PRB953_PM1-1:G10;PRB952_PM1-1:C3'
```

Of course, you can pass these arguments simultaneously.

```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53'
	-f 'PRB953_PM1-1:G10;PRB952_PM1-1:C3'
```

If you want to test a specific hypothesis with GP Regression, you can call it as follows with the *hypothesis* argument (```-h``` or ```-hypothesis```). This assumes a the null hypothesis (```OD ~ f(Time)```) and an alternative hypothesis (```OD ~ f(Time + Substrate)```).
```python

python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53;Substrate:Negative Control;alpha-D-Glucose'
	-h 'H0:Time;H1:Time+Substrate'
```

### Acknowledgements

Many thanks to the rest of the Biolog team in the Britton lab for help in designing and building this workflow: James Collins, Ph.D., Heather Danhof, Ph.D., Colleen Brand, and Robert Britton, Ph.D. This work was supported by the National Institutes of Health (U01AI124290). 



