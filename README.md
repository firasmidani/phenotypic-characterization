# phenotypic-characterization

repository for wrangling and analysing data from biolog-based phenotypic characterization


## [Required] Download repository/code 

```git clone https://github.com/firasmidani/phenotypic-characterization.git```

or simply download as zip folder and extract. 

## [Required] Python

If you have not previously worked with python, I would recommending a python distribution such as <a href="http://docs.continuum.io/anaconda/">Anaconda</a> or <a href="https://www.spyder-ide.org/">Spyder</a>. See this useful <a href="https://fangohr.github.io/blog/installation-of-python-spyder-numpy-sympy-scipy-pytest-matplotlib-via-anaconda.html">guide</a> on installation of Python.

## [Optional] Set-up a local python environment 

**Make sure your computer has virtual environments (e.g. virtualenv) for Python (see <a href="https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/">here</a>)**
Virtual environments allow you to create a virtual copy of your machine’s Python without affecting the set-up the native python. This way you can download modules/packages without affecting the dependencies for other applications that require python.

on macOs and Linux: 

```python -m pip install —user virtualenv```

on Windows: 

```py -m pip install —user virtualenv```

**Setup the environment**

```virtualenv .```

```source bin/activate``` 

**Install requirements**

```pip install -r requirements.txt```

## [Required] Package dependencies or requirements

See requirements.txt for full list of dependencies. 

If you have `matplotlib`, `seaborn`, `pandas`, `numpy`, `scipy`, `GPy`, you should be able to test AMiGA right away. The other packages in requirements.txt are dependencies for these main ones. Anaconda distributions typically have all of these except for GPy. You can install packages with conda as follows (see <a href="https://docs.anaconda.com/anaconda/user-guide/tasks/install-packages/">documentaiton</a>):

```python
conda install package-name
```

## How to set-up your working directory?

See readme_metadata_and_parameters.pdf for information on how to pass arguments via text file. 

## How to run AMiGA and pass arguments via text files

See readme_metadata_and_parameters.pdf for information on how to format your input data and pass arguments via text file. 

Call ```amiga.py``` with python and provide the only required argument that points to the working directory.

```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
```

## How to run AMiGA and pass arguments via the command line

See readme_metadata_and_parameters.pdf for information on how to format your input data and more details on the different parameters that AMiGA accepts.

Call ```amiga.py``` with python and provide the only required argument of input (```-i``` or ```--input```) that points to the working directory

```python
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
```

Let's say you have many plates in your data directory, but you only want to  analyze a specifc subset of your data set. You can use the *subset* (```-s``` or ```-subset```) argument to specify the desired conditions. For example, if you are using Biolog plates, you can restrict analysis to speific set of isolates and substrates.

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
	-s 'Isolate:PRB952,PRB53;Substrate:Negative Control;D-Trehalose'
	-h 'H0:Time;H1:Time+Substrate'
```


