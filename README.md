# phenotypic-characterization


*Installation Note 1:* This repository is still under construction and the package has not been fully completed. Several features and applications have not been fully tested. Please use with caution and consult with me if you have any questions.

*Installation Note 2:* This repository works only with Python 2. A future release will modernize it for use in Python 3, but for now use Python 2.

*Installation Note 3:* If you are having trouble installing and running AMiGA on your machine, please let me know and I may be able to help. 

AMiGA is a python-based program that facilitates the high-throughput analysis of microbial growth data. It models growth curves with Gaussian Processes (GP) to infer microbial growth parameters such as maximum specific growth rate, doubling time, lag phase, and carrying capacity. It is especially useful for the analysis of Biolog Phenotypic Microarray (PM) data. The flexibility and utility of GP regression enables:
1. the analysis of microbial growth data that does not follow standard logistic or sigmoidal growth,
2. inference of non-standard microbial dynamics such as diauxic shifts, and
3. hypothesis-driven statistical testing of differences in microbial growth under different environmental conditions. 

AMiGA has been designed to be a minimalist, modular, and user-friendly program that allows for the analysis of single or multiple files in a single batch. It requires a single command line in the terminal. User arguments can be passed via the terminal or simply using the text-based parameter files described below.

## [Required] Download repository/code 

```git clone https://github.com/firasmidani/phenotypic-characterization.git```

or simply download as zip folder and extract. See green button on top right. 

## [Required] Python 2

If you are a Mac user, your machine will have Python installed. You can proceed to the following section. 

If you have not previously worked with python, I would recommending a python distribution such as <a href="http://docs.continuum.io/anaconda/">Anaconda</a>. See this useful <a href="https://fangohr.github.io/blog/installation-of-python-spyder-numpy-sympy-scipy-pytest-matplotlib-via-anaconda.html">guide</a> on installation of Python. If you run into difficulties with installing ```GPy``` using Anaconda, please try to set-up a local python environment. 

Note: AMiGA is written for use in Python 2. Future release will modernize it so that it can be run in Python 3. So please make sure that you are using Python 2. You can check the python version in your terminal:

> python --version

## [Optional] Set-up a local python environment 

Virtual environments allow you to create a virtual copy of your machine’s Python without affecting the set-up of the native Python. Accordingly, you can download modules/packages without affecting the dependencies for other applications that require Python. For more info, see (see <a href="https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/">here</a>). __I recommend that you follow the instructions on this linked page instead of mine below. The linked instructions are more thorough, up-to-date, and delineate differences between installations for Python 2 vs Python 3 and Windows vs Unix.__

**for macOS and linux users of Python 2**

1) Install virtualenv.

>```python -m pip install —-user virtualenv``` # if you are using the native Python on your machine
>
>or 
>
>```conda install virtualenv``` # if you are using Anaconda or Miniconda for Python

2) Setup the environment in the folder where you would like to save it. Here, I name the environment amiga.

>```python -m virtualenv /Users/firasmidani/example/amiga```

3) Activate the environment (you will need to do this everytime you run AMiGA)

>```source  /Users/firasmidani/example/amiga/bin/activate``` 

**for Windows users of Python 2** 

1) Install virtualenv.

>```python -m pip install —-user virtualenv``` # if you are using the native Python on your machine
>
>or 
>
>```conda install virtualenv``` # if you are using Anaconda or Miniconda for Python

2) Setup the environment in the folder where you would like to save it. Here, I name the environment amiga.

>```python -m virtualenv C:\\Users\firasmidani\example\amiga``` 

3) Activate the environment (you will need to do this everytime you run AMiGA)

>``` C:\\Users\firasmidani\example\amiga\Scripts\activate```

## [Required] Package dependencies or requirements


Install requirements.

```pip install -r requirements.txt``` if you are using the native Python on your machine

or 

```conda install --file requirements.txt```  # if you are using Anaconda or Miniconda for Python

See `requirements.txt` for full list of dependencies. 

If you have `matplotlib`, `seaborn`, `pandas`, `numpy`, `scipy`, `GPy`, you should be able to test AMiGA right away. The other packages in requirements.txt are dependencies for these main ones. Anaconda distributions typically have all of these except for `GPy`. You can try to install `GPy` in Anaconda as follows

```conda install gpy``` 

If this fails, you can install `GPy` in Anaconda with conda as follows (see <a href="https://docs.anaconda.com/anaconda/user-guide/tasks/install-packages/">documentation</a>):

```python
conda install -c conda-forge gpy
```

## How to set-up your working directory and format data?

See `instructions.pdf`. At the very bare minimum, you need a ```data``` folder and your data files should be saved inside it. Each data file should be structured as wells x time. The first column must be your Well ID (i.e. A1, B1, ... H11, H12).

## How to run AMiGA and pass arguments via text files?

See `instructions.pdf` for information on how to format your input data and pass arguments via text file. 

Call ```amiga.py``` with python and provide the only required argument that points to the working directory or individual filename in the working directory.

```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
```
or
```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/data/od_bacteria.asc
```
or to only plot the raw data
```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/data/od_bacteria.asc --plot-plate-only
```

## How to run AMiGA and pass arguments via the command line?

See ```instructions.pdf``` for information on how to format your input data and more details on the different parameters that AMiGA accepts.

Call ```amiga.py``` with python and provide the only required argument of input (```-i``` or ```--input```) that points to the working directory

```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
```
or
```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/data/od_bacteria.asc
````

Let's say you have many plates in your data directory, but you only want to  analyze a specifc subset of your data set. You can use the *subset* (```-s``` or ```-subset```) argument to specify the desired conditions. For example, if you are using Biolog plates, you can restrict analysis to speific set of isolates and substrates.

```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53;Substrate:Negative Control;alpha-D-Glucose'
```

Maybe some of the wells in your data were noisy, you can flag those wells with the *flag* argument (```-f``` or ```--flag```) as follows. 

```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-f 'PRB953_PM1-1:G10;PRB952_PM1-1:C3'
```

Of course, you can pass these arguments simultaneously.

```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53'
	-f 'PRB953_PM1-1:G10;PRB952_PM1-1:C3'
```

If you want to test a specific hypothesis with GP Regression, you can call it as follows with the *hypothesis* argument (```-h``` or ```--hypothesis```). This assumes a the null hypothesis (```OD ~ f(Time)```) and an alternative hypothesis (```OD ~ f(Time, Substrate)```).
```sh
python amiga.py 
	-i /Users/firasmidani/tecan/xra/ 
	-s 'Isolate:PRB952,PRB53;Substrate:Negative Control;alpha-D-Glucose'
	-h 'H0:Time;H1:Time+Substrate'
```
### What are the information provided by amiga output?

| Column header | Description                                                                                             |
|---------------|---------------------------------------------------------------------------------------------------------|
| Plate_ID      | Unique ID for each Biolog plate                                                                         |
| PM            | Biolog PM 1 or 2                                                                                        |
| Replicate     | Technical replicate                                                                                     |
| Min_OD        | Minimal OD of raw data                                                                                  |
| Max_OD        | Maximum OD of raw data                                                                                  |
| Baseline_OD   | OD at first time point of raw data                                                                      |
| Fold Change   | Maximum OD in a well divided by maximum OD in the negative control “A1” well                            |
| GP_r          | Maximum specific growth rate (i.e. exponential growth rate)                                             |
| GP_K          | Carrying capacity (should be close to GP_max)                                                           |
| GP_d          | Growth lag time ** I have yet to verify the validity of this parameters, so don’t dwell on it too much. |
| GP_AUC        | Area Under the Curve                                                                                    |
| GP_td         | Doubling time (in minutes)                                                                              |
| GP_max        | Max OD after log-transformation and subtraction of log OD at T=0                                        |

Note that GP_* indicates variables that were inferred after natural-log tranformation of OD data followed by subtraction of log OD at the first time point of each curve (i.e. log(OD(t)) / log(OD(0)) which is equivalent to OD(t) - OD(0)).

### Acknowledgements

Many thanks to the rest of the Biolog team in the Britton lab help in designing and building this workflow: James Collins, Ph.D., Heather Danhof, Ph.D., Colleen Brand, and Robert Britton, Ph.D. This work was supported by the National Institutes of Health (U01AI124290). 



