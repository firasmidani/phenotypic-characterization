# Copyright 2014-2016 by Marco Galardini.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# MODIFIFED BY FIRAS MIDANI ON 05-03-2019

"""Growth curves fitting and parameters extraction for phenotype data.

This module provides functions to perform sigmoid functions fitting to
Phenotype Microarray data. This module depends on scipy curve_fit function.
If not available, a warning is raised.

Functions:
logistic           Logistic growth model.
gompertz           Gompertz growth model.
richards           Richards growth model.
guess_plateau      Guess the plateau point to improve sigmoid fitting.
guess_lag          Guess the lag point to improve sigmoid fitting.
fit                Sigmoid functions fit.
get_area           Calculate the area under the PM curve.
"""

import numpy as np
import pandas as pd

try:
    from scipy.optimize.minpack import curve_fit
    from scipy.integrate import trapz
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        'Install scipy to extract curve parameters.')


def logistic(x, A, u, d, v, y0):
    """Logistic growth modelself.

    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A / (1 + np.exp((((4 * u) / A) * (d - x)) + 2))) + y0
    return y


def gompertz(x, A, u, d, v, y0):
    """Gompertz growth model.

    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A * np.exp(-np.exp((((u * np.e) / A) * (d - x)) + 1))) + y0
    return y


def richards(x, A, u, d, v, y0):
    """Richards growth model (equivalent to Stannard).

    Proposed in Zwietering et al., 1990 (PMID: 16348228)
    """
    y = (A * pow(1 + (v + (np.exp(1 + v) * np.exp((u / A) *
                                                  (1 + v) * (1 + (1 / v)) * (d - x)))), -(1 / v))) + y0
    return y


def guess_lag(x, y):
    """Given two axes returns a guess of the lag point.

    The lag point is defined as the x point where the difference in y
    with the next point is higher then the mean differences between
    the points plus one standard deviation. If such point is not found
    or x and y have different lengths the function returns zero.
    """
    if len(x) != len(y):
        return 0

    diffs = []
    indexes = range(len(x))

    for i in indexes:
        if i + 1 not in indexes:
            continue
        diffs.append(y[i + 1] - y[i])
    diffs = np.array(diffs)

    flex = x[-1]
    for i in indexes:
        if i + 1 not in indexes:
            continue
        if (y[i + 1] - y[i]) > (diffs.mean() + (diffs.std())):
            flex = x[i]
            break

    return flex


def guess_plateau(x, y):
    """Given two axes returns a guess of the plateau point.

    The plateau point is defined as the x point where the y point
    is near one standard deviation of the differences between the y points to
    the maximum y value. If such point is not found or x and y have
    different lengths the function returns zero.
    """

    if len(x) != len(y):
        return 0

    diffs = []
    indexes = range(len(y))

    for i in indexes:
        if i + 1 not in indexes:
            continue
        diffs.append(y[i + 1] - y[i])
    diffs = np.array(diffs)

    ymax = y[-1]
    for i in indexes:
        if y[i] > (ymax - diffs.std()) and y[i] < (ymax + diffs.std()):
            ymax = y[i]
            break

    return ymax

def guess_rate(x,y):

	total_time = x[-1]-x[0];
	total_growth = max(y)-min(y);

	ini_r = float(total_growth) / float(total_time);

	return ini_r


def fit(function, x, y):
    """Fit the provided function to the x and y values.

    The function parameters and the parameters covariance.
    
    NOTE: curve_fit is very sensitive to initial parameters. for growth curves, ini_u is key
          modify so that you try multiple initial guess until you find the best fit --> see optimize_initial_u()
    """

    # Compute guesses for the parameters
    # This is necessary to get significant fits

    ini_u = 1e-4 # biopython default of 4 leads to the wrong minimal optima (I use 1e-4 sometimes)
    ini_v = 0.1; # neither gompertz or logisitc use this parameters, so meh

    p0 = [guess_plateau(x, y), optimize_initial_u(function,x,y), guess_lag(x, y), ini_v, min(y)]
    #p0 = [guess_plateau(x, y), guess_rate(x,y), guess_lag(x, y), ini_v, min(y)]
    #p0 = [guess_plateau(x, y), ini_u, guess_lag(x, y), ini_v, min(y)]

    # often, y0 is estimated to be really low while is estimated to be really high, 
    #        visual fit is good but parameter fit is awful
    p0_bounds = ([-np.inf,-np.inf,-np.inf,-np.inf,min(y)],
                 [np.inf,np.inf,np.inf,np.inf,np.inf]);

    params, pcov = curve_fit(function, x, y, p0=p0,bounds=p0_bounds,maxfev=10000,check_finite=True)
    
    return params, pcov

def optimize_initial_u(function,x,y):
    '''Try pre-determined initial guess of u (i.e. growth rate) and return the one that results in lowest SSE
    '''
    ini_u_list = [1e-3,1e-2,1e-1,1e0,1e1,1e2];
    sse_list = [];

    ini_v = 0.1;

    for ini_u in ini_u_list:
        
        p0 = [guess_plateau(x, y), ini_u, guess_lag(x, y), ini_v, min(y)]

        try: 
            params, pcov = curve_fit(function, x, y, p0=p0,maxfev=10000)
            y_true = y;
            y_pred = [function(xx,*params) for xx in x];
            sse_list.append(SumSquaredError(y_true,y_pred));
        except:
            sse_list.append(np.inf)

    ind = np.where(sse_list==np.min(sse_list))[0];
    ini_u = ini_u_list[ind];

    return ini_u

def SumSquaredError(y_true,y_pred):
    '''Returns sum of squared erros

    '''
    
    df = pd.DataFrame([y_true,y_pred])
    e = df.diff()
    se = e.apply(lambda x: x**2)
    sse = np.sqrt(se.sum(1).sum(0))

    return sse

def get_area(y, x):
    """Get the area under the curve."""
    return trapz(y=y, x=x)

