#!/usr/bin/env python

# Author: Firas Midani
# Date Initial: 2019-04-22
# Date Last Modified: 2019-05-07

# DESCRIPTION Framework for processing plate reader data at the Britton Laboratory (BCM) and inferring growth dynamics. 

# CLASS GrowthPlate
#     | -- init
#         | -- key
#         | -- control
#         | -- data
#         | -- time
#         | -- input_time
#         | -- input_data
#         | -- mods
#     | -- convertTimeUnits
#     | -- extractGrowthData
#     | -- logData
#     | -- plot
#     | -- smoothData
#     | -- subtractControl
#     | -- subtractBaseline

# CLASS GrowthData
#     | -- init
#         | -- key
#         | -- data
#         | -- time 
#         | -- mods 
#     | -- log
#     | -- plot

# CLASS GrowthMetrics
#     | -- init
#         | -- key
#         | -- data
#         | -- time
#         | -- mods
#         | -- classical_model
#     | -- fitClassical
#     | -- inferClassicalDynamics
#     | -- inferClassicalAUC
#     | -- inferDoublingTime
#     | -- inferGP_r
#     | -- inferGP_d (In Progress)
#     | -- inferGP_K
#     | -- inferGP_AUC
#     | -- inferGPDynamics
#     | -- fitGP
#     | -- predictClassical
#     | -- predictGP
#     | -- plot
#
# DEF gpDerivative(x,gp)

# ACKNOWELDGMENTS

# Parts of this framework include code snippets are inspired from work by Peter Tonner (github.com/ptonner/gp_growth_phenotype/)

# NOTES ON STRUCTURE OF THIS FRAMEWORK

# Why is GrowthMetrics a class? By making GrowthMetrics a separate class from GrowthData, you can avoid unintential mutability of raw grwoth data. You can also make different GrowthMetrics objects for the same GrowthData (e.g. both gompertz and logistical) and still maintain the parent GrowthData in its unmodified form. 

# Why does the library not include funtions for input/output of data? This framework is motivated by analysis of microbial cultures grown on Biolog Phenotypic Characterizaiton PM1 and PM2 plates (with pre-determined layout of substrates). Still, this framework needs to be applicable for non-Biolog based analyses. By de-coupling input/output of data and dat analysis, I can enable easier modularity of code. 

# TO DO

# 1. predictClassical() should not use gompertz by default

# IMPORT NECESSARY LIBRARIES

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import GPy
import imp
import os

from scipy.signal import savgol_filter
import scipy.stats as stats

# IMPORT IN-HOUSE LIBRARIES

sys.path.append('..')

from classical import fit, gompertz, logistic
from plates import plotPlateGrowth

# UTILITY FUNCTIONS 

def gpDerivative(x,gp):

    if isinstance(x,pd.DataFrame):
        x = x.values[:,np.newaxis]
    elif isinstance(x,list):
        x = np.array(x)[:,np.newaxis];
    elif isinstance(x,np.ndarray) & (x.ndim ==1):
        x = x[:,np.newaxis];

    # from Solak et al. <-- from Tonner et al. (github.com/ptonner/gp_growth_phenotype/)
    mu,_ = gp.predictive_gradients(x);
    _,cov = gp.predict(x,full_cov=True);
    mult = [[((1./gp.kern.lengthscale)*(1-(1./gp.kern.lengthscale)*(y-z)**2))[0] for y in x] for z in x];
    
    return mu,mult*cov   

# BEGIN FRAMEWORK

class GrowthPlate(object):
    
    def __init__(self,data=None,key=None,time=None,mods=None,control='A1',):
        
        '''
        Data structure for handling growth data for a microtiter plate (biolog)
        
        Attributes:
        
        data (pd.DataFrame): n by (p+1) DataFrame (for n timepoints and p wells) that stores the raw optical density data.
                             The first column is assumed to be time. 
                             Column headers are assumed to be well positions (e.g. A1)
        key (pd.DataFrame): p by k DataFrame (for k experimental variables) that stores experimental variables for each well. 
                             The index column is assumed to be well position (e.g. A1)
        control (str): position of control well in the plate. 
        mods (pd.DataFrame): 4 by 1 DataFrame that stores the status of transformation or modifications of the object data set. 
        time (pd.DataFrame): n by 1 DataFrame (for n timepoints) that stores the raw time points
        '''
        
        data = data.copy().sort_values(['Time'],ascending=True);

        self.key = key.copy();
        self.control = control;

        if time is None:
            self.data = data.iloc[:,1:];
            self.time = pd.DataFrame(data.iloc[:,0]);
        else:
            self.data = data.copy()
            self.time = time.copy();

        self.input_time = self.time.copy()
        self.input_data = self.data.copy()

        self.mods = pd.DataFrame(columns=['smoothed','floored','controlled','logged'],index=['status']);
        self.mods = self.mods.apply(lambda x: False);

        assert type(key) == pd.DataFrame, "key must be a pandas dataframe"
        assert data.columns[0]=='Time', "first data column must be Time"
        assert (data.shape[1]-1) == (key.shape[0]), "key must contain metadata for each sample"
        
    #enddef 

    def logData(self):
        '''Transform with a natural logarithm all data points'''

        self.data = self.data.apply(lambda x: np.log(x))
        self.mods.logged = True
        
    def smoothData(self,window=19,polyorder=3):
        '''Smooth each array (column) with a Savitzky-Golay filter'''

        self.data = self.data.apply(lambda x: savgol_filter(x,window,polyorder), axis=0)
        self.mods.smoothed = True

    def subtractBaseline(self):
        '''Subtract first value in each array (column) from all elements of the array

        NOTE: if performed after logging, this is equivalent to scaling relative to OD at T0!
        '''

        self.data = self.data.apply(lambda x: x-self.data.iloc[0,:],axis=1)
        self.mods.floored = True

    def subtractBaseline(self):
        '''Subtract first value in each array (column) from all elements of the array

        NOTE: if performed after logging, this is equivalent to scaling relative to OD at T0!
        '''

        self.data = self.data.apply(lambda x: x-self.data.iloc[0,:],axis=1)
        self.mods.floored = True

    def subtractControl(self):
        '''Subtract array (column) belonging to control well from all wells'''

        self.data = self.data.apply(lambda x: x-self.data.loc[:,self.control],axis=0)
        self.mods.controlled = True

    def convertTimeUnits(self):
        '''Convert time array (column) from units of seconds to hours
        '''
        
        self.time = self.time.astype(float)/3600        
                        
    def extractGrowthData(self,arg_dict={}):
        '''
        NOTE: argument is framed as a dictionary to allow for non-biolog formats of plates but could simply be substrate
        
        '''

        # make sure variable selections are not empty
        if not bool(arg_dict):

            print("Error: Selection of variables must be defined by a dictionary.")
            
            return None
        
        # make sure variable selections are formated as a dictionary of lists or np.arrays
        for dict_key,value in arg_dict.iteritems():

            if (len(value)==1) and not (isinstance(value,(list,np.ndarray))):

                arg_dict[dict_key] = [value];


        sub_key = self.key[self.key.isin(arg_dict).sum(1)==len(arg_dict)];
        sub_key_idx = sub_key.index;
        
        sub_data = self.data.loc[:,sub_key_idx];
        sub_time = self.time;

        sub_mods = self.mods
                        
        return GrowthData(sub_time,sub_data,sub_key,sub_mods)

    def plot(self):
        '''
        only works if data is a whole 96-well plate
        '''

        df = self.data.copy().T;
        df.columns = np.ravel(self.time.copy().values);
        summary = self.key;

        fig,axes = plotPlateGrowth(df,summary,threshold=1.5,title="",savefig=0,filepath="",logged=self.mods.logged);

        return fig,axes

class GrowthData(object):
    
    def __init__(self,time=None,data=None,key=None,mods=None):
        
        '''
        Data structure for handling growth data for a single sample
        
        Attributes:
        
        time (pd.DataFrame): n by 1 DataFrame (for n timepoints) that stores time.
        data (pd.DataFrame): n by 1 DataFrame (for n timepoints) that stores the raw optical density data.
        key (pd.DataFrame): 1 by k DataFrame (for k experimental variables) that stores experimental variables for each well. 
        mods (pd.DataFrame): 4 by 1 DataFrame that stores the status of transformation or modifications of the object data set. 
        '''
        
        self.time = time.copy();
        self.data = data.copy();

        self.input_time = time.copy();
        self.input_data = data.copy();

        self.key = key.copy();

        if mods is None:
            self.mods = pd.DataFrame(columns=['smoothed','floored','controlled','logged'],index=['status']);
            self.mods = self.mods.apply(lambda x: False);
        else:
            self.mods = mods

        assert type(time) == pd.DataFrame, "time must be a pandas dataframe"
        assert type(data) == pd.DataFrame, "data must be a pandas dataframe"
        assert type(key) == pd.DataFrame, "key must be a pandas dataframe"

    #enddef 
    
    def plot(self,fig=None,ax=None):
        
        if not ax:
            fig,ax = plt.subplots(figsize=[4,4]);
        
        ax.plot(self.time,self.data,lw=5,color=(0,0,0,1));
        
        [ii.set(fontsize=20) for ii in ax.get_xticklabels()+ax.get_yticklabels()];
        
        ax.set_xlabel('Time',fontsize=20);

        if self.mods.logged:
            ax.set_ylabel('log(OD)',fontsize=20);
        else:
            ax.set_ylabel('Optical Density',fontsize=20);
        
        if 'Substrate' in self.key.keys():
            ax.set_title(self.key.Substrate[0],fontsize=20);
       
        return fig,ax

    def logData(self):
        '''
        natural logarithm transform 
        '''

        if not self.mods.logged:
            self.data = self.data.apply(lambda x: np.log(x+1e-3))
            self.mods.logged = True
        else: 
            print 'WARNING: data has already been log-transformed.'

    def smoothData(self,window=19,polyorder=3):
        '''Smooth each array (column) with a Savitzky-Golay filter'''

        if not self.mods.smoothed:
            self.data = self.data.apply(lambda x: savgol_filter(x,window,polyorder), axis=0)
            self.mods.smoothed = True
        else: 
            print 'WARNING: data has already been smoothed.'

    def subtractBaseline(self):
        '''Subtract first value in each array (column) from all elements of the array

        NOTE: if performed after logging, this is equivalent to scaling relative to OD at T0!
        '''

        if not self.mods.floored:
            self.data = self.data.apply(lambda x: x-self.data.iloc[0,:],axis=1)
            self.mods.floored = True
        else: 
            print 'WARNING: baseline has already been subtracted'
   
class GrowthMetrics(object):
  
    def __init__(self,growth=None,classical_model=None,gp_model=None,params=None,pred=None,mods=None):
        '''
        Data structure for summarizing bacterial growth curves
        
        Attributes:
        
        data (pd.DataFrame): n by 1 DataFrame (for n time points) for optical density data
        time (pd.DataFrame): n by 1 DataFrame (for n time points) for time points
        key (pd.DataFrame): p by k DataFrame (for k experimental variables) that stores experimental variables for each well. 
                             The index column is assumed to be well position (e.g. A1)
        classical_model (function): classical growth model function (logistic, gompertz, or richards), see growth_fitting_library.py
        gp_model (?)
        mods (pd.DataFrame): 4 by 1 DataFrame that stores the status of transformation or modifications of the object data set. 
        pred (pd.DataFrame): n by 1 DataFrame (for n time points) for model predicted optical density data
        params (np.array): array of floats for parameters inferred by classical model fitting (K,r,d,v,y0)
            where K is carrying capacity, r is growth rate, d is lag itme, v is a shape parameter, and y0 is intiial 
        '''  

        self.time = growth.time.copy();
        self.data = growth.data.copy();
        self.key = growth.key.copy();
        self.mods = growth.mods.copy();
        self.classical_model = classical_model;

    def fitClassical(self,classical_model=None):
        '''Fits OD to a classical model of bacterial growth'''
       
        if classical_model==None:
            print("User must define choice of classical model from logistic, gompertz, or richards")
            pass
            
        x = np.ravel(self.time.values);
        y = np.ravel(self.data.values);
        
        try:
        	params, pcov = fit(classical_model,x,y)
        
      		self.classical_model = classical_model;
        	self.params = params;
        except:
        	self.classical_model = classical_model;
        	self.params = [np.nan,np.nan,np.nan,0.1,y[0]]
       
    def fitGP(self):
        '''Fit a Gaussian Process Regression model'''
        x = self.time; 
        y = self.data.values
        
        k = GPy.kern.RBF(x.shape[1],ARD=True)

        m = GPy.models.GPRegression(x,y,k)
        m.optimize()
        #print m
        
        self.gp_model = m;
        #self.params = params;

    def inferClassicalDynamics(self):
        '''Infers then stores growth parameters (r,K,d,AUC,td) in key'''
        
    	model = 'classical'

        self.key['%s_r' % model] = self.params[1]
        self.key['%s_K' % model] = self.params[0]#+self.params[4] #A-y0 so assumes that y0 is zero ?
        self.key['%s_d' % model] = self.params[2]
        self.key['%s_AUC' % model] = self.inferClassical_AUC();
        self.key['%s_td' % model] = self.inferDoublingTime(mtype=model);
        
    def inferClassical_AUC(self):
        '''Infers area under the curve using classical approach'''

        x = np.ravel(self.time);
        y = np.ravel(self.data);
        
        return np.trapz(y,x)

    def inferDoublingTime(self,mtype='classical'):
        '''
        assumes growth rate is per hour

        # reaffirmed by Swain et al. 2016. Nature Comm.
        #  The doubling time is ln(2) times the inverse of the growth rate. 
        '''

        r = self.key['%s_r' % mtype]

        #r = (np.log10(2)/r)*60; 
        r = (np.log(2)/r)*60; # doubling time in minutes

        return r

    def inferGP_r(self):
        '''Infers growth rate with GP regression'''
        
        x = self.time.values
        gp = self.gp_model;        
        
        mu,cov = gpDerivative(x,gp)
        
        ind = np.where(mu==mu.max())[0];

        return mu[ind,0,0][0],np.diag(cov)[ind][0],ind
    
    def inferGP_K(self):
        '''Infers carrying capacity with GP regression'''
        x = self.time.values
        gp = self.gp_model; 
        
        mu,cov = self.gp_model.predict(x,full_cov=True);
        ind = np.where(mu==mu.max())[0];

        return mu[ind,0][0],np.diag(cov)[ind][0],ind
    
    def inferGP_AUC(self):
        '''Infers area under the curve with GP regression'''
        x = self.time.values
        gp = self.gp_model; 
        
        mu,cov = gp.predict(x,full_cov=True);
        ind = np.where(mu==mu.max())[0];
        
        dt = np.mean(x[1:,0]-x[:-1,0]);
        D = np.repeat(dt,mu.shape[0]).T;
        
        mu = np.dot(D,mu)[0];
        var = np.dot(D,np.dot(cov,D));
        
        return mu,var
    
    def inferGP_d(self,threshold=.95):
        '''
        works but delta is based on time.values (i.e. discretized, not continous) - 2019.05.07
        '''
        x = self.time.values
        gp = self.gp_model; 
        
        mu,var = gpDerivative(x,gp);
            
        prob = np.array([stats.norm.cdf(0,loc=m,scale=np.sqrt(v))[0] for m,v in zip(mu[:,:,0],var[:,0])]);

        ind = 0;
        while ind < prob.shape[0] and prob[ind] > threshold: 
            ind += 1
        if ind == prob.shape[0]:
            ind -= 1
            
        return x[ind]
    
    def inferGPDynamics(self):
        '''Infers then stores growth parameters (r,K,AUC,td) in key'''
        #print self.inferGP_r()
        #print self.inferGP_K()
        #print self.inferGP_AUC()
        #print self.inferGP_d()
        
        self.key['GP_r'] = self.inferGP_r()[0]
        self.key['GP_K'] = self.inferGP_K()[0]
        self.key['GP_d'] = self.inferGP_d()[0]
        self.key['GP_AUC'] = self.inferGP_AUC()[0]
        self.key['GP_td'] = self.inferDoublingTime(mtype='GP');
        
    def predictClassical(self,classical_model=gompertz):
        '''Predict OD using classical model'''
        
        x = np.ravel(self.time.values);

        #K,r,d,v,y0 = self.params
        
        #### THIS SHOULDNOT PREDETERMINED AS GOMPERTZ !!!!!! ####
        self.pred = [classical_model(xx,*self.params) for xx in x];
        self.key['classical_max'] = np.max(self.pred)

    def predictGP(self):
        '''Predict OD using GP regression'''
        x = self.time.values;
        
        #K,r,d,v,y0 = self.params
        
        self.pred = np.ravel(self.gp_model.predict(x)[0])
        self.key['GP_max'] = np.max(self.pred)
               
    def plot(self,ax=None):
        
        if not ax:
            fig,ax = plt.subplots(figsize=[4,4]);
         
        ax.plot(self.time,self.data,lw=5,color=(0,0,0,0.65));
        ax.plot(self.time,self.pred,lw=5,color=(1,0,0,0.65));
        
        [ii.set(fontsize=20) for ii in ax.get_xticklabels()+ax.get_yticklabels()];
        
        ax.set_xlabel('Time',fontsize=20);
        ax.set_ylabel('Optical Density',fontsize=20);
        
        ax.set_title(self.key.Substrate[0],fontsize=20);
       
        return fig,ax    