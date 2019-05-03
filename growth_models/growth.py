#!/usr/bin/env python

# Author: Firas Midani
# Date Initial: 2019-04-22
# Date Last Modified: 2019-04-30

# DESCRIPTION Framework for processing plate reader data at the Britton Laboratory (BCM) and inferring growth dynamics. 

# CLASS GrowthPlate
#     | -- init
#         | -- key
#         | -- control
#         | -- data
#         | -- time
#     | -- smoothData
#     | -- subtractControl
#     | -- subtractBaseline
#     | -- convertTimeUnits
#     | -- extractGrowthData

# CLASS GrowthData
#     | -- init
#         | -- key
#         | -- data
#         | -- time 
#     | -- plot

# CLASS GrowthMetrics
#     | -- init
#         | -- key
#         | -- data
#         | -- time
#     | -- Classical
#     | -- inferClassicalDynamics
#     | -- inferClassicalAUC
#     | -- inferDoublingTime
#     | -- inferGP_r
#     | -- inferGP_d (In Progress)
#     | -- inferGP_K
#     | -- inferGP_AUC
#     | -- inferGPDynamics
#     | -- GP
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

# 1. key for objects should store whether modificiations (e.g. subtractControl or smoothData) were applied to object


# IMPORT NECESSARY LIBRARIES

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import GPy
import imp
import os

from scipy.signal import savgol_filter
import scipy.stats as stats

# IMPORT IN-HOUSE LIBRARIES

lib_path = '/Users/firasmidani/Downloads/phenotypic-characterization/growth_fitting_library.py'

foo = imp.load_source('growth_fitting_library',lib_path);

from growth_fitting_library import fit, gompertz, logistic

# UTILITY FUNCTIONS 

def gpDerivative(x,gp):

    # from Solak et al. <-- from Tonner et al. (github.com/ptonner/gp_growth_phenotype/)
    mu,_ = gp.predictive_gradients(x);
    _,cov = gp.predict(x,full_cov=True);
    mult = [[((1./gp.kern.lengthscale)*(1-(1./gp.kern.lengthscale)*(y-z)**2))[0] for y in x] for z in x];
    
    return mu,mult*cov   

# BEGIN FRAMEWORK

class GrowthPlate(object):
    
    def __init__(self,data=None,key=None,control='A1'):
        
        '''
        Data structure for handling growth data for a microtiter plate (biolog)
        
        Attributes:
        
        data (pd.DataFrame): n by (p+1) DataFrame (for n timepoints and p wells) that stores the raw optical density data.
                             The first column is assumed to be time. 
                             Column headers are assumed to be well positions (e.g. A1)
        key (pd.DataFrame): p by k DataFrame (for k experimental variables) that stores experimental variables for each well. 
                             The index column is assumed to be well position (e.g. A1)
        control (str): position of control well in the plate. 

        '''
        
        data = data.copy().sort_values(['Time'],ascending=True);

        self.key = key.copy();
        self.control = control;
        self.data = data.iloc[:,1:];
        self.time = pd.DataFrame(data.iloc[:,0]);
        self.input_time = self.time.copy()
        self.input_data = self.data.copy()

        self.mods = pd.DataFrame(columns=['smoothed','floored','controlled','logged'],index=['status']);
        self.mods = self.mods.apply(lambda x: False);

        assert type(key) == pd.DataFrame, "key must be a pandas dataframe"
        assert data.columns[0]=='Time', "first data column must be Time"
        assert (data.shape[1]-1) == (key.shape[0]), "key must contain metadata for each sample"
        
    #enddef 

    def logData(self):

        self.data = self.data.apply(lambda x: np.log(x))
        self.mods.logged = True
        
    def smoothData(self,window=19,polyorder=3):
        
        self.data = self.data.apply(lambda x: savgol_filter(x,window,polyorder), axis=0)
        self.mods.smoothed = True

    def subtractBaseline(self):
                        
        self.data = self.data.apply(lambda x: x-self.data.iloc[0,:],axis=1)
        self.mods.floored = True

    def subtractControl(self):
                
        self.data = self.data.apply(lambda x: x-self.data.loc[:,self.control],axis=0)
        self.mods.controlled = True

    def convertTimeUnits(self):
        '''        
        converts seconds to hours
        
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
        
    #enddef

class GrowthData(object):
    
    def __init__(self,time=None,data=None,key=None,mods=None):
        
        '''
        Data structure for handling growth data for a single sample
        
        Attributes:
        
        time (pd.DataFrame): n by 1 DataFrame (for n timepoints) that stores time.
        data (pd.DataFrame): n by 1 DataFrame (for n timepoints) that stores the raw optical density data.
        key (pd.DataFrame): 1 by k DataFrame (for k experimental variables) that stores experimental variables for each well. 
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
    
    def plot(self):
        
        fig,ax = plt.subplots(figsize=[4,4]);
        
        ax.plot(self.time,self.data,lw=5,color=(0,0,0,1));
        
        [ii.set(fontsize=20) for ii in ax.get_xticklabels()+ax.get_yticklabels()];
        
        ax.set_xlabel('Time',fontsize=20);

        if self.mods.logged:
            ax.set_ylabel('log(OD)',fontsize=20);
        else:
            ax.set_ylabel('Optical Density',fontsize=20);
        ax.set_title(self.key.Substrate[0],fontsize=20);
       
        return fig,ax

    def log(self):
        '''
        natural logarithm transform 
        '''
        self.data = self.data.apply(lambda x: np.log(x+1e-3))
        self.mods.logged = True
    
   
class GrowthMetrics(object):
    
    def __init__(self,growth=None,model=None,params=None,pred=None):
        
        
        self.time = growth.time.copy();
        self.data = growth.data.copy();
        self.key = growth.key.copy();
        self.model = str(model)

    def Classical(self,model=None):
        
        if model==None:
            print("User must define choice of classical model from logistic, gompertz, or richards")
            pass
            
        x = np.ravel(self.time.values);
        y = np.ravel(self.data.values);
        
        try:
        	params, pcov = fit(model,x,y)
        
      		self.model = model;
        	self.params = params;
        except:
        	self.model = model;
        	self.params = [np.nan,np.nan,np.nan]
        
    def inferClassicalDynamics(self):
        
    	model = 'classical'

        self.key['%s_r' % model] = self.params[1]
        self.key['%s_K' % model] = self.params[0]
        self.key['%s_d' % model] = self.params[2]
        self.key['%s_AUC' % model] = self.inferClassical_AUC();
        self.key['%s_td' % model] = self.inferDoublingTime(mtype=model);
        
    def inferClassical_AUC(self):

        x = np.ravel(self.time);
        y = np.ravel(self.data);
        
        return np.trapz(y,x)

    def inferDoublingTime(self,mtype='classical'):
        '''
        assumes growth rate is per hour

        # reaffirmed by Swain et al. 2016. Nature Comm.
        #  The doubling time is ln(2) times the inverse of the growth rate. 
`
        '''

        r = self.key['%s_r' % mtype]

        #r = (np.log10(2)/r)*60; 
        r = (np.log(2)/r)*60; # doubling time in minutes

        return r

    def inferGP_r(self):
        
        x = self.time.values
        gp = self.model;        
        
        mu,cov = gpDerivative(x,gp)
        
        ind = np.where(mu==mu.max())[0];

        return mu[ind,0,0][0],np.diag(cov)[ind][0],ind
    
    def inferGP_K(self):
    
        x = self.time.values
        gp = self.model; 
        
        mu,cov = self.model.predict(x,full_cov=True);
        ind = np.where(mu==mu.max())[0];

        return mu[ind,0][0],np.diag(cov)[ind][0],ind
    
    def inferGP_AUC(self):
        
        x = self.time.values
        gp = self.model; 
        
        mu,cov = gp.predict(x,full_cov=True);
        ind = np.where(mu==mu.max())[0];
        
        dt = np.mean(x[1:,0]-x[:-1,0]);
        D = np.repeat(dt,mu.shape[0]).T;
        
        mu = np.dot(D,mu)[0];
        var = np.dot(D,np.dot(cov,D));
        
        return mu,var
    
    def inferGP_d(self,threshold=.95):
        
        x = self.time.values
        gp = self.model; 
        
        mu,var = gpDerivative(x,gp);
        #mu,cov = gp.predict(x,full_cov=True)
        
        #print mu
        #print cov
        #print mu[:,:,0]
        #print var
        #print var[:,0]
        #print mu[:,0]
        #print cov[:,0]
        #print np.sqrt(cov[:,0][0])
        #print stats.norm.cdf(0,mu[:,0][0],np.sqrt(cov[:,0][0]))
        
#         for ii,t in enumerate(zip(mu[:,0],cov[:,0])):
#             print ii,
#             m,v = t[0],t[1]
#             print m,v
#             print x[ii]
#             z = stats.norm.cdf(x[ii],loc=m,scale=np.sqrt(v))
#             print z
#             z = z
            
        prob = np.array([stats.norm.cdf(0,loc=m,scale=np.sqrt(v))[0] for m,v in zip(mu[:,:,0],var[:,0])]);
        #prob = np.array([stats.norm.cdf(0,loc=m,scale=np.sqrt(v))[0] for m,v in zip(mu[:,0],cov[:,0])]);

        ind = 0;
        while ind < prob.shape[0] and prob[ind] > threshold: 
            ind += 1
        if ind == prob.shape[0]:
            ind -= 1
            
        return x[ind]
    
    def inferGPDynamics(self):
        
        #print self.inferGP_r()
        #print self.inferGP_K()
        #print self.inferGP_AUC()
        #print self.inferGP_d()
        
        self.key['GP_r'] = self.inferGP_r()[0]
        self.key['GP_K'] = self.inferGP_K()[0]
        #self.key['GP_d'] = self.inferGP_d()[0]
        self.key['GP_AUC'] = self.inferGP_AUC()[0]
        self.key['GP_td'] = self.inferDoublingTime(mtype='GP');
            
    def GP(self):
                
        x = self.time; 
        y = self.data.values
        
        k = GPy.kern.RBF(x.shape[1],ARD=True)

        m = GPy.models.GPRegression(x,y,k)
        m.optimize()
        #print m
        
        self.model = m;
        #self.params = params;
        
    def predictClasscial(self):
        
        x = np.ravel(self.time.values);

        #K,r,d,v,y0 = self.params
        
        self.pred = [gompertz(xx,*self.params) for xx in x];
    
    def predictGP(self):
        
        x = self.time.values;
        
        #K,r,d,v,y0 = self.params
        
        self.pred = self.model.predict(x)[0]
        
                
    def plot(self):
        
        fig,ax = plt.subplots(figsize=[4,4]);
        
        ax.plot(self.time,self.data,lw=5,color=(0,0,0,0.65));
        ax.plot(self.time,self.pred,lw=5,color=(1,0,0,0.65));
        
        [ii.set(fontsize=20) for ii in ax.get_xticklabels()+ax.get_yticklabels()];
        
        ax.set_xlabel('Time',fontsize=20);
        ax.set_ylabel('Optical Density',fontsize=20);
        
        ax.set_title(self.key.Substrate[0],fontsize=20);
       
        return fig,ax    