# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 15:27:36 2017

@author: yuan
"""


#standard modules
#import math
import os
#import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
#from scipy.special import erf
#from scipy.interpolate import interp1d #lowess regression


#personal modules
#class: regression
#class: sigAnalysis
import myIO
import myPlot
import myDataframe
import myStat

#######################################
#
class linear:
    def __init__(self, data=None, xy=None, outdir=None):
        #for linearReg_1: data should be DataFrame
        #for calculate response using linearReg_1: data should be Series
        if data is not None:
            self.data = pd.DataFrame(data)
            #default y~x: y is response. x is predictor
            if xy is None:
                self.x_name, self.y_name = list(self.data)[:2]
            else:
                self.x_name, self.y_name = xy
            self.py_name = 'pred_' + self.y_name
        self.outdir = outdir
        #default export
        self.lm={'df':None, 'fit': None, 'params':None, 'slope':0, 'intercept':0, 'R2':0}
        
#export csv and picture after regression 
    def export_regress(self, file_prefix):
        #export fitting of std
        plot_par={'df':self.data[[self.x_name, self.y_name]], 
                  'xlim':(-.5, np.nanmax(self.data[self.x_name])), 
                  'ylim':(-.5, np.nanmax(self.data[self.x_name])), 
                  'text':self.lm['params'] }
        if self.outdir is not None:
            file_prefix = '{}{}_{}'.format(file_prefix, self.x_name, self.y_name)
            self.data.to_csv(file_prefix+'.csv', header=True, index_label='row_names')
            plot_par['picfile']=file_prefix+'.png'
        #draw graph
        try:
            myPlot.plot(plot_par).regressionP(self.reg_df[self.x_name], self.reg_df[self.py_name])
        except ValueError:
            print('Failed to drawding {}'.format(file_prefix))
#
#ordinary linear regression
    def linear(self):
        #default two columns: y~x
        #always the first column is x, and the second is y
        formula = self.y_name + '~' + self.x_name
        
        #####fit models
        try:
            #filler out zero points
            self.reg_df=self.data.loc[self.data[self.x_name]>0,:].copy()
            self.reg_df=self.reg_df.sort_values([self.x_name], ascending=True)
        
            #fit models
            model = smf.ols(formula=formula, data=self.reg_df)
            result = model.fit()
            self.lm['fit']=result
            #print 'Linear regression', result.summary()
            
            #get essential parameters: slope and intercept
            self.lm['params']=result.params
            if len(result.params) == 1:
                self.lm['slope'] = result.params
            else:
                self.lm['intercept'], self.lm['slope'] = result.params
            #correlation
            self.lm['R2'] = result.rsquared_adj
            #print 'R(correlation coefficient)=', R
            
            #predict
            self.reg_df[self.py_name] = result.predict()
            self.data[self.py_name] = result.predict({self.x_name:self.data[self.x_name]})
            self.lm['df'] = self.data
            #export
            self.export_regress(self.outdir+'linear_')
        except ValueError:
            pass
        return self.lm
        
#Robust linear regression. xy is tuple of x_name and y_name
    def Robust_linear(self):
        #print self.data.head
        if self.data.shape[0] < 3 and self.data.shape[1] < 3:
            print('Error: Points are fewer for modeling. ')
            return self.lm
        
        #####fit models
        #print self.data
        try:
            #fit robust linear model
            X=sm.add_constant(self.data[self.x_name])
            Y=self.data[self.y_name]
            model = smf.RLM(Y, X, M=sm.robust.norms.HuberT() )
            result = model.fit()
            self.data[self.py_name]=model.predict(result.params)
            self.lm['fit']=result
            self.lm['df']=self.data
            self.lm['params']=result.params
            #print 'Linear regression', result.summary()
            #get essential parameters: slope and intercept
            if len(result.params) == 1:
                self.lm['slope'] = result.params[0]
                self.lm['intercept'] = 0
            else:
                self.lm['intercept'] = result.params[0]
                self.lm['slope'] = result.params[1]
            #correlation
            self.lm['R2'] = result.rsquared
            #print 'R(correlation coefficient)=', R
        except ValueError:
            pass
        #
        return self.lm

#NC regression based one phipseq experiement
#xy is x_name and y_name for RLM
    def RLM(self):
        ##1: regressed std of negative controls
        #print('\tRobust Linear Regression of {}~{}.'.format(self.y_name, self.x_name))
        #print NC_df
        #filler out zero points
        self.reg_df=self.data.loc[self.data[self.x_name]>0,:].copy()
        self.reg_df=self.reg_df.sort_values([self.x_name], ascending=True)
        
        #fit robust linear model
        X = sm.add_constant(self.reg_df[self.x_name])
        Y = self.reg_df[self.y_name]
        std_model = smf.RLM(Y, X, M=sm.robust.norms.HuberT() )
        sfit = std_model.fit()
        self.lm['fit']=sfit
        #print 'Linear regression:', std_result.summary()
        #predict
        self.reg_df[self.py_name] = std_model.predict(sfit.params)
        self.data[self.py_name] = std_model.predict(sfit.params, sm.add_constant(self.data[self.x_name]))
        self.data[self.py_name][self.data[self.py_name]==0]=np.inf
        self.lm['df']=self.data
        #export
        self.export_regress(self.outdir+'RLM_')
        #
        return self.data, sfit


#lowess regression
#frac is between 0 and 1. The fraction of the data used when estimating each y-value.
    #def lowess(self, frac=0.6):
            