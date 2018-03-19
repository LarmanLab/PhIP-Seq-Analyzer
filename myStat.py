# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 09:38:51 2015

@author: yuan
"""

#standard modules
import math
import os
import sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.special import erf
#from scipy.interpolate import interp1d #lowess regression


#personal modules
#class: regression
#class: sigAnalysis
import myIO
import myPlot
import myDataframe
import myRegression


############
##################################
class sigAnalysis:
    def __init__(self,data=None):
        self.data = data
    
#'Cumulative distribution function for the standard normal distribution'
    def phi(self):
        phi_value = (1.0 + math.erf(float(self.data)/math.sqrt(2.0))) / 2.0
        if phi_value == 1.0:
            phi_adj = 1.0 - sys.float_info[3]
        elif phi_value == 0:
            phi_adj = sys.float_info[3]
        else:
            phi_adj = phi_value
        return phi_value, phi_adj
#         
    def median_zscores(self):
        m = np.median(self.data)
        s = np.std(self.data)
        zscores = [(x-m)/float(s) if s>0 else (x-m) for x in self.data]
        return zscores
#
    def mean_zscores(self):
        m = np.mean(self.data)
        s = np.std(self.data)
        zscores = [(x-m)/float(s) if s>0 else (x-m) for x in self.data]
        return zscores
        
#calculate z scores        
    def Z_test(self, mean, sd):
        self.data = list(self.data)
        #mean is mean of the population
        #sd is standard deviations
        #self.data should be list
        #calculate Z scores
        zscores = [(x-mean)/float(sd) for x in self.data]
        #print zscores

        #calculate p values with one tail
        #'Cumulative distribution function for the standard normal distribution'
        #left-tailed side
        leftpvals = [(1.0 + erf(x / math.sqrt(2.0))) / 2.0 for x in zscores]
        #print leftpvals
        #right-tailed side
        rightpvals = [1.0-x for x in leftpvals]
        leftlog = np.round(-np.log10(leftpvals),2)
        rightlog = np.round(-np.log10(rightpvals),2)
        #export
        zdf = pd.DataFrame({'x':self.data, 'zscores':np.round(zscores,2), 
                            'leftpvals':leftpvals, 'leftlog':leftlog, 
                            'rightpvals':rightpvals, 'rightlog':rightlog})
        #print 'Z scores(first 6 rows):\n', zdf.head(n=6)
        return zdf

        
#################################
#normalization methods
class normalization:
    #default: samples in columns, and genes/peptides/proteins in rows
    def __init__(self, par, infile, outfile, index_label='row_names'):
        self.par = par
        #read a file
        self.infile = infile
        insep = ',' if infile.endswith('.csv') else '\t'
        self.indf = pd.read_csv(infile, index_col=0, sep=insep, low_memory=False)
        #print self.indf
        #out file name
        self.outfile = outfile
        self.outsep = ',' if outfile.endswith('.csv') else '\t'
        self.index_label=index_label

#Formula normaliation by raw counts in columns
    def sd_func(self, x, scaling_factor):
        sd_x = np.std(x)
        if sd_x == 0: #outliers though seldemly happens
            norm_x = x
        else:
            norm_x = x*scaling_factor/sd_x
            norm_x = np.round(norm_x)
            norm_x = norm_x.astype(float)
        return norm_x
        
#Formula normaliation by raw counts in columns
    def norm_func(self, x, scaling_factor):
        sum_x = np.sum(x)
        if sum_x == 0: #outliers though seldemly happens
            norm_x = x
        else:
            norm_x = x*scaling_factor/sum_x
            norm_x = np.round(norm_x, 2)
            norm_x = norm_x.astype(float)
        return norm_x
        
#corrected RC    
    def negative_func(self, x):
            corrRC = x[1:] - x[0]
            corrRC[corrRC<0] = 0
            #print corrRC
            return corrRC

#divided by standard deviation of raw read counts and scaled by million
    def sd_scaling(self):
        #
        norm_df = self.indf.apply(self.sd_func, args=(self.par['scaling_factor'],), axis=0)
        #export norm_df
        print('export normalized data frame into ', self.outfile)
        norm_df.to_csv(self.outfile, index_label=self.index_label, sep=self.outsep)
        #print self.par['sample_names']
        return norm_df
            
#divided by raw read counts and scaled by a million reads
    def RC_scaling(self):
        #
        norm_df = self.indf.apply(self.norm_func, args=(self.par['scaling_factor'],), axis=0)
        #export norm_df
        print('export normalized data frame into ', self.outfile)
        norm_df.to_csv(self.outfile, index_label=self.index_label, sep=self.outsep)
        #print self.par['sample_names']
        return norm_df

#substract by negative controls followed by divided by raw read counts and then scaled by a million reads
    def rawcounts_NC_correction(self):
        #raw counts scaling
        sdf = self.indf.apply(self.norm_func, args=(self.par['scaling_factor'],), axis=0)
        
        #corrected by negative controls        #
        NC_df = sdf[self.par['NC_samples']].copy()
        NC_median = NC_df.median(axis=1)#by rows
        sdf.insert(0, 'NG_median', NC_median)
        #correction process
        cdf = sdf.apply(self.negative_func, axis=1)
        #print cdf
            
        #export norm_df
        print('export normalized data frame into ', self.outfile)
        sdf.to_csv(self.outfile, index_label=self.index_label, sep=self.outsep)
        return sdf

    def select_NC(self):
        #select df with control samples
        NC_df = self.indf[self.par['NC_samples']].copy()
        NC_df['median'] = self.indf[self.par['NC_samples']].median(axis=1)
        NC_df['std'] = self.indf[self.par['NC_samples']].std(axis=1)
        NC_df['mean'] = self.indf[self.par['NC_samples']].mean(axis=1)
        return NC_df
    
#calculate z scores against negative controls
    def NC_zscores(self):
        #select df with control samples
        df0 = self.select_NC()
        #remove 0 mean and highest outlier
        median_99 = np.percentile(NC_median,99)
        #print mean_99
        df1 = df0.loc[(df0['median']>0) & (df0['median']<median_99),:].copy()
        #print df1
        #fit linear model
        lm = myRegression.linear(df1, ('median','std'), self.par['dir_stat']).linear()
        #print 'Linear regression:', lm['params']
        
        print('calculate z scores:')
        self.indf.insert(0, 'predicted_std', lm['df']['pred_std'])
        self.indf.insert(0, 'median', lm['df']['median'])
        #loops of data frame
        zscores_df = self.indf.apply(lambda x: sigAnalysis(x[2:]).Z_test(x[0],x[1])['zscores'], axis=1)
        zscores_df.columns = list(self.indf)[2:]
        #export
        myDataframe.basic(zscores_df).export_df(self.outfile, self.par['zscore_threshold'], self.index_label)
        return zscores_df

#calculate z scores against negative controls
    def NCPHIPzscores_linear(self):
        ##1: regressed std of negative controls
        print('\tLinear regression of std ~ median of negative controls')
        sdf = self.indf[self.par['sample_names']].copy()
        #select df with control samples
        NC_df = self.select_NC()
        #print NC_df
        #remove 0 mean and highest outlier
        median99 = np.percentile(NC_df['median'],99)
        NC_regress = NC_df.loc[(NC_df['median']>0)&(NC_df['median']<median99),:].copy()
        #print NC_regress
        #fit linear model
        lm = myRegression.linear(NC_regress, ('median','std'), self.par['dir_stat']).linear()
        #print 'Linear regression:', lm['params']
        
        #2:regressed PHIP values against negative controls
        print('\tLinear regression of medians of phip sample ~ medians of negative controls')
        reg_df = pd.DataFrame(np.zeros(shape=sdf.shape), columns=list(sdf), index=list(sdf.index))
        for sample in self.par['sample_names']:
            subdf = pd.DataFrame({'NC': NC_df['median'], 'phip':sdf[sample]})
            median99 = np.percentile(subdf['phip'], 99)
            #remove top 1% phipseq values
            sub_regress = subdf.loc[(subdf['NC']>0)&(subdf['phip']<median99),:]
            #linear regression
            #print sample
            #print sub_regress
            sample_dir = self.par['dir_result'] + sample + '/'
            phip_lm = myRegression.linear(sub_regress, ('NC', 'phip'), sample_dir).linear()
            #predicted phipRC
            reg_df[sample] = phip_lm['df']['pred_phip']
        #3:
        print('\t calculate z scores against NC:')
        #RC minus regress RC
        residuals_df = sdf - reg_df
        #then divied by regressed std of negative control in rows
        divisor=pd.Series(NC_df['pred_std'])
        divisor[divisor==0]=1
        zscores_df=residuals_df.div(divisor, axis=0)
        zscores_df=np.round(zscores_df,1)
        #print residuals_df.ix[1:4, 8:10]
        
        #export
        myDataframe.basic(zscores_df).export_df(self.outfile, self.par['zscore_threshold'], self.index_label)
        #
        return zscores_df
#
#calculate z scores against negative controls based on robust linear regression
    def NCPHIPzscores_RLM(self):
        ##1: regressed std of negative controls
        print('\tLinear regression of std~median of negative controls')
        #select df with control samples
        NC_df = self.select_NC()
        #print NC_df
        #fit robust linear model
        NC_df, NC_fit=myRegression.linear(NC_df, ('median', 'std'), self.par['dir_stat']).RLM()
        divisor=NC_df['pred_std']
        divisor[divisor<=0]=1
        
        #2:regressed PHIP values against negative controls
        print('\tLinear regression of sample-specific medians of negative controls')
        sdf = self.indf[self.par['sample_names']].copy()
        med_df = pd.DataFrame(np.zeros(shape=sdf.shape), columns=list(sdf), index=list(sdf.index))
        zscores_df = pd.DataFrame(np.zeros(shape=sdf.shape), columns=list(sdf), index=list(sdf.index))
        for sample in self.par['sample_names']:
            sample_dir = self.par['dir_result'] + sample + '/'
            phip_df=pd.DataFrame({'median':NC_df['median'], 'phip':sdf[sample]})
            #fit robust linear model
            phip_df, phip_fit=myRegression.linear(phip_df, ('median', 'phip'), sample_dir).RLM()
            #calculate z-score
            med_df[sample] = phip_df['pred_phip']
            zscores_df[sample] = (sdf[sample]-med_df[sample])/divisor
        #3:
        #RC minus regress RC
        #residuals_df=sdf-med_df        
        #then divied by regressed std of negative control in rows
        #zscores_df=residuals_df.div(divisor, axis=0)
        zscores_df = np.round(zscores_df, 1)
        #print residuals_df.ix[1:4, 8:10]
        myDataframe.basic(zscores_df).export_df(self.outfile, self.par['zscore_threshold'], self.index_label)
        return zscores_df
        
########
#calculate z scores against negative controls based on polynomial regression
    def NCPHIPzscores_PN(self):
        ##1: regressed std~mean of negative controls
        #regression of logstd~logmedian across 261 Beads-only file
        #self.par['file_NC'], self.par['scaling_factor']
        wNC, wNC_fit = self.NC_whole_std()
        #RLM of std~median of beads only of this dataset
        NC_df = pd.DataFrame({'wNC_median':wNC['median'],'wNC_std':wNC['std'],\
                    'mean':self.indf[self.par['NC_samples']].mean(axis=1), \
                    'median':self.indf[self.par['NC_samples']].median(axis=1), \
                    'std':self.indf[self.par['NC_samples']].std(axis=1)})
        pNC, pNC_fit = myRegression.linear(NC_df, ('median','std'), self.par['dir_QC']).RLM()
        
        #2:regressed PHIP values against negative controls
        print('\tLinear regression of sample-specific medians of negative controls')
        #reg_x = NC['median'].drop_duplicates()
        zscores_df = self.indf.copy()
        zscores_df[:]=0.0
        for sample in self.par['sample_names']:
            zdf=NC_df.copy()
            zdf['phip'] = self.indf[sample]
            #fit robust linear model
            sample_dir = '{}{}/'.format(self.par['dir_result'], sample)
            mNC, mfit = myRegression.linear(zdf, ('mean','phip'), sample_dir).RLM()
            #zscores
            zdf['pred_phip']=mNC['pred_phip']
            pred_phip=zdf['pred_phip']
            pred_phip[pred_phip<=0] = np.nan
            zdf['pred_logphip'] = np.log10(pred_phip)# work as x value of NCstd~median
            pred_std = 10**wNC_fit.predict({'logmedian':zdf['pred_logphip']})
            zdf['pred_std']=pred_std
            zscores_df[sample] = (zdf['phip'] - zdf['pred_phip'])/pred_std
            zdf['zscores'] = zscores_df[sample]
            #export zscore
            zdf.to_csv(sample_dir+'polynomial_median.csv', header=True, index_label=self.index_label)
        #3:export z scores into self.outfiles
        zscores_df.replace([np.nan, np.inf, -np.inf], 0, inplace=True)
        zscores_df = np.round(zscores_df, 1)
        myDataframe.basic(zscores_df).export_df(self.outfile, self.par['zscore_threshold'], self.index_label)
        return zscores_df
    
#read all beads-only, which are 261 till now
#self.par['file_NC'], self.par['scaling_factor']
    def NC_whole_std(self):
        print('\tPolynomial regression of std~median across ALL BEADS-ONLY.')
        file_prefix = '{}{}_'.format(self.par['dir_result'], myIO.file_os(self.par['file_NC']).name_prefix())
        norm_ncfile = file_prefix+'scalingRC.txt'
        if os.path.isfile(norm_ncfile):
            phip_nc = pd.read_csv(norm_ncfile, sep='\t', index_col=0, low_memory=False)
        else:
            phip_nc = normalization(self.par, self.par['file_NC'], norm_ncfile).RC_scaling()
        #print(phip_nc.shape)
        
        #summary of nc: mean and std
        NC=pd.DataFrame({'mean':phip_nc.mean(axis=1), 'median':phip_nc.median(axis=1), \
                            'std':phip_nc.std(axis=1), 'sum':phip_nc.sum(axis=1)})
        NC['median'][NC['median']==0] = np.nan
        NC['std'][NC['std']==0] = np.nan
        NC['logmedian'] = np.log10(NC['median'])
        NC['logstd'] = np.log10(NC['std'])
        #NC=NC.replace([np.inf, -np.inf], -10) #an extreme small value
        #
        #initiate reg_df for regression
        #fill out outliers
        reg_df = NC.loc[(NC['median']>0),:].copy()
        #order for polynomial regression
        reg_df = reg_df.sort_values(['logmedian'], ascending=True)

        #polynomial regression
        formula = 'logstd~logmedian+I(logmedian**2)+I(logmedian**3)'
        pn_model = smf.ols(formula, data=reg_df)
        pn_fit = pn_model.fit()
        #print(pn_fit.params)
        reg_df['pred_logstd'] = pn_fit.predict()
        reg_df['pred_std'] = 10**pn_fit.predict()
        NC['pred_logstd'] = pn_fit.predict({'logmedian':NC['logmedian']})
        NC['pred_std'] = 10**NC['pred_logstd']
         
        #refresh total log
        #params=dict(pn_fit.params)
        #NC_dict = dict([('polynomial_NC_std:' + x, params[x]) for x in params.keys()])
        #myIO.file_os(self.par['file_total_log'], '=').line_replace(NC_dict)
        #export fitting of std
        NC.to_csv(file_prefix+'polynomial_std.csv', header=True, index_label='row_names')
        #draw graph
        xm=round(np.nanmax(list(NC['logmedian'])))
        ym=round(np.nanmax(list(NC['logstd'])))
        plot_par={'df': NC[['logmedian','logstd']], 'xlim':(-.5,xm), 'ylim':(-.5,ym),\
                  'picfile':file_prefix+'polynomial_std.png', 'text':pn_fit.params }
        try:
            myPlot.plot(plot_par).regressionP(reg_df['logmedian'], reg_df['pred_logstd'])
        except ValueError:
            print('Failed to drawing pic and save into {}'.format(plot_par['picfile']))
        
        #return fitting model object
        return NC, pn_fit
#
#calculate z scores against median/sd of samples
    def sample_zscores(self):
        #
        def median_zscores(x):
            y = x[x>0]
            m = np.median(y)
            s = np.std(y)
            zscores = [(r - m)/float(s) if s > 0 else (r - m) for r in x]
            return zscores
        #calculate z-scores
        zscores_df = self.indf.apply(median_zscores, axis=0)
        zscores_df = np.round(zscores_df, 1)
        #export
        myDataframe.basic(zscores_df).export_df(self.outfile, self.par['zscore_threshold'], self.index_label)
        return zscores_df
##########
#end