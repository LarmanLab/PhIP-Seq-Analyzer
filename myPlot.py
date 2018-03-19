# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 13:34:53 2015

@author: yuan
"""
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#in-house models
import myList
plt.rcParams['agg.path.chunksize'] = 20000

############################33
class plot:
    def __init__(self, par):
        #data frame
        if 'df' in par.keys():
            self.data = par['df']
            self.nrow = self.data.shape[0]
            self.ncol = self.data.shape[1]
            self.colnames = self.data.columns
            self.xlabel = par['xlabel'] if 'xlabel' in par.keys() else list(self.data)[0]
            self.ylabel = par['ylabel'] if 'ylabel' in par.keys() else list(self.data)[1]
            self.xlim = par['xlim'] if 'xlim' in par.keys() else myList.basic(self.data.ix[:,0]).min_max()
            self.ylim = par['ylim'] if 'ylim' in par.keys() else myList.basic(self.data.ix[:,1]).min_max()
            #initiate plot window
            plt.clf()
        elif 'list' in par.keys():
            self.data = par['list']
            self.xlabel = par['xlabel'] if 'xlabel' in par.keys() else 'x'
            self.ylabel = par['ylabel'] if 'ylabel' in par.keys() else 'y'
            self.xlim = par['xlim'] if 'xlim' in par.keys() else None
            self.ylim = par['ylim'] if 'ylim' in par.keys() else None
        else:
            self.data=None
            print('Error:No data frame input as long as drawing a plot!')
        #file
        self.picfile = par['picfile'] if 'picfile' in par.keys() else None
        #colors
        self.col = 'bgrcmykw'
        #line styles:solid, dashed, dotted, dashdot
        self.title = par['title'] if 'title' in par.keys() else 'Plot'
        self.text = par['text'] if 'text' in par.keys() else None
        self.legend = par['legend'] if 'legend' in par.keys() else False
        self.pch = par['pch'] if 'pch' in par.keys() else 'o'
        self.lty=par['lty'] if 'lty' in par.keys() else 'solid'
        self.lwd=par['lwd'] if 'lwd' in par.keys() else 1

#default is dataframe with at least two columns
#lty is linestyle, lwd=line width
#x_value=0, data of x axis always come from 0 column no matter how many y axis
    def lineP(self, x_value=0):
        plt.title(self.title) 
        if x_value == 0:
            for col in range(1, self.ncol):
                 plt.plot(self.data.ix[:,0], self.data.ix[:,col], color=self.col[col], 
                          linestyle=self.lty, linewidth=self.lwd, label=self.colnames[col])
        else:
            for col1 in range(0, self.ncol, 2):
                 col2 = col1+1
                 plt.plot(self.data.ix[:,col1], self.data.ix[:,col2], linestyle=self.lty, linewidth=self.lwd)
        plt.xlabel(self.colnames[0])
        plt.ylabel('-'.join(self.colnames[1:]) )
        if self.legend is not None:
            #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.legend(loc=self.legend)
        if self.xlabel is not None:
            plt.xlabel(self.xlabel)
        if self.ylabel is not None:
            plt.ylabel(self.ylabel)
        #export
        self.show_pic()

    #default is dataframe with at least two columns
    def dotP(self):
        plt.title(self.title) 
        #pch:1:-,2:|,4:<,5:>6:,8,o:solid circle;
        for col in range(1, self.ncol):
            plt.scatter(self.data.ix[:,0], self.data.ix[:,col], color=self.col[col], marker=self.pch)
        #plt.xlabel(self.colnames[0])
        #plt.ylabel('-'.join(self.colnames[1:]) )
        if self.legend is True:
            plt.legend()
        if self.text is not None:
            y_max=np.max(np.max(self.data.ix[:,1:]))
            plt.text(0, y_max, self.text)
        if self.xlabel is not None:
            plt.xlabel(self.xlabel)
        if self.ylabel is not None:
            plt.ylabel(self.ylabel)
        #export
        self.show_pic()
            
#draw a plot with observed and linear regression dots.
    def regressionP(self, regressed_x, regressed_y):
        plt.clf()
        #observed dots
        plt.plot(self.data.ix[:,0], self.data.ix[:,1],'b.')
        #regressed lines
        plt.plot(regressed_x, regressed_y,'r-')
        #other parameters
        #OverflowError: Allocated too many blocks if missing the below
        plt.xlim(self.xlim)
        plt.ylim(self.ylim)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        if self.text is not None:
            x=self.xlim[1]*.6
            y=self.ylim[1]*.6
            plt.text(x, y, self.text)
        #export
        self.show_pic()
            
#draw a pie plot
    def pieP(self, labels=True, values=True):
        plt.title(self.title) 
        #label
        if labels is True:
            labels = list(self.data.index)
        if values is True:
            values = '%.2f'
        
        #
        plt.pie(self.data.ix[:,0],labels=labels,autopct=values) 
        #if legend is True:
        #    plt.legend(loc=3, labels=labels)
        #export
        self.show_pic()
            
#draw a barplot: transfer data frame into series type
    #verticle bars
    def simple_barv(self, col_index=0):
        #parameters of barplot
        labels=list(self.data.index)
        y_pos = np.arange(len(labels))
        #
        plt.bar(y_pos, self.data.ix[:,col_index], align='center', alpha=0.5, color=self.col[0])
        plt.set_xlim = ([0,len(labels)])
        plt.xticks(y_pos, labels)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        #export
        self.show_pic()
            
#horizontal bars
    def simple_barh(self, col_index=0):
        #parameters of barplot
        labels = list(self.data.index)
        x_pos = np.arange(len(labels))
        #
        plt.barh(x_pos, self.data, align='center', alpha=0.5, color=self.col[0])
        plt.set_ylim = ([0,len(labels)])
        plt.yticks(x_pos, labels)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        #export
        self.show_pic()

#
    def show_pic(self):
        #export
        if self.picfile is None:
            plt.show()
        else:
            print("Draw and save plots into ", self.picfile)
            plt.savefig(self.picfile, dpi=600)
            plt.close() 
        
################
#end