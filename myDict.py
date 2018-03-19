# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 15:24:45 2015

@author: yuan
"""
import collections
import itertools
import string
import pandas as pd
import numpy as np
#
import myList
#
class basic:

    def __init__(self, dictionary=None):
        self.dict = dictionary
        if isinstance(self.dict, dict):
            self.sorted_keys=myList.basic(self.dict.keys()).sort_list()
        self.out_dict={}

#initiate nested dictionary as all zeros all something else    
    def init_dict2(self, keys1, keys2, default_value=0):
        for k1 in keys1:
            self.out_dict[k1] = {}
            for k2 in keys2:
                self.out_dict[k1][k2] = default_value
        return self.out_dict
        
    def print_dict(self):
        n = 1
        for key in self.sorted_keys:
            print('{:5}: {:10}\t{}'.format(n, key, self.dict[key]))
            n += 1
            
    def print_dict2(self):
        keys2 = self.nested_keys()
        #header line
        print('{:10}\t{}'.format('row_name', '\t'.join(map(str, keys2))))
        #rows loop
        for key1 in self.sorted_keys:
            items = [key1]
            sub_dict = self.dict[key1]
            for key2 in keys2:
                value = sub_dict[key2] if key2 in sub_dict else '0'
                items.append(value)
            print('{:10}\t{}'.format(key1, '\t'.join(map(str, items))))
        
    def dict_to_file(self, out_file, pattern='='):
        out_obj = open(out_file, 'wt')
        for key in self.sorted_keys:
            line = ''.join([str(key),str(pattern),str(self.dict[key]),'\n'])
            out_obj.write(line)
        out_obj.close()
        print('Write a dictionary to ', out_file)

    def nested_keys(self):
        nkeys = []
        for key in self.dict.keys():
            sub_dict = self.dict[key]
            sub_keys = sub_dict.keys()
            #combine lists and get unique elements
            nkeys = list(set(nkeys) | set(sub_keys))
        nkeys = sorted(list(map(str, nkeys)))
        #print nkeys
        return nkeys
        
    #dict2 is nested dictionary    
    def dict2_to_file(self, out_file, pattern='\t', NA='0', index_label='row_name', row_names=None):
        keys2 = self.nested_keys()
        #print keys2
        #
        out_obj = open(out_file, 'wt')
        #header line
        header = index_label+pattern+pattern.join(map(str, keys2))+'\n'
        out_obj.write(header)
        #rows loop
        if row_names is None:
            row_names = self.sorted_keys
        #print row_names
        for key1 in row_names:
            if isinstance(key1, str):
                key1 = key1.rstrip(string.whitespace) #remove white space characters
            items = [key1]
            sub_dict = self.dict[key1]
            for key2 in keys2:
                value = sub_dict[key2] if key2 in sub_dict else '0'
                items.append(value)
            line = pattern.join(map(str, items))+'\n'
            #print key1, sub_dict
            out_obj.write(line)
        out_obj.close()
        print('write a nested dictionary to ', out_file)

    #dict2 is nested dictionary, but the file only three columns
    #key1, key2, and value
    #not header
    def flat_dict2_to_file(self, out_file, pattern=','):
        #open file
        out_obj = open(out_file, 'wt')
        for key1 in self.dict.keys():
            if isinstance(key1, str):
                key1 = key1.rstrip(string.whitespace) #remove white space characters
            sub_dict = self.dict[key1]
            for key2 in sub_dict.keys():
                if isinstance(key2, str):
                    key2 = key2.rstrip(string.whitespace) #remove white space characters
                value = sub_dict[key2]
                line = pattern.join(map(str, [key1, key2, value]))+'\n'
                #print key1, sub_dict
                out_obj.write(line)
        out_obj.close()
        print('write a nested dictionary to ', out_file)
        
    #row to column and column to row
    def transform_dict2(self):
        #transform
        for key1, sub_dict in self.dict.items():
            for key2, value in sub_dict.items():
                if key2 not in self.out_dict:
                    self.out_dict[key2] = {}
                self.out_dict[key2][key1] = value
        #print self.out_dict
        return self.out_dict
        
    #key to value and value to key           
    def reverse_dict(self):
        self.out_dict = dict( (value, key) for key, value in self.dict.items())
        return self.out_dict

    #key to value and multiple values to multiple keys
    def reverse2_dict(self,sep=','):
        for key, values_str in self.dict.items():
            values = values_str.split(sep)
            for v in values:
                self.out_dict[v] = key
                #print v+key
        return self.out_dict
        
    def counting_reversed_dict(self):
        for key, value in self.dict.items():
            if type(value) is list:
                value = ','.join(value)
            self.out_dict[value] = self.out_dict.setdefault(value, 0)+1
        return self.out_dict
        
    #normalization scaling by the sum read counts in row
    def normalize_df(self):
        self.dict = pd.DataFrame(self.dict)
        #sum_RC = self.dict.apply(np.sum, axis=0)
        #
        def norm_func(x):
            sum_x = np.sum(x)
            norm_x = x*1e6/sum_x
            norm_x = np.round(norm_x)
            norm_x = norm_x.astype(int)
            return norm_x
        norm_df = self.dict.apply(norm_func, axis=0)
        #print norm_df
        return norm_df
        
#elements frquency of compound dict(value is list)  
#such as input: {'a':'04','b':'053','c':'04','d':'04','e':'1'}, output: {'04':3,'053':1,'1':1}
#return is frequency dict
    def elements_frequency(self, selected_keys=None, sep=','):
        counts = {} # frequency values only
        details = {} # items
        if selected_keys is None:
            for key, value_list in self.dict.items():
                for value in value_list:
                    counts[value] = counts[value] + 1 if value in counts else 1
                    details[value] = details[value] + sep + str(key) if value in details else str(key)
        else:
            #x=len([k for k,v in Counter(selected_keys).items() if v>1] )
            #print "%s: %s" % (len(selected_keys),  x)
            for key in selected_keys:
                value_list = self.dict[key]
                for value in value_list:
                    counts[value] = counts[value]+1 if value in counts else 1
                    details[value] = details[value] + sep + str(key) if value in details else str(key)
        #print details
        return counts,details
        
#update nested dictionary
    def update_dict2(self, dictB):
        new_d = collections.defaultdict(dict)
        #
        for k,v in itertools.chain(self.dict.items(),dictB.items()):
            new_d[k].update(v)
        #print new_d
        return new_d

#combine nested dictionary: the two keys should be the same
#the values would be a list of the values from the two dict
#tag is used for recrusive combination more than 2 times
    def combine_dupdict2(self, dictB, tag=0):
        #self.dict is not empy
        if bool(self.dict):
            if tag >= 2:
                for key1,nested_dict in dictB.items():
                    for key2, value in nested_dict.items():
                        self.dict[key1][key2].append(value)
            else:
                for key1,nested_dict in dictB.items():
                    for key2, value in nested_dict.items():
                        self.dict[key1][key2] = [self.dict[key1][key2], value]
        #the first dict is empty
        else:
            self.dict = dictB
        #print self.dict
        return self.dict
        
#the keys and length of the two dict can be different
    def combine_dict2(self, dictB):
        #
        for key1, nested_dict in itertools.chain(self.dict.items(), dictB.items()):
            #initiate 
            if key1 not in self.out_dict: 
                self.out_dict[key1] = {}
            for key2, value in nested_dict.items():
                if key2 not in self.out_dict[key1]: 
                    self.out_dict[key1][key2] = [value]
                else:
                    self.out_dict[key1][key2].append(value)
        #print self.out_dict
        return self.out_dict
            
#end