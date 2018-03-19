#!/usr/bin/env
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 10:12:34 2015

@author: yuan
"""
#standar modules
#import multiprocessing as mp   # multi process
#import multiprocessing.dummy as mpd #multi-threading
import os
import re
#personal modules
import myDict
import myIO
#class
#class:samples

class samples:
    def __init__(self, par):
        self.par = par
        #print self.par['dir_raw_data']
        #print self.par['file_sample_info']
        #print self.par['dir_result']
        #print self.par['dir_result_array']
        if 'dir_result' in self.par and 'dir_result_array' not in self.par:
            self.par['dir_result_array'] = self.par['dir_result']
        self.raw_sample = {}
        self.sample_raw = {}
        self.sample_names = []
        self.sample_dirs = {}

    def seek_fq(self,dir_raw_data):
        print('Retrieve all *.fastq files under', dir_raw_data)
        raw_files = []
        #get all files
        all_files = myIO.dir_os(dir_raw_data).recrusive_files() 
        #find file with .fastq or .fq
        for af in all_files:
            m = re.search(r'fastq$|fq$', af)
            if m:
                #print 'raw data:',af
                raw_files.append(af)
        return raw_files
    
    #sample names based on raw files    
    def raw_to_samples(self):
        #get all fastq files
        raw_files = self.seek_fq(self.par['dir_raw_data'])
        #print raw_files
        
        #connect raw file to sample name
        for raw_file in raw_files:
            sample_name = myIO.file_os(raw_file).name_prefix()
            self.raw_sample[raw_file] = sample_name
            if sample_name in self.sample_raw:
                self.sample_raw[sample_name].append(raw_file)
            else:
                self.sample_raw[sample_name] = [raw_file]
        #print dict
        #myDict.basic(self.sample_raw).print_dict()
        #myDict.basic(self.raw_sample).print_dict()

    #sample names based on sample_info.csv
    def file_to_samples(self):
        #get all fastq files
        raw_files = self.seek_fq(self.par['dir_raw_data'])
        print('Number of raw files:', len(raw_files))
        #read sample info file
        print('Read sample file: ', self.par['file_sample_info'])
        in_obj = open(self.par['file_sample_info'], 'rt')

        #set connections between raw data and sample_name
        for line in in_obj:
            line = line.rstrip("\n")
            items = line.split(',')
            raw_file_name = items[0]
            sample_name = items[1]
            #print prefix
            for raw_file in raw_files:
                file_name = myIO.file_os(raw_file).file_name()
                #print file_name
                if file_name.find(raw_file_name) == 0:
                    #dict: raw_sample
                    self.raw_sample[raw_file] = sample_name
                    #dict: sample_raw
                    if sample_name in self.sample_raw:
                        self.sample_raw[sample_name].append(raw_file)
                    else:
                        self.sample_raw[sample_name] = [raw_file]   
        in_obj.close()
        #print dict
        #myDict.basic(self.sample_raw).print_dict()
        #myDict.basic(self.raw_sample).print_dict()
        
    def sample_info(self):
        sample_pairs = {}
        for raw_file, sample_name in self.raw_sample.items():
            raw_file_name = myIO.file_os(raw_file).file_name()
            group = 'NC' if 'BEADS' in raw_file_name.upper() else 'PhIP'
            if not 'unassigned' in raw_file_name:
                sample_name = re.sub('_R1', "", sample_name)
                pair = '{},{}'.format(raw_file_name, sample_name)
                sample_pairs[pair]=group
        #export dict to file
        print('Generate sample file: ', self.par['file_sample_info'])
        #order per record: fastq file name, sample_name, phip_group
        myDict.basic(sample_pairs).dict_to_file(self.par['file_sample_info'], ',')

#storage given sample_dir        
    def sample_storage(self):
        result_dirs = self.par['dir_result_array'].split(',')
        result_dirs_pool = result_dirs*len(self.sample_names)
        #print '===', result_dirs_pool
        #    
        for sample_name in self.sample_raw:
            #search if sample dir exists
            flag = 0
            for result_dir in result_dirs:
                sample_dir = result_dir+sample_name+'/'
                if os.path.isdir(sample_dir):
                    self.sample_dirs[sample_name] = sample_dir
                    flag=1
                    break
            #assign this sample with a directory
            if flag == 0:
                given_result_dir = result_dirs_pool.pop()
                sample_dir = given_result_dir+sample_name+'/'
                self.sample_dirs[sample_name] = sample_dir
        #myDict.basic(self.sample_dirs).print_dict()
    
    def group_names(self,col=3):
        if col < 3: col = 3
        #
        group_samples = {}
        sample_groups = {}
        in_obj = open(self.par['file_sample_info'],'rt')
        for line in in_obj:
            line = line.rstrip("\n")
            items = line.split(',')
            if len(items) < col:
                print('No group name in Column ', col)
                break
            sample_name, group_name = items[1], items[col-1]
            sample_groups[sample_name]=group_name
            if group_name in group_samples:
               group_samples[group_name] += ',' + sample_name
            else:
               group_samples[group_name] = sample_name
        in_obj.close()
        return group_samples, sample_groups
        
    #TOP function        
    def export_sample_info(self):
        print('get sample_raw and raw_sample')
        if os.path.isfile(self.par['file_sample_info']):
            self.file_to_samples()
        else:
            self.raw_to_samples()
            #generate sample file
            self.sample_info()
        #
        print("\nSample and raw files:")
        myDict.basic(self.sample_raw).print_dict()
        self.par['raw_to_sample'] = self.raw_sample # one raw file vs one sample name
        self.par['sample_to_raw'] = self.sample_raw #one sample name vs a list of raw files
        #get sample names
        self.sample_names = sorted(self.sample_raw.keys())
        self.par['sample_names'] = self.sample_names
        #get sample_dirs
        self.sample_storage()
        self.par['sample_dirs'] = self.sample_dirs

        print('get group names if exists')
        flag = 1
        while flag > 0:
            group_samples, sample_groups=self.group_names(flag+2)
            if group_samples == {}:
                flag = 0
            else:
                key = 'group' + str(flag)
                self.par[key] = group_samples
                flag += 1
                print('Groups of {}: {}'.format(key,self.par[key].keys()))
        return self.par
    
########
#end    
