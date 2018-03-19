# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 14:03:33 2015

@author: yuan
"""
#myIO means operations of input-output related to files or directories
#modules
import gzip #compress/uncompress file
import os
import shutil
import sys
import subprocess   #multi-threads
import re       #regular expression
import urllib   #web service
#import xlwt     #excel file reading
import numpy as np
import pandas as pd

#personal modules
import myDict
#
#class file_os
#class dir_os
#class file_IO
############################################################################
#
class dir_os:
    def __init__(self, indir=None):
        #format input dir
        self.indir=indir
        self.format_dir()
        self.outdir=''
        self.out=[]
        #print(self.indir)
        
#format directory
    def format_dir(self):
        #judge if absolute dir
        if self.indir.find('/')==0:
            pass
        elif self.indir.find('~')==0:
            self.indir = re.sub('~', os.getenv('HOME'), self.indir)
        elif self.indir is None:
            self.indir = os.getcwd()
        else: # ./test or test
            self.indir = os.path.abspath(self.indir)
        #alway followed by '/'
        if not self.indir[-1]=='/':
            self.indir += '/'
        return self.indir
#
    def create_dir(self, permission=0o755):
        #create directory
        if os.path.isdir(self.indir):
            print('{} exits! No action was done.'.format(self.indir))
        else:
            #
            try:
                os.makedirs(self.indir, permission)
                print('{} was created.'.format(self.indir))
            except OSError:
                self.indir=''
                print('Failed to creating {}.'.format(self.indir))
        return self.indir

    def delete_dir(self):
        try:
            os.rmdir(self.indir)
        except OSError:
            if not os.path.isdir(self.indir):
                raise

#
    def child_items(self):
        items=[]
        for entry in os.listdir(self.indir):
            if entry[0] != '.' and entry != 'lost+found':
                items.append(entry)
        return items
        
#list all directories with a given directory
    def child_dirs(self, pattern=None):
        #get all files and dirs
        entries=self.child_items()
        #filter out including characters
        if pattern != None:
            entries = [x for x in entries if pattern in x]
        #
        for entry in entries:
            #print(entry)
            abs_entry = self.indir + entry + '/'
            #not hidden dirctories: /.
            if os.path.isdir(abs_entry):
                self.out.append(abs_entry)
                #print(abs_entry)
        return self.out

#list all files with a given directory
    def child_files(self, suffix=None):
        entries=self.child_items()
        #filter out including characters
        if suffix != None:
            entries = [x for x in entries if x.endswith(suffix)]
        #
        for entry in entries:
            abs_entry = self.indir + entry
            #not hidden dirctories: /.
            if os.path.isfile(abs_entry):
                self.out.append(abs_entry)
                #print(abs_entry)
        return self.out

#list all files with a given directory and sub directories
#get all files
    def recrusive_files(self): 
        for root, dirs, files in os.walk(self.indir):
            #print root, files[:4]
            for filename in files:
                out_file = os.path.join(root, filename)
                if os.path.isfile(out_file) and out_file.find('/.') == -1:
                    #print(out_file)
                    self.out.append(out_file)
        return self.out
#
    def recrusive_suffix(self, suffix): 
        for root, dirs, files in os.walk(self.indir):
            #print root, files[:4]
            for filename in files:
                out_file = os.path.join(root, filename)
                if os.path.isfile(out_file) and out_file.find('/.') == -1:
                    if re.search(r''+suffix+r'$', out_file):
                        #print(out_file)
                        self.out.append(out_file)
        return self.out

#get directory from stdin:
    def stdin_dir(self, words='Enter directory'):
        print('{} (default:{}):.'.format(words, self.indir))
        stdin_dir = sys.stdin.readline()
        stdin_dir = stdin_dir.strip()
        if stdin_dir == '':
            self.outdir = self.indir
        else:
            self.outdir = dir_os(stdin_dir).create_dir()
        #print(self.outdir)
        return self.outdir
        
######################################################
#suggest infile with its absolute path
class file_os(dir_os):
    def __init__(self, infile, sep='='):
        self.file = str(infile)
        self.sep = sep
        self.outfile=''
        self.outdict={}
        
#'/home/yuan/test.txt': dir_name /home/yuan/
    def dir_name(self):
        #default is with absolute directory
        prefix = os.path.dirname(self.file)
        #print(prefix)
        #only file name, default as current folder
        if prefix.startswith('/'):
            file_dir = prefix
        elif prefix.startswith('.'):
            file_dir = re.sub('\\.', os.getcwd(), prefix, 1)
        elif prefix.startswith('~'):
            file_dir = re.sub('~', os.getenv('HOME'), prefix, 1)
        elif prefix == '':
            file_dir = os.getcwd()
        else:
            file_dir = os.getcwd() + '/' + prefix
        file_dir += '/'
        #print(file_dir)
        return file_dir
    
#'/home/yuan/test.txt': file_name is test.txt
    def file_name(self):
        #infile should be with absolute path
        file_name = os.path.basename(self.file)
        #print(file_name)
        return file_name
    
#'/home/yuan/test.txt': name_prefix is test
    def name_prefix(self):
        #file name
        name = self.file_name()
        name_prefix = name[:name.rfind(".")] if name.rfind(".") > 0 else name
        #print(name_prefix)
        return name_prefix

#'/home/yuan/test.txt': name_suffix is txt
    def name_suffix(self):
        #file name
        name = self.file_name()
        name_suffix = name[name.rfind(".")+1:] if name.rfind(".") > 0 else ''
        #print(name_suffix)
        return name_suffix

#'/home/yuan/test.txt': file_prefix is /home/yuan/test
    def file_prefix(self):
        file_prefix = self.dir_name() + self.name_prefix()
        #print(file_prefix)
        return file_prefix
        
    def change_suffix(self, name_suffix, file_name=None):
        name = self.file_name() if file_name is None else file_name
        
        #change
        new_name = name[:name.rindex('.')] if name.find('.') > 0 else name
        new_file_name = new_name + '.' + name_suffix
        #print(new_file_name)
        return new_file_name
        
    def file_size(self, unit='bytes'):
        #size by byte
        file_size = os.stat(self.file).st_size #bytes
        if unit == 'KB':
            file_size = file_size/1024
        elif unit == 'GB':
            file_size = file_size/(1024**2)
        elif unit == 'TB':
            file_size = file_size/(1024**3)
        #file_size='{:.2f} {}'.format(file_size, unit)
        file_size=round(file_size, 2)
        #print(file_size)
        return(file_size)

    def copy(self, out_dir):
        self.outdir=dir_os(out_dir).create_dir()
        self.outfile = self.outdir + self.file_name()
        if os.path.isfile(self.outfile):
            print(self.outfile, 'will be covered.')
        else:
            print(self.outfile,'will be created.')
        #copy file
        shutil.copy(self.file, self.outdir)
        return self.outfile

#no copy if file exists
    def soft_cp(self, out_dir):
        self.outdir=dir_os(out_dir).create_dir()
        self.outfile = self.outdir + self.file_name()
        if os.path.isfile(self.outfile):
            print('NO file is copied.')
        else:
            shutil.copy(self.file, self.outdir)
            return self.outfile
        
    def move(self, out_dir):
        #make directory
        self.outdir = dir_os(out_dir).create_dir()
        #move file
        self.outfile = self.outdir + self.file_name()
        #rename and move file
        self.rename(self.outfile)
        return self.outfile
        
    def rename(self, new_file):
        try:
            os.rename(self.file, new_file)
        except OSError:
            raise
        return new_file
        
    def download(self, out_dir, cover=False):
        self.outfile = dir_os(out_dir).create_dir() + self.file_name()
        #
        if os.path.isfile(self.outfile) and cover is False:
            print('Skip this download because {} exists in the local.'.format(self.outfile))
        else:
            web_obj = urllib.URLopener()
            web_obj.retrieve(self.file, self.outfile)
            print('Download {}, and store it as {}.'.format(self.file, self.outfile))
        return self.outfile
        
    #lines unit is million, default lines=4 million lines
    def split_txt(self, lines=4e6, out_dir=None):
        #out_dir
        self.outdir = self.dir_name() if out_dir is None else dir_os(out_dir).create_dir()
        #make temporary dir
        tmp_dir = dir_os(self.outdir+'tmp/').create_dir()
        prefix = tmp_dir + self.name_prefix() + '_'

        #print shell_command        
        command = 'split -l {} -d -a 3 {} {}'.format(int(lines), self.file, prefix)
        print('Split the file {}: {}'.format(self.file, command))
        output = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True).stdout.read()
        
        #move and rename files
        tmp_files = dir_os(tmp_dir).recrusive_files()
        for tmp_file in tmp_files:
            tmp_file_name = file_os(tmp_file).file_name()
            new_file = '{}{}.{}'.format(self.outdir, tmp_file_name, self.name_suffix())
            file_os(tmp_file).rename(new_file)
            print('\t', new_file)
        #remove tmp dir
        dir_os(tmp_dir).delete_dir()
        return output
    
    def readonly_handle(self):
        in_obj = gzip.open(self.file, 'rt') if self.name_suffix() == 'gz' else open(self.file, 'rt')
        return in_obj
    
    def write_handle(self):
        in_obj = gzip.open(self.file, 'wt') if self.name_suffix() == 'gz' else open(self.file, 'wt')
        return in_obj
    
    def to_dict(self):
        try: 
            in_obj = self.readonly_handle()
            for line in in_obj:
                line = line.strip()
                if 0 < line.find(self.sep) < len(line):
                    key, value = line.split(self.sep)
                    self.outdict[key] = value
                    #print key, value
            in_obj.close()
        except FileNotFoundError:
            print(self.file, 'didnot exit!')
            pass
        #print(self.outdict)
        return self.outdict
    
    def first_line(self):
        in_obj = self.readonly_handle()
        first_line = in_obj.readline()
        first_line = first_line.rstrip()
        in_obj.close()
        #print(first_line)
        return first_line
    
    #dict2 is nested dict
    def to_dict2(self):
        in_obj = self.readonly_handle()
        first_line = in_obj.readline()
        first_line = first_line.strip()
        header = first_line.split(self.sep)#colnames
        header.pop(0)#remove the first item
        for line in in_obj:
            line = line.strip()
            items = line.split(self.sep)
            rowname = items.pop(0)
            #print rowname
            dict2 = {}
            for index, value in enumerate(items):
                colname = header[index]
                dict2[colname] = value
            self.outdict[rowname] = dict2
            #print rowname, dict2, '\n'
        in_obj.close()
        #print(self.outdict)
        return self.outdict 

#dict2 is nested dictionary, but the file only two or three columns
#key1, key2, and value
#the types of key or value are string type
#not header
    def flat_file_to_dict2(self, value_index=2, mydict={}):
        in_obj = self.readonly_handle()
        for line in in_obj:
            line = line.strip()
            items = line.split(self.sep)
            key1, key2, value = items[0], items[1], items[value_index]
            if key2 not in  mydict:
                mydict[key1] = {}
            if not value:
                value = ""
            mydict[key1][key2] = value
        in_obj.close()
        return mydict 
        
#read flat file into a data frame
#at least three columns, and key1, key2, and value
#not header
    def flat_file_to_df(self, df_index=[0,1,2], value_numeric=True):
        #
        dict2 = {}
        #get index of rows, columns, and values
        row_index, col_index, value_index = df_index 
        #read file into nested dict
        in_obj=self.readonly_handle()
        for line in in_obj:
            line = line.strip()
            items = line.split(self.sep)
            row, col, value = items[row_index], items[col_index], items[value_index]
            if col not in dict2: dict2[col] = {}
            #assign value
            dict2[col][row] = float(value) if value_numeric == True else value
        in_obj.close()
        #convert dict2 to df and fill NA with 0
        mydf = pd.DataFrame(dict2).fillna(0)

        return mydf
    #
    def to_dbm(self):
        dbm_file = self.file+'.db'
        if not os.path.isfile(dbm_file):
            #generate db_file
            db_dict = anydbm.open(dbm_file, 'c')
            #read lines into dictionary
            db_dict = self.flat_file_to_dict2()
            db_dict.close()
        #tie db_dict with db_file
        db_dict = anydbm.open(dbm_file, 'c')
        return db_dict
        
    #read txt/csv file into data frame of pandas 
    def to_df(self, header=False, rowname=False):
        header = None if header is False else 'infer'
        rowname = None if rowname is False else 0
        try:
            df = pd.read_table(self.file, sep=self.sep, header=header, index_col=rowname)
        except ValueError:
            print ('Failed to opening ', self.file)
            return np.nan
        else:
            #print df
            print ('(nrow, ncol):', df.shape)
            return df
            
    #convert text file to xls file
    def to_xls(self, out_dir=None):
        self.outdir = self.dir_name() if out_dir is None else dir_os(out_dir).create_dir()
        self.outfile = self.outdir + self.name_prefix() + '.xls'
            
        #read file
        file_obj = self.readonly_handle()
        row_list = []
        for row in file_obj:
            row_list.append(row.split(self.sep))
        column_list = zip(*row_list)
        #set workbook and add sheet
        workbook = xlwt.Workbook()
        worksheet = workbook.add_sheet('Sheet1')
        i = 0 
        for column in column_list:
            for item in range(column.__len__()):
                worksheet.write(item, i, column[item])
            workbook.save(self.outfile)
            i += 1
        return(self.outfile)
        
    def line_replace(self, new_dict={}):
        #get old variables
        var_dict = self.to_dict()
        #refresh new
        for name in new_dict.keys():
            var_dict[name] = new_dict[name]
        #export to file
        myDict.basic(var_dict).dict_to_file(self.file, self.sep)
        
    def line_add(self, new_dict={}):
        #get old variables
        var_dict = self.to_dict()
        #refresh new
        for name in new_dict.keys():
            if name in var_dict:
                var_dict[name] = int(var_dict[name]) + int(new_dict[name])
            else:
                var_dict[name] = new_dict[name]
        #export to file
        myDict.basic(var_dict).dict_to_file(self.file, self.sep)
    
##
#end