# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 10:00:42 2016


@author: yuan
"""
#
import glob
import multiprocessing as mpd
import os
import subprocess
import sys
import time
#
import myIO


def phip_thread(file_var):
    #get current _dir   
    current_dir = os.path.dirname(os.path.abspath(__file__))
    script = current_dir + '/bioPHIPseq.py'
    #print script   
    #
    shell_command = 'python {} {} > {}.log'.format(script, file_var, file_var)
    print("@@@@@{}@@@@@@@\n{}\n".format(time.ctime(), shell_command))
    subprocess.Popen(shell_command, stdout=subprocess.PIPE, shell=True).stdout.read()
####
    
    
#main program
if __name__=="__main__":
    try:
        file_type = sys.argv[1]   #should be human or virus
    except IndexError:
        file_type = 'human'
    #############
    #revise 
    par = {'file_type':file_type, 'threads_num':24, 
         #'specieZ_threshold':10,'align_score':80,'sim_threshold':0.8,
         'phip_alignment':'yes', 'phip_counting':'yes', 'phip_merge':'yes', 
         'phip_GP':'yes', 'phip_zscores':'yes', 'phip_enrichment':'yes'}
    #input dir storing results 
    dirs = ['/home/yuan/results_phip', '/home-4/tyuan10@jhu.edu/work/yuan/results_phip']
    #get all variables.txt
    files_var = []
    for d in dirs:
        if os.path.isdir(d):
            files_formula = os.path.join(d, '*'+file_type+'*/variables.txt')
            #print files_formula
            files_var += glob.glob(files_formula)
            #sub=myIO.dir_os(d).recrusive_files('variables.txt')
            #for s in sub:
            #    if file_type in s: files_var.append(s)

    #get all command lines
    for index, file_var in enumerate(files_var):
        print(index+1, file_var)
        #revise the parameters of variables.txt
        myIO.file_os(file_var, '=').line_replace(par)
    
    #parallel processing
    #threads number
    #pool=mpd.Pool(processes=8)
    #pass one argument at a time 
    #pool.map(phip_thread, files_var)  
    #pool.close()
    #pool.join()    

    print('\n\n\n\nGreat! The batch running is done!\n\n\n')
#end