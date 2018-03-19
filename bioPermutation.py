# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 15:03:36 2017

@author: yuan
"""

import os
import sys
import time
import multiprocessing.dummy as mpd

#personal modules
import myCommon
import myDataframe
import myIO
import mySystem




#main program
if __name__ == "__main__":
    start_time = time.time()
    dir_bin,script_name = os.path.split(os.path.realpath(__file__))
    dir_home = os.path.dirname(dir_bin)
    
    #pass arguments
    start, end = sys.argv[1].split('-')
    par = {'specie_permutation':'yes', 'organism_permutation':'yes', 'threads_num':24,
         'start':int(start), 'end':int(end)+1, 'align_score':80, 'sim_threshold':0.8,
         'dir_bin':dir_bin+'/', 'dir_home':dir_home+'/', 'permutation_times':100}
    par['dir_permutation'] = myIO.dir_os(par['dir_home']+'permutation/').create_dir()
    
    
    print('###permutation procedure\n\n')
    pool = mpd.Pool(processes=par['threads_num'])

    #permuation of organism alignment
    if par['organism_permutation'] == 'yes':
        #read aln file
        file_aln = par['dir_home']+'ref_seq/organism_blast.txt'
        par['binary_aln_df'] = myDataframe.basic().aln_df(file_aln, par['align_score'])
        par['type'] = myIO.file_os(file_aln).name_prefix()
        par['dir_out'] = myIO.dir_os(par['dir_home']+'permutation/'+par['type']).create_dir()
        #
        for hits_num in range(par['start'], par['end']):
            pool.apply_async(myCommon.basic(par).permute_taxon_blast, args=(hits_num,))
            time.sleep(1)
            
    #permuation of specie alignment
    if par['specie_permutation'] == 'yes':
        #read aln file
        file_aln = par['dir_home']+'ref_seq/specie_blast.txt'
        par['binary_aln_df'] = myDataframe.basic().aln_df(file_aln, par['align_score'])
        par['type'] = myIO.file_os(file_aln).name_prefix()
        par['dir_out'] = myIO.dir_os(par['dir_home']+'permutation/'+par['type']).create_dir()
        #
        for hits_num in range(par['start'], par['end']):
            pool.apply_async(myCommon.basic(par).permute_taxon_blast, args=(hits_num,))
            time.sleep(1)
        

            
    #return process
    pool.close()
    pool.join()  


    #end
    times_dict = mySystem.system().get_time(start_time)
    print('\n\nEnding at ', time.ctime(), 'Duration is ', times_dict['duration'])
    print('\n\nGreat! It is done!!!\n\n\n')
    
