# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Tue Sep 22 07:47:14 2015

@author: yuan
"""
import os
import sys
import time
#personal modules
import myAlign
import myCommon
import myDataframe
import myDict
import myGenome
import myIO
import myParallel
import mySystem
import bioPHIPfunc

class launch_phip:
    def __init__(self, par):
        #default var is variables.txt
        self.par = par

#initiate default dir and file
    def init_dir_file(self):
        self.par['dir_home'] = myIO.dir_os(self.par['dir_home']).create_dir()
        print('home directory of phip tool:', self.par['dir_home'])
        #dir_home = /home/yuan/phip/
        
        #alignment related
        self.par['dir_aligner'] = self.par['dir_home']+'bowtie1/'
        self.par['aligner_options'] = '{}bowtie {}'.format(self.par['dir_aligner'], self.par['aligner_options'])
        self.par['genome_index'] = self.par['dir_aligner'] + self.par['genome_index_name']
        self.par['dir_ref_seq'] = self.par['dir_home']+'ref_seq/'
        self.par['file_ref_fa'] = '{}{}.fa'.format(self.par['dir_ref_seq'], self.par['genome_index_name'])
        if 'file_annotation' in self.par.keys():
            self.par['file_annotation'] = self.par['dir_ref_seq'] + self.par['file_annotation']
        #
        #judge ref library human or virus
        if 'VirScan' in self.par['genome_index_name']:
            self.par['lib'] = 'virus'
            self.par['file_NC'] = self.par['dir_ref_seq'] + 'virus_BeadsOnly.txt'
        elif 'human' in self.par['genome_index_name']:
            self.par['lib'] = 'human'
            self.par['file_NC'] = self.par['dir_ref_seq'] + 'human_BeadsOnly.txt'
        elif 'PublicEpitope' in self.par['genome_index_name']:
            self.par['lib'] = 'PE'
        elif 'LISH' in self.par['genome_index_name']:
            self.par['lib'] = 'LISH'
                
        #dir of raw data
        if 'dir_raw_data' not in self.par.keys():
            self.par['dir_raw_data'] = myIO.dir_os(self.par['dir_home']+'raw_data').create_dir()
        #results related
        if 'dir_result' not in self.par.keys():
            self.par['dir_result'] = myIO.dir_os(self.par['dir_home']+'result').create_dir()
        #print('Result directory', self.par['dir_result'])
        if 'dir_result_array' not in self.par.keys():
            self.par['dir_result_array'] = self.par['dir_result']
        
        #dir of statistics
        self.par['dir_stat'] = myIO.dir_os(self.par['dir_result']+'statistics').create_dir()
        self.par['dir_QC'] = myIO.dir_os(self.par['dir_stat']+'QC').create_dir()
        self.par['dir_enrichment'] = myIO.dir_os(self.par['dir_stat']+'enrichment').create_dir()
        
        #sample info
        self.par['file_sample_info'] = self.par['dir_result']+'sample_info.csv'
        self.par['dir_log'] = self.par['dir_result'] + 'sample_log/'
        self.par['file_log'] = self.par['dir_result'] + 'output.log'
        self.par['file_total_log'] = self.par['dir_result'] + 'Total.log'
        self.par['file_stat'] = self.par['dir_QC'] + 'statistics.csv'
        self.par['file_ref_txt'] = self.par['dir_result'] + 'references.txt'
        self.par['file_pro_pep'] = self.par['dir_result'] + 'protein_peptides.txt'
        #raw data related
        #print(self.par['dir_raw_data'])
        #
        self.par['RC_levels']=['lowRC']#lowRC, midRC, highRC
        self.par['phip_levels'] = ['pep','promax','prosum']
        files_dict = {}
        for pl in self.par['phip_levels']:
            file_head = '{}{}_'.format(self.par['dir_stat'], pl)
            #raw reads
            files_dict[pl+'_RC'] = file_head+'RC.txt'
            #noramlized by total raw counts
            files_dict[pl+'_scalingRC'] = file_head+'scalingRC.txt'
            files_dict[pl+'_scalingRC_prosum'] = file_head+'scalingRC_prosum.txt'
            files_dict[pl+'_scalingRC_promax'] = file_head+'scalingRC_promax.txt'
            #scalingRC against regressed median of phip sample and regressed sd of negative controls
            files_dict[pl+'_NCPHIPzscores'] = file_head+'NCPHIPzscores.txt'
            files_dict[pl+'_NCPHIPzscores_prosum'] = file_head+'NCPHIPzscores_prosum.txt'
            files_dict[pl+'_NCPHIPzscores_promax'] = file_head+'NCPHIPzscores_promax.txt'
        self.par['files_dict'] = files_dict
        
        #default parameters
        self.par['specieZ_threshold'] = int(self.par['specieZ_threshold']) if 'specieZ_threshold' in self.par.keys() else 10
        self.par['align_score'] = float(self.par['align_score']) if 'align_score' in self.par.keys() else 80
        #p value cutoff for binomial testing
        self.par['p_threshold'] = float(self.par['p_threshold']) if 'p_threshold' in self.par.keys() else .001
        #x value is observed successes cutoff for binomial test
        self.par['x_threshold'] = float(self.par['x_threshold']) if 'x_threshold' in self.par.keys() else 1
        self.par['sim_threshold'] = float(self.par['sim_threshold']) if 'sim_threshold' in self.par.keys() else 0.8
        self.par['zscore_threshold'] = int(self.par['zscore_threshold']) if 'zscore_threshold' in self.par.keys() else 10
        self.par['permutation_times'] = int(self.par['permutation_times']) if 'permutation_times' in self.par.keys() else 100
        self.par['threads_num'] = int(self.par['threads_num']) 
        self.par['scaling_factor'] = int(self.par['scaling_factor']) if 'scaling_factor' in self.par.keys() else 1e6
        
        #print self.par
        myDict.basic(self.par).print_dict()
        #
        return(self.par)

#initiate sample information, bowtie index, and genome annotations
#initiate analysis 
    def init_analysis(self):
        #1: read annotation file
        if 'file_annotation' in self.par.keys():
            self.par['annot_df'] = myDataframe.basic().annot_df(self.par['file_annotation'])
            #genome annotation: associations of protein-peptides 
            self.par['dict_pro_pep']=myCommon.basic(self.par).protein_peptides()
            #virus only
            if 'VirScan' in self.par['file_annotation']:
                #extract aa stretch
                #get dependent petides that two peptides shared at least  7-aa.
                self.par['dependent_pep']=myCommon.basic(self.par).taxon_dependent_peptides()
        
        #2: check bowtie or build bowtie index
        myAlign.alignment(self.par).build_bowtie_index()
        
        #3: sample info
        self.par = myParallel.samples(self.par).export_sample_info()
        #samples of negative controls
        group1 = self.par['group1']
        if 'NC' in group1.keys():
            self.par['NC_samples'] = group1['NC'].split(',')
            self.par['phip_samples'] = list(set(self.par['sample_names'])-set(self.par['NC_samples']))
            print('\nNumber of negative Controls (Beads only): ', self.par['NC_samples'].__len__())
            print('Number of PhIP samples: ', self.par['sample_names'].__len__())
            #myDict.basic(self.par['sample_dirs']).print_dict()
        
        #read reference sequence file (*.fa)
        ref_dict, ref_ids = myGenome.genome(self.par['file_ref_fa']).read_fa()
        self.par['ref_dict'] = ref_dict
        

####
    def init_par(self):
        self.init_dir_file()
        self.init_analysis()
        return self.par
    
#################################################################################
#main program
if __name__ == "__main__":
    start_time = time.time()
    
    #check python version
    if int(sys.version[0])<=2:
        print(sys.version)
        print('\nError: The version of python required for the pipeline running is at least v3.4.~\n')
        sys.exits(2)
    
    ########################################
    #read variables.txt
    #var_file = '/home/yuan/rawdata_phip/phipseq17_virus_variables.txt'
    var_file = os.path.abspath(sys.argv[1])
    #print(var_file)
    par = myIO.file_os(var_file, '=').to_dict()
    par['file_var'] = var_file
    #initiate parameters, directories and files
    par = launch_phip(par).init_par()
    
    ######################################
    #main loop
    bioPHIPfunc.phip(par).main_loop()
    ######################################
    
    #end
    times_dict = mySystem.system().get_time(start_time)
    myIO.file_os(par['file_total_log'], '=').line_replace(times_dict)
    print('\n\nDuration: ', times_dict['duration'])
    
    print('\n\nGreat! It is done!!!\n\n\n')
    
