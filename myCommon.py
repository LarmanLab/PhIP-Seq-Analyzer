# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 09:36:58 2016

@author: yuan
"""
#standard modules
import collections
import multiprocessing as mp   # multi process
import multiprocessing.dummy as mpd #multi-threading
import os
import random
import numpy as np
import pandas as pd

#personal modules
import myDataframe
import myDict
import myIO
import myList
import myPlot




###############
class basic:
    def __init__(self, par=None):
        self.par = par

#multiple threads by threads pool
    def pp_map_threads(self, func, args_list):
        t = self.par['threads_num']
        print("Multiple threads = ", t)
        #t=1
        if t < 2:
            for par in args_list:
                func(par)
        else:
            #threads number
            pool = mpd.Pool(t)
            #pass one argument at a time 
            pool.map(func, args_list) 
            pool.close()
            pool.join()  

#multiple threads with mixed functions input
    def pp_apply_threads(self, args_list):
        t = self.par['threads_num']
        if t < 2: 
            for args in args_list:
                func, args_tuple = args[0], args[1:]
                func(args_tuple)
        else:
            pool = mpd.Pool(processes=t)
            for args in args_list:
                func, args_tuple = args[0], args[1:]
                print(func, args_tuple)
                #func(args_tuple)
                pool.apply(func, args=(args_tuple,) )
            pool.close()
            pool.join()
        
#multiple process by process pool
#Note: pass a function - and not a method - to the worker processes
    def pp_map_process(self, func, args_list):
        t = self.par['threads_num']
        print("Multiple processes: ", t)
        #t=1
        if t < 2:
            for par in args_list:
                func(par)
        else:
            #process number
            pool = mp.Pool(processes=t)
            #pass one argument at a time 
            pool.map(func, args_list)
            pool.close()
            pool.join()  
    
#peptides that a protein belongs to
#three parameters from annotation file required:
#pro_id, pep_id, pep_rank
    def protein_peptides(self):
        pro_pep = {}
        if os.path.isfile(self.par['file_pro_pep']):
            ##read pro_pep file
            pro_pep = myIO.file_os(self.par['file_pro_pep'], "\t").to_dict()
        else:
            #read annotation file
            annot_dict = myIO.file_os(self.par['file_annotation'], "\t").to_dict2()
            if 'Rnl2_SPIKEIN' in annot_dict: annot_dict['Rnl2_SPIKEIN']['pep_rank'] = 0
            in_pro = [ annot_dict[p][self.par['protein_assoc']] for p in annot_dict.keys()]
            in_pro = list(set(in_pro))
            print('In proteins:{}, In peptides:{}'.format(in_pro.__len__(), annot_dict.keys().__len__()))
            
            ##
            pro_rank_pep = {}
            for pep_id in annot_dict.keys():
                pro_id = annot_dict[pep_id][self.par['protein_assoc']]
                pep_rank = int(annot_dict[pep_id]['pep_rank'])
                if pro_id in pro_rank_pep:
                    pro_rank_pep[pro_id][pep_id] = pep_rank
                else:
                    pro_rank_pep[pro_id] = {pep_id:pep_rank}
                    #print pro_rank_pep[pro_id]
            #
            pep_num = 0
            for pro_id, pep_dict in pro_rank_pep.items():
                #print sorted(pep_dict.keys())
                peps = sorted(pro_rank_pep[pro_id], key = pro_rank_pep[pro_id].get)
                pep_num += len(peps)
                pro_pep[pro_id] = ','.join(peps)
            #export
            print("Number of protein:{}\tNumber of peptides:{}.".format(len(pro_pep.keys()), pep_num))
            myDict.basic(pro_pep).dict_to_file(self.par['file_pro_pep'], "\t")
        #
        return pro_pep
    
#dep_pep is nested dictionary, values are the list of pep_ids
    def taxon_dependent_peptides(self):
        dependent_pep = {}
        file_dependent = self.par['dir_ref_seq'] + 'virus_dependent_peptides.csv'
        if not os.path.isfile(file_dependent):
            print("Error: Failled to detect ", file_dependent)
            sys.exit(2)
        else:
            in_obj = open(file_dependent, 'rt')
            for line in in_obj:
                #line = line.strip("\n")
                items = line.split(',')
                pep1 = str(items[1])
                pep2 = str(items[3])
                if pep1 in dependent_pep:
                    dependent_pep[pep1].append(pep2)
                else:
                    dependent_pep[pep1] = [pep2]
            in_obj.close()
            #
            #print(dependent_pep['54664'])
            #print(dependent_pep['53343'])
        return dependent_pep
    
#collpase peptide matrix into protein matrix
    def collapse_matrix(self, pars):
        infile, outfile, collapse_func=pars
        print('t:', outfile)
        #read counts_file: pep_df
        file_sep = '\t' if myIO.file_os(infile).name_suffix() == 'txt' else ','
        pep_df = pd.read_csv(infile, header=0, index_col=0, sep=file_sep, low_memory=False)
        #pep_df.index=pep_df.index.astype(str)
        sample_names = list(pep_df)
        #print(sample_names)
        #both column and row names should be string type
        pep_df.columns = pep_df.columns.astype(str)
        pep_df.index = pep_df.index.astype(str)
        #combined pep_df with protein annotation
        pep_pro = self.par['annot_df'][['pep_id', self.par['protein_assoc']]]
        #pep_pro['pep_id']=pep_pro['pep_id'].astype(str)
        combined_df = pd.merge(pep_df, pep_pro, how='inner', left_index=True, right_on='pep_id')
        #group by protein id
        group_dict = combined_df.groupby([self.par['protein_assoc']], as_index=False).groups
        collapse = dict()
        for protein_id, row_names in group_dict.items():
            subdf = pep_df.ix[row_names]
            collapse[protein_id] = subdf.apply(collapse_func, axis=0)
            #if protein_id =='A0A126':
                #print subdf[['CTLA4.BEADS_ONLY.BEADS_ONLY.BEADS_ONLY.20A20G.1', 'CTLA4.BEADS_ONLY.BEADS_ONLY.BEADS_ONLY.20A20G.2']]
        #convert to data frame and transpose
        cdf = pd.DataFrame(collapse).transpose()
        #remove the first column
        cdf = np.round(cdf[sample_names], 1)
        #print cdf
        #export
        file_sep = '\t' if myIO.file_os(outfile).name_suffix() == 'txt' else ','
        cdf.to_csv(outfile, index_label=self.par['protein_assoc'], sep=file_sep)
        return cdf
    
#merge RC files        
    def QC_statistics(self):
        print("###Quality control: statistics summary")
        #print(self.par['sample_names'])
        #print(self.par['dir_result'])
        stat_dict = collections.defaultdict(dict)
        for sample_name in self.par['sample_names']:
            sample_log = '{}{}/{}.log'.format(self.par['dir_result'], sample_name, sample_name)
            stat_dict[sample_name] = myIO.file_os(sample_log, '=').to_dict()
        #convert to data frame
        stat_df = pd.DataFrame(stat_dict)
        stat_df = stat_df.transpose()
        
        #1: scatter plot 1
        sub_df = stat_df[['raw_reads_num', 'unique_aligned_reads_num']].astype(float)/1e6
        #print sub_df
        plot_par={'df':sub_df, 'title':'raw_reads_vs_aligned_reads', 
                  'picfile':self.par['dir_QC'] + 'raw_reads_vs_aligned_reads.png',
                  'pch':'o', 'text':'million reads'}
        myPlot.plot(plot_par).dotP()
        #2: scatter plot 2
        stat_df['unique_aligned_percentage'] = sub_df['unique_aligned_reads_num']*100/sub_df['raw_reads_num']
        plot_par['df'] = stat_df[['raw_reads_num','unique_aligned_percentage']].astype(float)
        plot_par['title'] = 'percentage_aligned_reads'
        plot_par['picfile'] = self.par['dir_QC'] + 'percentage_aligned_reads.png'
        myPlot.plot(plot_par).dotP()
        #3: export to csv file
        print('\tSave statistical summary into {}.'.format(self.par['file_stat']))
        stat_df.to_csv(self.par['file_stat'], index_label='sample_names')
        #
        
#saturation analysis    
    def QC_saturation(self):
        print("###saturation analysis\n")
        combined_df = {} 
        combined_dynamics = {}
        #plot suaturation curve per sample
        #n=1
        for sample_name in self.par['sample_names']:
            file_head = '{}{}/'.format(self.par['dir_result'], sample_name)
            #read saturation file
            df = pd.read_table(file_head+'QC_saturation.txt', sep="\t", index_col=False)
            #print list(df)
            #print list(df.index)
            
            #saturation curves
            saturation_df = df[['row_name', '1', '5', '10']]
            #shrink dict
            shrinked_index = myList.basic( list(saturation_df.index) ).interval_list()
            #print shrinked_index
            sample_df = saturation_df.ix[shrinked_index] #select rows
            #sample_df=sample_df.transpose().astype(float)
            sample_df.ix[:,0] = sample_df.ix[:,0]/1e6
            #print sample_df
            #scatter plot
            plot_par={'df':sample_df, 'legend':'upper left', 
                    'title': 'Saturation analysis (Sequencing depth)',
                      'picfile': file_head+'QC_saturation_analysis.png',
                      'xlabel':'Number of raw reads (million)',
                      'ylabel':'Number of references'}
            myPlot.plot(plot_par).lineP()
            #combine data frame
            sample_df.index = range(sample_df.shape[0])
            for cutoff in ['1','5','10']:
                sub_df = sample_df[['row_name',cutoff]].copy()
                sub_df.columns = ['raw_reads:'+sample_name, sample_name]
                if cutoff in combined_df:
                    combined_df[cutoff] = pd.merge(combined_df[cutoff], sub_df, 
                        left_index=True, right_index=True, how='outer')
                else:
                    combined_df[cutoff] = sub_df.copy()
            
            #dynamics analysis
            dynamics_df = df[['row_name', 'max']] #select df
            #shrink dict
            shrinked_index = myList.basic( list(dynamics_df.index) ).interval_list()
            sample_df = dynamics_df.ix[shrinked_index]#select rows
            sample_df.ix[:,0] = sample_df.ix[:,0]/1e6 #divided by millions
            sample_df.reset_index(drop=True, inplace=True)
            #combined
            combined_dynamics[sample_name] = sample_df
            #plot
            plot_par={'df':sample_df, 'legend':'upper left', 
                      'title': 'Saturation analysis:dynamics of read conts',
                      'picfile': file_head+'QC_read_counts_dynamics.png',
                       'xlabel': 'Number of raw reads (million)',
                       'ylabel':'Maximum read counts'}
            myPlot.plot(plot_par).lineP()
        #export saturated curves
        for cutoff in ['1','5','10']:
            plot_par={'df': combined_df[cutoff], 'legend':None,
                      'title': 'samples={}, RC-cutoff={}'.format(len(self.par['sample_names']), cutoff),
                      'picfile': '{}saturation_cuttoff_{}.png'.format(self.par['dir_QC'], cutoff),
                      'xlabel':'Number of raw reads (million)', 'ylabel':'Number of references'}
            myPlot.plot(plot_par).lineP(x_value=1)

        #export dynamics curves
        combined_dynamics = pd.concat(combined_dynamics, axis=1)
        combined_dynamics.columns = [':'.join(x) for x in list(combined_dynamics)]
        #print combined_dynamics.shape
        #print combined_dynamics
        plot_par={'df':combined_dynamics, 'legend':None,
                  'title': 'Sequencing depth,sample={}'.format(len(self.par['sample_names'])),
                  'picfile': '{}saturation_dynamics.png'.format(self.par['dir_QC']),
                  'xlabel':'Number of raw reads (million)', 'ylabel':'Maximum read counts' }
        myPlot.plot(plot_par).lineP(x_value=1)

#relationship between significant hits and raw read num
    def QC_hits(self, infile, threshold=None):
        print('###Relationship between significant hits and raw read num of ', infile)
        file_prefix = '{}{}_'.format(self.par['dir_QC'], myIO.file_os(infile).name_prefix())
        if threshold is None: threshold = float(self.par['zscore_threshold'])
        #read statistics file
        stat_df = pd.read_table(self.par['file_stat'], sep=",", index_col=0, low_memory=False)
        stat_df.index = stat_df['sample_name'] #assign row names
        stat_df = stat_df.ix[self.par['sample_names']]#order rows by sample_names
        raw_reads = stat_df['raw_reads_num']/1e6
        #print stat_df[['sample_name','raw_reads_num']]
        #read values file
        in_df = pd.read_table(infile, sep="\t", index_col=0, low_memory=False)#rownames and colnames
        order_df = in_df[self.par['sample_names']].copy()#order columns
        #print(order_df.shape)
        
        #plot of raw reads vs number of hits
        #print list(order_df)
        def func1(x,y=threshold):
            sig = x[x>=y]
            return len(sig)
        hits_num = order_df.apply(func1, axis=0)
        #get compared df
        comp_df = pd.DataFrame({'A': raw_reads, 'B': hits_num})
        comp_df.to_csv(file_prefix+'raw_vs_sighits.csv', sep=',')
        #plot
        plot_par={'df':comp_df, 'legend':None,
                  'title': 'Effects of sequencing depth on significant hits',
                  'picfile': file_prefix + 'raw_vs_sighits.png',
                  'xlabel':'Number of raw reads (million)',
                  'ylabel':'Number of signficant hits'}
        myPlot.plot(plot_par).dotP()
        
        #plot of raw reads vs mean values of hits
        #print list(order_df)
        def func2(x,y=threshold):
            x = pd.Series(x)
            #print list(x)
            sig = x[x>=y]
            #print list(sig)
            sig_mean = np.mean(sig)
            return sig_mean
        hits_mean = order_df.apply(func2, axis=0)
        #print hits_mean
        #get compared df
        comp_df = pd.DataFrame({'A': raw_reads, 'B': hits_mean})
        outfile=file_prefix+'raw_vs_mean_significant_hits.csv'
        print('\texport QC to {}.'.format(outfile))
        comp_df.to_csv(outfile, sep=',')
        #plot
        plot_par={'df':comp_df, 'legend':None,
                  'title': 'Effects of sequencing depth on significant hits',
                  'picfile': file_prefix + 'raw_vs_mean_significant_hits.png',
                  'xlabel':'Number of raw reads (million)',
                  'ylabel':'Mean values of signficant hits'}
        myPlot.plot(plot_par).dotP()
        
#combine value matrix and annotation matrix
    def combine_df(self, counts_file, annot_index='pep_id'):
        #read count file
        file_sep = '\t' if myIO.file_os(counts_file).name_suffix() == 'txt' else ','
        counts_df = pd.read_table(counts_file, sep=file_sep, index_col=0, low_memory=False)
        counts_df.index = [str(x) for x in counts_df.index]
        #print 'counts:', counts_df.shape
        #print list(counts_df.index)[:20]
        
        #read annotation file
        file_sep = '\t' if myIO.file_os(self.par['file_annotation']).name_suffix() == 'txt' else ','
        annot_df = pd.read_table(self.par['file_annotation'], sep=file_sep, index_col=None, low_memory=False)
        annot_df.index = [str(x) for x in annot_df[annot_index]]
        #print 'annot:', annot_df.shape
        #print list(annot_df[annot_index])[:20]

        #combine by rows
        comb_df = pd.merge(annot_df,counts_df, left_index=True, right_index=True, how='inner')
        comb_df.index = list(comb_df[annot_index])
        #comb_df=comb_df.rename(columns={self.par['protein_assoc']:'pro_id'})
        #print comb_df[['pep_id','row_name']]
        #print comb_df.shape
        #sample df
        sample_df = comb_df[self.par['sample_names']]
        sample_df.index = list(comb_df[annot_index])
        return (comb_df, sample_df)

#extract annotation from annotation file, return dict
#left_column should be unique
    def extract_annot(self, left_column, right_column, FUN):
        annot_dict = {}
        #read annotation file
        annot_df = myIO.file_os(self.par['file_annotation'], sep="\t").to_df(header=True, rowname=False)
        for index, row in annot_df.iterrows():
            key = row[left_column]
            value = FUN(row[right_column])
            annot_dict[key] = value
            #print "%s:%s" % (key, value)
        return annot_dict

#permutation of key, and return the permuatated matrix of value 
#Note: key vs multiple values, default sep is comma
    def hits_permutation1(self, in_dict, sample_size=10):
        #get the pool for sampling
        pool = in_dict.keys()
        #
        permute_dict = {}
        for i in range(self.par['permutation_times']):
            #random select some keys from the pool
            random_keys = random.sample(pool,sample_size)
            random_values = {}
            for k in random_keys:
                values_list = in_dict[k].split(',')
                for v in values_list:
                    if v in random_values:
                        random_values[v] += 1
                    else:
                        random_values[v] = 1
            #
            permute_dict[i] = random_values
        #transform dict: times in columns, value of in_dict is in rows
        permute_dict = myDict.basic(permute_dict).transform_dict2()
        return permute_dict

#permutation of key, and return the permuatated times with a given target
#Note: key vs multiple values, default sep is comma
    def hits_permutation2(self, in_dict, target, sample_size=10):
        #get the pool for sampling
        pool = in_dict.keys()
        #
        permute_nums = []
        for i in range(self.par['permutation_times']):
            num = 0
            #random select some keys from the pool
            random_keys = random.sample(pool,sample_size)
            for k in random_keys:
                values_list=in_dict[k].split(',')
                for v in values_list:
                    if v == target: num += 1
            permute_nums.append(num)
        return permute_nums

#input is pep_id list, return pro_peps dict: pro_id vs peps list
    #def pep_to_pro_dict(self, pep_pro_dict, pep_list):

#permutation: select m elements from n elements (n>m) with x times
#sampling pool include n elements as list type
#sampling_type might be transcirpt_id or unipro_acc
    def permute(self, sampling_pool, sampling_number, sampling_type):
        file_permutation = self.par['dir_permutation'] + sampling_type + str(self.par['permutation_times']) + sampling_number
        #permutate data if no such file
        if not os.path.isfile(file_permutation):
            df = pd.DataFrame()
            for i in range(self.par['permutation_times']): 
                df[i] = random.sample(sampling_pool, sampling_number)            
            #export
            df.to_csv(file_permutation, sep='\t', index=False, header=False)
        #or directly read samplings
        permuted_df = pd.read_csv(file_permutation, sep="\t", header=None, index_col=None, lower_memory=False)
        return permuted_df

#permutation of specie alignment
    def permute_taxon_blast(self, hits_num):
        print('permutation of viral blast:{}\t{}'.format(self.par['type'], hits_num))
        #
        counts_df = pd.DataFrame()
        outfile = '{}{}.txt'.format(myIO.dir_os(self.par['dir_out']).create_dir(), hits_num)
        if os.path.isfile(outfile):
            print('Read file: ', outfile)
            counts_df = pd.read_csv(outfile, header=0, index_col=0, sep="\t", low_memory=False)  
        else:
            #1: permutated peptides
            pep_names = list(self.par['binary_aln_df'].index)
            pep_df = myList.basic(pep_names).permute_list(self.par['permutation_times'], hits_num)
            #2: permutation based on the non-overlapped hits num
            for col, perm_pep in pep_df.items():
                perm_zb = self.par['binary_aln_df'].ix[perm_pep]
                p_collapse_zb, p_sim_tag = myDataframe.basic(perm_zb).unispecie(self.par['sim_threshold'])
                counts_df[col] = p_collapse_zb.apply(sum,axis=0) + p_sim_tag
                #print list(perm_tmp[col])
            #export
            counts_df.to_csv(outfile, sep='\t', header=True, index_label=self.par['type'])
        #combine permuated counts
        #print counts_df.shape
        perm_mean = counts_df.apply(lambda x: np.mean(np.floor(x)), axis=1).round()
        #print perm_mean
        return perm_mean
        
#enrichment of protein motifs: infile is usually zscores_file of peptides
    def enrich_pro(self, infile, annot_A, annot_B, sep1, sep2):
        if annot_A is None: annot_A='transcript_id'
        if annot_B is None: annot_B='pro_motifs'
        print("Enrichment analysis of {} => {} : {}".format(annot_A, annot_B, infile))
        #read data frame
        file_sep = ',' if infile.endswith('.csv') else '\t'
        counts_df = pd.read_csv(infile, index_col=0, sep=file_sep, low_memory=False)
        #get all ids connect counts_df with annot_df
        A_ids = list(self.par['annot_df'][annot_A])
        #get all ids based on annot_type in list formate
        B_ids = myDataframe.basic(self.par['annot_df']).df_list(annot_B, sep1, sep2)
        #get A_ids vs list of b_ids in dict formate
        AB_dict = myDataframe.basic(self.par['annot_df']).list_dict(annot_A, annot_B, sep1, sep2)
        
        #initiate: #frequency of observed enriched motifs
        hits_observed = myDict.basic().init_dict2(B_ids, list(counts_df), 0)
        #initiate: zscores of obs based on permutation models
        hits_zscores = myDict.basic().init_dict2(B_ids, list(counts_df), 0)
        #initiate: detect bugs
        debugging = myDict.basic().init_dict2(B_ids+['hits_counts','interact_counts'], {}, 'NA')
        #loop of data frame by columns
        for sample_name, zscores in counts_df.items():
            #print sample_name
            zscores = pd.Series(zscores)
            zscores.index = list(counts_df.index)
            #1: get ids of significant hits
            sig_zscores = zscores[zscores>=self.par['zscore_threshold']]
            obs_ids = list(sig_zscores.index)
            sig_num = len(obs_ids)
            #print annot_B, sample_name,sig_num
            #2: count frequency of enriched annotations, namely motifs
            obs_freq, obs_details = myDict.basic(AB_dict).elements_frequency(obs_ids)
            #print obs_freq.values()
            #debugging
            debugging['hits_counts'][sample_name] = sig_num
            debugging['interact_counts'][sample_name] = sum(obs_freq.values())
            
            #3: permute samples
            #print "\tenrichment: %s\t%s\t%s" % (sample_name, sig_num, len(obs_freq.keys()))
            perm_dict = {}
            for i in range(self.par['permutation_times']):
                perm_peps = random.sample(A_ids, sig_num)
                tmp_perm, tmp_details = myDict.basic(AB_dict).elements_frequency(perm_peps)  # frequency dict
                for key, value in tmp_perm.items():
                    if key in perm_dict:
                        perm_dict[key].append(value)
                    else:
                        perm_dict[key] = [value]
            #print perm_dict
           
            #4: calcuate z-scores of observed counts
            for enriched_id, obs_num in obs_freq.items():
                #update hit_observed
                hits_observed[enriched_id][sample_name] = obs_num #frequency of observed enriched annot
                #update debugging
                debugging[enriched_id][sample_name] = '{}:{}'.format(obs_num, obs_details[enriched_id])
                #update zscores_dict
                if enriched_id in perm_dict:
                    perm_pools = perm_dict[enriched_id]
                    #append zero and all pools are the same length
                    perm_pools = perm_pools + [0] * (5-len(perm_pools))
                    perm_mean = np.mean(perm_pools)
                    perm_sd = np.std(perm_pools)
                    #zscores of observed hits against the null model
                    zscore = (obs_num-perm_mean)/perm_sd if perm_sd > 0 else (obs_num-perm_mean)
                    hits_zscores[enriched_id][sample_name] = round(zscore, 2)
                else:
                    hits_zscores[enriched_id][sample_name] = obs_num
            #print hits_zscores
        
        #export
        file_head = '{}{}_{}_'.format(self.par['dir_enrichment'], myIO.file_os(infile).name_prefix(), annot_B)
        myDict.basic(hits_observed).dict2_to_file(out_file=file_head+'counting.txt', index_label=annot_B)
        myDict.basic(hits_zscores).dict2_to_file(out_file=file_head+'zscores.txt', index_label=annot_B)
        myDict.basic(debugging).dict2_to_file(out_file=file_head+'debugging.txt', index_label=annot_B, NA='NA')

    #combine counting files into a matrix
    def f(self, args_tuple):
        #row_names should be None or list type
        infile_tail, RC_level, out_file, row_names = args_tuple
        #
        counting_dict2 = {}
        for sample_name in self.par['sample_names']:
            #get read counts of a given sample
            counting_file = self.par['dir_result'] + sample_name + '/' + sample_name + infile_tail
            sample_dict2 = myIO.file_os(counting_file, '\t').to_dict2()
            for ref in sample_dict2.keys():
                #print ref
                counts = sample_dict2[ref][RC_level]
                if ref in counting_dict2:
                    counting_dict2[ref].update({sample_name:counts})
                    #print '=='+ref+'=='
                else:
                    counting_dict2[ref] = {sample_name:counts}
                #print sample_name, ref,counting_dict2[ref]
        #export counting_dict
        myDict.basic(counting_dict2).dict2_to_file(out_file=out_file, row_names=row_names)
        #return RC_files
#end            
