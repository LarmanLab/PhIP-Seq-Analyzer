# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 07:46:46 2015

@author: yuan
"""
#standard modules
#import multiprocessing as mp
import multiprocessing.dummy as mpd
import time
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests
#personal modules
import myAlign
import myCommon
import myDataframe
import myDict
import myGenome
import myIO
import myList
import mySystem
import myStat

def mp_alignment(args_list):
    instance_phip, sample_name = args_list
    #print(instance_method, args_z)
    return instance_phip.phipseq_alignment(sample_name)
    
#class
class phip:
    def __init__(self, par):
        self.par = par
        #self.par['binary_aln_df'] = myDataframe.basic().aln_df(self.par['file_aln'],self.par['align_score'])

#main programm
    def main_loop(self):
        print("\n\n####Parameters of PHIP: \n")
        #parallel procesing
        if self.par['phip_alignment'] == 'yes' or self.par['phip_counting'] == 'yes':
            sample_names = self.par['sample_names']
            print(sample_names.__len__(), ' samples will be analyzed.\n')
            #multi-threads
            #myCommon.basic(self.par).pp_map_threads(self.phipseq_alignment, sample_names)
            #multi-processes
            myCommon.basic(self.par).pp_map_process(mp_alignment, [(self, s) for s in sample_names])
            
        #combine RC and statistics file
        if self.par['phip_merge'] == 'yes':
            pep_names = myGenome.genome(self.par['file_ref_fa']).read_fa()[1]
            #1: combine RC files into RC matrix
            print('\n\n\n###Combine RC files (phip_merge)\n')
            #get arguments
            args_list = []
            RC_level='lowRC'
            #peptide level: lowRC
            out_file = self.par['files_dict']['pep_RC']
            arg_tuple = ('_RC.txt', RC_level, out_file, pep_names) 
            args_list.append(arg_tuple)
            if 'file_annotation' in self.par:
                #promax level
                out_file = self.par['files_dict']['promax_RC']
                arg_tuple = ('_pro_maxRC.txt', RC_level, out_file, None)
                args_list.append(arg_tuple)
                #prosum level
                out_file = self.par['files_dict']['prosum_RC']
                arg_tuple = ('_pro_sumRC.txt', RC_level, out_file, None)
                args_list.append(arg_tuple)
            #multi-threads
            myCommon.basic(self.par).pp_map_threads(myAlign.alignment(self.par).combine_countfiles, args_list)
            #myCommon.basic(self.par).pp_apply_threads(args_list)
            #2: generate statistics.csv
            myCommon.basic(self.par).QC_statistics()
            
        #significance analysis using Z score
        if self.par['phip_zscores'] == 'yes':
            print('\n\n\n###normalization of RC (phip_zscores)\n')
            #peptides level
            RC_file = self.par['files_dict']['pep_RC'] #infile
            #1: scaling RCs
            sRC_file = self.par['files_dict']['pep_scalingRC']  # outfile
            myStat.normalization(self.par, RC_file, sRC_file, 'pep_id').RC_scaling()
            #2: z-scores of scaling RCs against negative controls and phipseq samples
            zfile = self.par['files_dict']['pep_NCPHIPzscores'] #outfile
            if 'file_NC' in self.par.keys():
                myStat.normalization(self.par, sRC_file, zfile, 'pep_id').NCPHIPzscores_PN()
            else:
                myStat.normalization(self.par, sRC_file, zfile, 'pep_id').NCPHIPzscores_RLM()
            
            #3:collpase matrix
            if 'file_annotation' in self.par:
                print("\t######collapse peptide matrix into protein matrix")
                pars = []
                for name in ['scalingRC', 'NCPHIPzscores']:
                    pep_file = self.par['files_dict']['pep_'+name] #infile
                    sum_file = self.par['files_dict']['pep_'+name+'_prosum']#outfile
                    pars.append((pep_file, sum_file, sum))
                    max_file = self.par['files_dict']['pep_'+name+'_promax']  #outfile
                    pars.append((pep_file, max_file, max))
                #multiple-threading
                myCommon.basic(self.par).pp_map_threads(myCommon.basic(self.par).collapse_matrix, pars)

        #Functional analysis after normalization and correction
        #parallel processing
        print('\n\n\n###Functional Analysis (phip_GP and phip_enrichment)\n')
        pool = mpd.Pool(processes = self.par['threads_num'])
        #set the list of parameters
        pep_zfile = self.par['files_dict']['pep_NCPHIPzscores'] #infile
        promax_zfile = self.par['files_dict']['pep_NCPHIPzscores_promax']
        prosum_zfile = self.par['files_dict']['pep_NCPHIPzscores_prosum']
        if self.par['phip_GP'] == 'yes':
            #1: polyclonal of signficant peptides
            pool.apply_async(self.sig_polyclonal, args = (pep_zfile,) )

            #virus only
            if 'VirScan' in self.par['file_annotation']:
                #5: inter/intra specie searching only for virus library#####
                pool.apply_async(self.taxon_spec, args=(pep_zfile, 'phip_taxon', 'pep_id',))
                #6: specie alignment of virus only
                file_aln = self.par['dir_ref_seq'] + 'specie_blast.txt'
                pool.apply_async(self.taxon_blast, args=(file_aln, pep_zfile,) )
                #7: organism alignment of virus only
                file_aln = self.par['dir_ref_seq'] + 'organism_blast.txt'
                pool.apply_async(self.taxon_blast, args=(file_aln, pep_zfile,) )

            ##quality control
            #1: relationship between significant hits and raw read num
            pool.apply_async(myCommon.basic(self.par).QC_hits, args=(pep_zfile,))
            pool.apply_async(myCommon.basic(self.par).QC_hits, args=(prosum_zfile,))
            pool.apply_async(myCommon.basic(self.par).QC_hits, args=(promax_zfile,))
            #2:saturation analysis
            pool.apply_async(myCommon.basic(self.par).QC_saturation)
            
        if self.par['phip_enrichment'] == 'yes':
            #5:Detection of enriched protein motifs
            E = myCommon.basic(self.par)
            if 'pro_motifs' in list(self.par['annot_df']):
                pool.apply_async(E.enrich_pro, args=(pep_zfile, 'pep_id', 'pro_motifs', ';', ',',))
            #6:GO,loci,PPI,KEGG,InterPro, multifunctional scaffold protein enrichment analysis
            terms = set(['GO','map','PPI','KEGG','InterPro','MIM', 'autoantigen']) & set(list(self.par['annot_df']))
            for term in terms:
                pro = self.par['protein_assoc']
                pool.apply_async(E.enrich_pro, args=(prosum_zfile, pro, term, ',', None,))
                pool.apply_async(E.enrich_pro, args=(promax_zfile, pro, term, ',', None,))
        pool.close()
        pool.join()
################################

    
    #step 1: parallel processing
    def phipseq_alignment(self,sample_name):
        print('\n######Anslysis of {} will be trigerred!#####'.format(sample_name))
        #initiate sample par
        sample_var = dict(self.par)
        sample_var['start_time'] = time.time()
        #sample name
        sample_var['sample_name'] = sample_name
        #sample directory
        sample_dir = self.par['sample_dirs'][sample_name]
        sample_var['sample_dir'] = myIO.dir_os(sample_dir).create_dir()
        print('\tSample directory: ', sample_var['sample_dir'])
        #raw data
        sample_var['sample_raw_files'] = ','.join(sample_var['sample_to_raw'][sample_name])
        print('\tRaw files: ', sample_var['sample_raw_files'])
        #export
        sample_var['file_head'] = sample_var['sample_dir'] + sample_name
        #default same file
        sample_var['sample_sam_file'] = sample_var['file_head'] + '.sam'
        #file of read counts
        sample_var['sample_RC_file'] = sample_var['file_head'] + '_RC.txt'
        sample_var['sample_pro_sumRC_file'] = sample_var['file_head'] + '_pro_sumRC.txt'
        sample_var['sample_pro_maxRC_file'] = sample_var['file_head'] + '_pro_maxRC.txt'
        #file for saturation analysis
        sample_var['sample_saturation_file'] = sample_var['file_head'] + '_saturation.txt'
        #sample log
        sample_var['sample_log'] = sample_var['file_head'] + '.log'
    
        
        #sequence alignment
        if sample_var['phip_alignment'] == 'yes':
            print("\n###sequence alignment", sample_var['tool_aligner'])
            #output is sam file
            if sample_var['tool_aligner'] == 'bowtie1':
                myAlign.alignment(sample_var).bowtie1_alignment()
                
        #counts reads    
        if sample_var['phip_counting'] == 'yes':
            #RC matrix by peptides  
            myAlign.alignment(sample_var).count_reads()
            #RC matrix by proteins
            if 'file_annotation' in self.par.keys():
                self.combine_peptides(sample_var)
        
        #update sample log
        sample_times = mySystem.system().get_time(sample_var['start_time'])
        sample_times['sample_name'] = sample_name
        myIO.file_os(sample_var['sample_log'], '=').line_replace(sample_times)
    
#get RC matrix by proteins
    def combine_peptides(self,sample_var):
        #read annotation file
        annot_df = pd.read_table(self.par['file_annotation'], index_col=False)
        sub_annot_df = annot_df[['pep_id', self.par['protein_assoc']]]
        #print(sub_annot_df)
        #read RC file
        RC_df = pd.read_table(sample_var['sample_RC_file'], index_col=False)
        #print(list(RC_df))
        
        #combine two data frame
        select_col = list(self.par['RC_levels'])
        select_col.insert(0, self.par['protein_assoc'])
        comb_df = pd.merge(sub_annot_df, RC_df, left_on='pep_id', right_on='row_name')
        comb_df = comb_df[select_col]
        #print(comb_df)
        #group by sample_par['protein_assoc']
        pro_sumRC_df = comb_df.groupby(self.par['protein_assoc'], axis=0).sum()
        #print(pro_sumRC_df)
        #export to file
        pro_sumRC_df.to_csv(sample_var['sample_pro_sumRC_file'], sep="\t")
        #group by sample_par['protein_assoc']
        pro_maxRC_df = comb_df.groupby(self.par['protein_assoc'], axis=0).max()
        #export to file
        pro_maxRC_df.to_csv(sample_var['sample_pro_maxRC_file'], sep="\t", index_label=self.par['protein_assoc'])
        
         
#virus taxonomy specifications
    def taxon_spec(self, count_file,taxon_rank,annot_index):
        #combine two data frame
        combined_df, phip_df = myCommon.basic(self.par).combine_df(count_file, annot_index)
        #print(combined_df)
        #print(list(combined_df))

        #taxonomy names: 
        taxon_group = combined_df.groupby(taxon_rank).groups
        taxon_names = taxon_group.keys()
        taxon_names = [t for t in taxon_names if str(t) != 'nan'] #remove nan
        #print(taxon_names)
        taxon_pairs = {'phip_specie':'InterSpecie', 'phip_genus':'InterGenus', \
                     'phip_family':'InterFamily', 'phip_taxon':'InterTaxon'}
        taxon_inter = taxon_pairs[taxon_rank]
        
        #inter-score dict
        #taxon_inter should be pep_ids separated by comma
        pepid_taxoninter = pd.Series(combined_df[taxon_inter], index = list(phip_df.index) )
        inter_df = myDataframe.basic(phip_df).interact_df(pepid_taxoninter, max, count_file+taxon_inter)

        #make permutation of pep_ids
        #permute_dict = myList.basic(list(phip_df.index)).permute_Series(self.par['permutation_times'], slice_dict = taxon_group)
            
        
        #the hits of significant specie specific
        #rows are peptides, and columns are phip samples plus species names
        #z-scores matrix of specific peptides
        #initiate nested dict
        taxon_dict = dict([(s,{}) for s in list(phip_df)]) # number of hits
        taxon_dict['peptides'] = dict([ (a, len(b)) for a, b in taxon_group.items()])
        #taxon_pval_dict = dict([(s,{}) for s in list(phip_df)]) #pvalues of the hits by permutations
        taxon_pep_dict = dict([(s,{}) for s in list(phip_df)]) #pepid and zscores of hits
        debugging_dict = {} #for identify bugs
        for s in list(phip_df):
            debugging_dict[s+':all_hits'] = {}
            debugging_dict[s+':inter_hits'] = {}
            debugging_dict[s+':intra_hits'] = {}
            debugging_dict[s+':hits'] = {}
            debugging_dict[s+':counts'] = {}
            #debugging_dict[s+':pvals'] = {}
        #loop by sample_names
        for sample_name, col in phip_df.items():
            #print(sample_name)
            for s, indexs in taxon_group.items():
                #1: inter-taxon searching
                inter_list = inter_df.ix[indexs][sample_name]
                inter_dict = self.taxon_inter_searching(col[indexs], inter_list)
                #export
                debugging_dict[sample_name+':all_hits'][s] = inter_dict['all_hits']
                debugging_dict[sample_name+':inter_hits'][s] = inter_dict['inter_hits']
                #print(inter_dict)
                
                #2: intra-taxon searching
                intra_dict = self.taxon_intra_searching(col[inter_dict['other_hits']])
                #export
                debugging_dict[sample_name+':intra_hits'][s] = intra_dict['intra_hits']
                debugging_dict[sample_name+':hits'][s] = intra_dict['hits']
                all_hits = ['{}:{}'.format('all', len(inter_dict['all_hits'])),
                        '{}:{}'.format('inter', len(inter_dict['inter_hits'])),
                        '{}:{}'.format('intra', len(intra_dict['intra_hits'])),
                        '{}:{}'.format('hits', len(intra_dict['hits'])) ]
                debugging_dict[sample_name+':counts'][s] = ','.join(all_hits)
                hit_list = ['({},{})'.format(a, b) for a, b in col[intra_dict['hits']].items()]
                taxon_pep_dict[sample_name][s] = ','.join(hit_list)
                #counts matrix of taxonomy search
                taxon_dict[sample_name][s] = len(intra_dict['hits'])
                
                #3: permutation
                #hit_scores = col[intra_dict['hits']]
                #permuted_scores = permute_dict[s]#df, pepids in rows, permuted scores in columns
                #pval_dict = self.taxon_permutation(hit_scores, permuted_scores, col)
                #export
                #pval_list = [len(intra_dict['hits']), pval_dict['ttest_pval'], pval_dict['utest_pval']]
                #taxon_pval_dict[sample_name][s] = ','.join(map(str, pval_list))
                #pval_list = [ a+':'+str(b) for a,b in pval_dict.items()]
                #debugging_dict[sample_name+':pvals'][s] = ','.join(pval_list)
        #export to file
        file_head = '{}_{}_'.format(myIO.file_os(count_file).file_prefix(), taxon_rank)
        taxon_dict = myDict.basic(taxon_dict).transform_dict2()
        myDict.basic(taxon_dict).dict2_to_file(file_head+'counting.txt', "\t")
        taxon_pep_dict = myDict.basic(taxon_pep_dict).transform_dict2()
        myDict.basic(taxon_pep_dict).dict2_to_file(file_head+'peptides.txt', "\t", 'NA')
        debugging_dict = myDict.basic(debugging_dict).transform_dict2()
        myDict.basic(debugging_dict).dict2_to_file(file_head+'debugging.txt', "\t", 'NA')
        #myDict.basic(taxon_pval_dict).dict2_to_file(file_head+'pvalues.txt', "\t", 'NA')
    

#search inter-taxonomy speicific    
    def taxon_inter_searching(self, score_list, inter_list):
        inter_dict = {'all_hits':[], 'inter_hits':[], 'other_hits':[]}
        for name, score in score_list.items():
            #condition 1: >threshold
            if score >= self.par['zscore_threshold']:
                inter_dict['all_hits'].append(name)
                #condition 2: inter-overlapping
                #with a given pep_id, check if its inter-taxon peptides are signficant based on hashlf z-scores
                if inter_list[name] >= (self.par['zscore_threshold']/2):
                    inter_dict['inter_hits'].append(name)
                else:
                    inter_dict['other_hits'].append(name)
        #print [(key,len(ids)) for key,ids in inter_dict.items()]
        return inter_dict
        
#search intra-taxonomy specific
    def taxon_intra_searching(self, score_list):
        #print score_list
        intra_dict = {'intra_hits':[], 'hits':[]}
        intra_dict['input_hits'] = list(score_list.index)
        if not score_list.empty:
            #1: rank by scores from hightest to lowest
            score_list.sort_values(ascending=False, inplace=True)
            intra_dict['hits'].append(score_list.index[0]) #the top one is hit
            #2: reserve independent peptides of intra-species
            for i in range(1, len(score_list)):
                query_pep = score_list.index[i]
                #overlapped pepids against query_pep
                if query_pep in self.par['dependent_pep']:
                    dependent_peps = self.par['dependent_pep'][query_pep]
                    #compare with higher ranked peptides
                    high_peps = list(score_list.index[:i])
                    intersect = list(set(high_peps) & set(dependent_peps))
                    if len(intersect) > 0:
                        intra_dict['intra_hits'].append(query_pep)
                    else:
                        intra_dict['hits'].append(query_pep)
                else:
                    intra_dict['hits'].append(query_pep)
        #                    
        #print [(key,len(ids)) for key,ids in intra_dict.items()]
        return intra_dict

#hit_scores are pd.Series, values are numberic type
#permuted_scores is df
    def taxon_permutation(self, hit_scores, permuted_scores, all_scores):
        #print(hit_scores)
        hits_num = len(hit_scores)
        pval_dict = {'utest_pval':1, 'ttest_pval':1, 'null_mean':0, 'null_sd':1}
        #get permuted scores, type is df
        #hit_names are in columns, permutation times are in rows
        #permuted_hit_df = permuted_scores.ix[hit_scores.index].copy()
        
        #signficance analysis
        if hits_num > 2:
            utest_pvals = []
            permuted_hits_nums = []
            #loop
            for t, shuffled_index in permuted_scores.items():
                shuffled_scores = all_scores[shuffled_index].copy() #pd.Series
                shuffled_scores.index = permuted_scores.index
                scores = shuffled_scores[hit_scores.index]
                #u test
                try:
                    perm_p = stats.mannwhitneyu(hit_scores, scores)[1]
                except ValueError:
                    perm_p = 1
                utest_pvals.append(perm_p)
                #t-test
                permuted_hits = scores[scores >= self.par['zscore_threshold']]
                permuted_hits_nums.append(len(permuted_hits))
            
            #calculate mean of pvalues of u-test
            pval_dict['utest_pval'] = np.mean(utest_pvals)
            #calculate null model of t-test
            pval_dict['null_mean'] = np.mean(permuted_hits_nums)
            pval_dict['null_sd'] = np.std(permuted_hits_nums)
            #pval
            try:
                p_ttest = stats.ttest_1samp(permuted_hits_nums, hits_num)[1]
                #one-sideed tail
                pval_dict['ttest_pval'] = p_ttest/2 if hits_num>pval_dict['null_mean'] else 1-(p_ttest/2)
            except ValueError:
                pass
        #    
        return pval_dict

#polyclonal
    def sig_polyclonal(self, count_file):
        #count_file = args_tuple
        print("Polyclonal analysis of ", count_file)
        comb_df, pep_df = myCommon.basic(self.par).combine_df(count_file)
        #functions
        def hits_func(x, peps, threshold, pro_id):
            #signficant hits
            hits = x[x>= threshold]
            #non_overlapping peptides
            peps = [str(x) for x in peps]
            hit_peps = [str(x) for x in hits.index]
            none_overlapped_hits_num = myList.basic(peps).un_neighbours(hit_peps, return_type='hits_num')
            #if none_overlapped_hits_num>1: print "%d,%d" %(len(list(hits.index)), none_overlapped_hits_num)
            #if len(hit_peps)>0: print pro_id, peps, hit_peps
            #if pro_id == 'Q9YLJ1': print pro_id, peps, hit_peps
            return len(list(hits.index)), none_overlapped_hits_num, ','.join(hit_peps)
            
        #collapse by protein
        hits1 = {}
        hits2 = {}
        #n = 1
        for pro_id, row_index in comb_df.groupby(self.par['protein_assoc']).groups.items():
            #row is protein id
            ##get protein-peptides annotations
            peps_str = self.par['dict_pro_pep'][pro_id]
            peps = peps_str.split(',')
            #df by protein
            sub_df = pep_df.ix[row_index]
            #print pro_id, list(sub_df.index)
            #hits num beyond zscore threshold
            hits_num = sub_df.apply(hits_func, axis=0, args=(peps, self.par['zscore_threshold'], pro_id))
            #if pro_id == 'Q9YLJ1': print hits_num
            #all number of significant hits
            num1 = [ h[0] for h in hits_num]
            hits1[pro_id] = dict(zip(list(sub_df), list(num1)))
            #number of sig hits without overlapping
            num2 = [ h[1] for h in hits_num]
            hits2[pro_id] = dict(zip(list(sub_df), list(num2)))
            #if (np.sum(num1))>10:
                #pd.set_option('display.max_columns', None)
                #pd.set_option('display.max_rows', None)
                #print np.matrix(np.round(sub_df))
                #print num1
                #print num2
            #n+ = 1
            #if n == 10: break
            
        #export
        file_head = myIO.file_os(count_file).file_prefix()
        myDict.basic(hits1).dict2_to_file(file_head+'_polyclonal.txt', "\t")
        myDict.basic(hits2).dict2_to_file(file_head+'_polyclonal_nonoverlapped.txt', "\t")
        
#specie alignment to remove overlapped peptides
#the output matrix: specie in rows and sample in columns
    def taxon_blast(self, file_aln, zscore_file):
        print('###Signficant taxon by removing overlapped hits based on blast alignment.')
        taxon_type = myIO.file_os(file_aln).name_prefix()
        print('{}: {}'.format(taxon_type, zscore_file))
        #read zscore_df
        zdf = myDataframe.basic().standard_df(zscore_file)

        #match order of align score and zscore,replace na
        #read alignment file for specie alignment
        binary_b = myDataframe.basic().aln_df(file_aln,self.par['align_score'])
        binary_b = binary_b.reindex(zdf.index).fillna(0)
        #print binary_b
        
        #sample names in columns, and specie in rows
        sum_df = pd.DataFrame(0, index=list(binary_b), columns=list(zdf))
        pep_df = pd.DataFrame(np.nan, index=list(binary_b), columns=list(zdf))
        #perm_df = pep_df.copy()
        #print binary_z.apply(sum, axis = 0)
        #n = 1
        for sample_name, column in zdf.items():
            #1: select peptides
            #column = zscore_df.ix[:,20]
            #first remove all nont-hits 
            hits = column[column >= self.par['specieZ_threshold']].copy()#all hits
            hits.sort_values(axis = 0, ascending = False, inplace=True)
            #print hits
            #remove overlapped hits
            nonoverlap_hits,overlap_debug = myList.basic(hits).remove_overlap(self.par['dependent_pep'])
            input_num = len(nonoverlap_hits)
            print('{}: hits={}, nonoverlapped={}'.format(sample_name, len(hits), input_num))

            #2: remove overlap hits between species
            if input_num > 0:
                ###2-1: export peptides
                try:
                    outfile = '{}{}/{}.csv'.format(self.par['dir_result'], sample_name, taxon_type)
                    overlap_debug.to_csv(outfile, header=True, index_label='peptides')
                except FileNotFoundError:
                    pass
                ###2-2: specie-specific hits based on non-overlapped hits
                #sample zscore-alignscore matrix times by zscore
                #print(nonoverlap_hits.index)
                zb_df = binary_b.ix[nonoverlap_hits.index]
                #print(list(binary_b.apply(lambda x: sum(x), axis = 0)))
                #loop
                collapse_zb, sim_tag = myDataframe.basic(zb_df).unispecie(self.par['sim_threshold'])
                #counts of hits
                sum_df[sample_name] = collapse_zb.apply(sum,axis = 0) + sim_tag
                #print(list(sum_df[sample_name]))
                #high_sum = sum_df[sample_name]
                #print(high_sum[high_sum>0])
                #all peptide_id list
                pep_df[sample_name] = collapse_zb.apply(lambda x: myList.basic(x).names_string(0.001),axis=0)
                
                #2-3:permutation
                #perm_df[sample_name] = self.specie_alignment_permutation(input_num)
            #if n == 10: break
            #n+ = 1
        #export to file
        file_head = '{}_{}_'.format(myIO.file_os(zscore_file).file_prefix(), taxon_type)
        sum_df.to_csv(file_head+'counting.txt', sep='\t', header=True, index_label='Specie')
        pep_df.to_csv(file_head+'peptides.txt', sep='\t', header=True, index_label='Specie')
        #perm_df.to_csv(file_head+'permutation.txt', sep = '\t', header = True, index_label = 'Specie')
        #
        #return zscore_df, binary_b

#coding by Sanjay
#specie alignment to remove overlapped peptides
#the output matrix: specie in rows and sample in columns
#input functions: aln_df(), filter_aln(), and binom_unispecie in myDataframe
#input functions: gen_ind_hits() in myList
    def taxon_blast2(self, file_aln, zscore_file):
        taxon_type = myIO.file_os(file_aln).name_prefix()
        print("\n{}:{}\n".format(taxon_type, zscore_file))
        #read zscore_df
        zdf = myDataframe.basic().standard_df(zscore_file)

        #match order of align score and zscore,replace na
        #read alignment file for specie alignment
        binary_b = myDataframe.basic().aln_df(file_aln, self.par['align_score'])
        #binary_b = myDataframe.basic(binary_b).filter_aln()
        binary_b = binary_b.reindex(zdf.index).fillna(0)
        
        #print binary_b
        
        #sample names in columns, and specie in rows
        sum_df = pd.DataFrame(0, index=list(binary_b), columns=list(zdf))
        pep_df = pd.DataFrame(np.nan, index=list(binary_b), columns=list(zdf))
        p_df = pd.DataFrame(index=list(binary_b), columns=list(zdf))
        #perm_df=pep_df.copy()
        #print binary_z.apply(sum, axis=0)
        #n=0
        for sample_name, column in zdf.iteritems():
            #n += 1
            #1: select peptides
            #column=zscore_df.ix[:,20]
            #first remove all nont-hits 
            hits = column[column>=self.par['specieZ_threshold']].copy()#all hits
            hits.sort_values(axis=0, ascending=False, inplace=True)
            #print hits
            #remove overlapped hits
            nonoverlap_hits = myList.basic(hits).gen_ind_hits(self.par['dependent_pep'])
            input_num = len(nonoverlap_hits)
            print("{}:\thits={}, nonoverlapped={}".format(sample_name, len(hits),input_num))

            #2: remove overlap hits between species
            if input_num>0:
                zb_df = binary_b.loc[nonoverlap_hits.index]
                #print list(binary_b.apply(lambda x: sum(x), axis=0))
                #loop
                collapse_zb, sim_tag, p_series = myDataframe.basic(zb_df).binom_unispecie(self.par['dir_ref_seq'], input_num, self.par['p_threshold'], self.par['x_threshold'])
                #counts of hits
                sum_df[sample_name] = collapse_zb.apply(sum, axis=0)+sim_tag
                #all peptide_id list
                pep_df[sample_name] = collapse_zb.apply(lambda x: myList.basic(x).names_string(0.001), axis=0)
                p_df[sample_name] = p_series
                #padjust_df[sample_name]=p_adjust_series
            #if n==5: break
            #n+=1
        #export to file
        file_head = myIO.file_os(zscore_file).file_prefix()+'_'+taxon_type+'_'
        #file_head='random_min_HI_HC_'+taxon_type+'_'
        sum_df.to_csv(file_head+'counting.txt', sep='\t', header=True, index_label='Specie')
        pep_df.to_csv(file_head+'peptides.txt', sep='\t', header=True, index_label='Specie')
        p_df.to_csv(file_head+'p-values.txt', sep='\t', header=True, index_label='Specie')
        
        #Adjusted p-values using B-H
        '''
        stats = importr('stats')
        for i in p_df:
            pvalue_list = p_df[i].values
            p_adjust = list(stats.p_adjust(FloatVector(pvalue_list), method = 'BH'))
            padjust_df[i] = p_adjust
        padjust_df.to_csv(file_head+'p-adjusted.txt', sep='\t', header=True, index_label='Specie')
        '''
        padjust_df=pd.DataFrame(index=list(binary_b), columns=list(zdf))
        for i in p_df.columns:
            pvals=np.array(p_df[i].values)
            if not np.isnan(pvals).all():
                mask = [j for j in np.where(np.isfinite(pvals))[0]]
                pval_corrected = np.empty(pvals.shape)
                pval_corrected.fill(np.nan)
                pval_corrected[mask] = multipletests(pvals[mask], method='fdr_bh')[1]
                padjust_df[i] = pval_corrected
        padjust_df.to_csv(file_head+'padjusted.txt', sep='\t', header=True, index_label='Specie')
        #perm_df.to_csv(file_head+'permutation.txt', sep='\t', header=True, index_label='Specie')


 
#end