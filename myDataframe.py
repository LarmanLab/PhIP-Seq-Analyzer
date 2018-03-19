# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 14:45:17 2016

@author: yuan
"""
#standard modules
import re
import pandas as pd
import numpy as np
from scipy.stats import binom_test

#in-house modules
import myDict
import myIO
import myList
import myPlot

#
class basic:
    def __init__(self, df=None):
        self.df = df

#shrink df to Series
#axis=0, shrink values into a joined value str in each column 
    def to_series(self, sep=',', axis=0):
        self.df = pd.DataFrame(self.df)
        s = pd.Series()
        if axis == 0:
            for name, col in self.df.iteritems():
                s[name] = sep.join(str(x) for x in col)
        elif axis == 1:
            for name, row in self.df.iterrows():
                #print name
                s[name] = sep.join(str(x) for x in row)
        #print s
        return s
        
#read standard data frame: 
#rowname in the first  column, columnname in the first row, both names are the string type
#replace all NAN with zero
    def standard_df(self, infile, fill=True):
        #read txt file
        sep = ',' if myIO.file_os(infile).name_suffix() == 'csv' else '\t'
        stand_df = pd.read_table(infile,header=0, index_col=0,sep=sep,low_memory=False)
        #string of row names
        stand_df.index = [str(x) for x in list(stand_df.index)]
        stand_df.columns = [str(x) for x in list(stand_df)]
        #replace NAN
        if fill == True:
            stand_df.fillna(0, inplace=True)
        return stand_df

#read annotation file
    def annot_df(self, infile):
        file_sep = '\t' if myIO.file_os(infile).name_suffix() == 'txt' else ','
        annot_df = pd.read_csv(infile, header=0, index_col=None, sep=file_sep, low_memory=False)  
        #bother column and row names should be string type
        annot_df.index = annot_df['pep_id']
        annot_df.index = annot_df.index.astype(str)
        #annot_df.columns=self.par['annot_df'].columns.astype(str)
        return annot_df

#read alignment file for specific specie alingment searching
    def aln_df(self, infile, align_score):
        print("\tRead align file ", infile)
        file_sep = '\t' if myIO.file_os(infile).name_suffix() == 'txt' else ','
        aln_df = myIO.file_os(infile, file_sep).flat_file_to_df([0,1,11], True)
        #print aln_df
        #convert to binary matrix
        binary_b = pd.DataFrame(np.where(aln_df >= align_score, 1, 0))
        #print np.sum(binary_b)
        binary_b.index = [str(x).split('_')[1] for x in list(aln_df.index)]
        binary_b.columns = [re.sub(',', ';', str(x)) for x in list(aln_df)]
        binary_b.fillna(0, inplace=True)
        #print binary_b.ix['Rnl2_SPIKEIN']
        #print list(binary_b.index)[:10]
        '''
        ##remove subset species, of which alignment score <1
        #calculate shared probes
        taxon_sums = binary_b.apply(sum, axis=0)
        taxons = list(binary_b)
        shared_prob=pd.DataFrame(index=taxons, columns=taxons)
        for taxon, col in binary_b.iteritems():
            shared_prob[taxon]=binary_b.apply(lambda x, y=col: np.dot(y,x), axis=0)/taxon_sums
            shared_prob[taxon][taxon]=0
        #
        tags_df=pd.DataFrame(index=taxons, columns=taxons)
        for taxon, col in shared_prob.iteritems():
            row=shared_prob.loc[taxon]
            tags_df[taxon]=[1 if a!=1 and b==1 else 0 for a, b in zip(col,row)]
        sum_tags=tags_df.apply(max, axis=0)
        reserved_taxons=list(sum_tags[sum_tags==0].index)
        binary_b=binary_b[reserved_taxons]
        '''
        return binary_b

#remove subset species, of which alignment score <1
#input is binary_b 
    def filter_aln(self):
        print('remove subset species.')
        #aln_df = pd.read_csv("aln_nonsubset_matrix.csv", header=0, index_col=0)
        #binary_b=pd.DataFrame(np.where(aln_df > 0, 1, 0), columns=list(aln_df.keys()), index=list(aln_df.index))
        #del aln_df
        virus_sums = self.df.apply(sum, axis=0)
        #virus_sums.to_csv("noseg_virus_total_aligned.csv")
    
        viruses = list(self.df.columns)
        virus_intersections = pd.DataFrame(index=viruses, columns=viruses)
        for i in viruses:
            a = self.df[i]
            for j in viruses:
                b = self.df[j]
                virus_intersections.loc[i,j] = np.dot(a,b)
        #virus_intersections.to_csv("noseg_virus_intersections.csv", header=True, index=True)
        #print(virus_intersections)
        shared_prob = pd.DataFrame()
        for i in virus_intersections:
            shared_prob[i] = virus_intersections[i]/virus_sums[i]
            #shared_prob.to_csv("noseg_virus_shared_probabilities.csv", header=True, index=True)
        binary_copy = self.df.copy()
        for i in list(shared_prob.columns):
            for j in list(shared_prob.index):
                if i!=j and shared_prob.loc[j,i] == 1 and shared_prob.loc[i,j] != 1:
                    #print(i+","+j)
                    if i in binary_copy:
                        binary_copy.drop(i, axis=1, inplace=True)
        #binary_copy.to_csv("noseg_binary_aln_reduced.csv", index=True, header=True)
    
        virus_sums = binary_copy.apply(sum, axis=0)
        #virus_sums.to_csv("sums_test.csv")
        virus_sums_index = virus_sums.index
        for i in range(len(virus_sums)):
            if virus_sums.iloc[i] <= 10:
                binary_copy.drop(virus_sums_index[i], axis=1, inplace=True)
    
        #binary_copy.to_csv("noseg_binary_aln_reduced.csv", index=True, header=True)    
        return binary_copy
    
    #convert a string to list
    def string_to_list(self, string, sep1=None, sep2=None):
        outlist = []
        if sep1 is None:
            outlist = [string]
        else:
            if sep2 is None:
                outlist = string.split(sep1)
            else:
                tmplist = string.split(sep1)
                outlist = [ t.split(sep2)[0] for t in tmplist]
        return outlist
        
#extract unique items from a given column of a dataframe
    def df_list(self, col_name, sep1=None, sep2=None):
        outlist = []
        out_dict = {}
        try:
            inlist = list(self.df[col_name])
        except:
            print('No column name in the data frame:', col_name)
        else:
            if sep1 is None:
                outlist = list(set(inlist))
            else:
                #out_dict={}
                if sep2 is None: #sep1 only
                    for ele1 in inlist:
                        ele1 = str(ele1)# avoid null list export
                        ele1_list = ele1.split(sep1)
                        for item in ele1_list:
                            out_dict[item] = 0
                else:#sep1 and spe2 both
                    out_dict = {}
                    for ele1 in inlist:
                        ele1 = str(ele1) # avoid null list export
                        ele1_list = ele1.split(sep1)
                        for ele2 in ele1_list:
                            item = ele2.split(sep2)[0]
                            out_dict[item] = 0
                #unique  items
                outlist = myList.basic(out_dict.keys()).sort_list()
        finally:
            #print outlist
            pass
            return outlist
#old version
    def df_list0(self, col_name, sep1=None, sep2=None):
        outlist = []
        inlist = list(self.df[col_name])
        if sep1 is None:
            outlist = list(set(inlist))
        else:
            out_dict = {}
            if sep2 is None: #sep1 only
                for ele1 in inlist:
                    ele1_list = ele1.split(sep1)
                    for item in ele1_list:
                        out_dict[item] = 0
            else:#sep1 and spe2 both
                print('sep12')
                for ele1 in inlist:
                    ele1 = str(ele1)
                    ele1_list = ele1.split(sep1)
                    for ele2 in ele1_list:
                        item = ele2.split(sep2)[0]
                        out_dict[item] = 0
                    #print len(out_dict.keys())
            #unique  items
            #print(out_dict)
            outlist = myList.basic(out_dict.keys()).sort_list()
        return outlist
        
#multiple vs multiple by two given columns
#generate a dict: key(from colA) vs list(from colB)
    def list_dict(self, colA, colB, sep1=None, sep2=None):
        out_dict = {}
        if sep1 is None:
            for A, B in zip(self.df[colA], self.df[colB]):
                out_dict[A] = [B]
        else:
            listB = [str(x) for x in self.df[colB]]
            if sep2 is None: #sep1 only
                for A, B in zip(self.df[colA], listB):
                    out_dict[A] = B.split(sep1)
            else:#sep1 and sep1 both
                for A, B in zip(self.df[colA], listB):
                    tmp_list = B.split(sep1)
                    out_dict[A] = [t.split(sep2)[0] for t in tmp_list]
        #
        return out_dict

#export data frame to a file and generate plot of sample vs number of hits 
#Always samples in columns
    def export_df(self, outfile, threshold=10, index_label='row_names'):
        print('\texport data frame to ', outfile)
        outsep = ',' if outfile.endswith('.csv') else '\t'
        self.df.to_csv(outfile, sep=outsep, index_label=index_label) 
        
        #draw a scatterplot
        counts = self.df.apply(lambda x, y=threshold: len(x[x>=y]), axis=0)
        #print counts
        plot_par={'list':counts, 'ylabel':'Sample_names', 'xlabel':'Number of hits', 
                  'picfile': myIO.file_os(outfile).file_prefix()+'.png',
                  'title': 'Number of hits, threshold='+str(threshold) }
        myPlot.plot(plot_par).simple_barh()

#embed data frames, output is list-nested df
#input is list of df, self.df is none
    def embed_dflist(self, df_list):
        #combine data frame by rows
        concat_df = pd.concat(df_list, axis=0)
        #print concat_df.shape
        colnames = list(concat_df)
        #reset index
        concat_df = concat_df.reset_index()
        #print concat_df.shape
        #print list(concat_df)
        #group by index column
        rows_obj = concat_df.groupby(list(concat_df)[0])            
        
        def flat_list(x):
            l = []
            for i in x:
                #print i
                l = np.append(l,i)
            #print l
            return list(l)
        
        #embed values into an array        
        embed_dict = {}
        for row_name, indexs in rows_obj.groups.iteritems():
            sub_df = concat_df.ix[indexs][colnames].copy()
            #print sub_df
            values = np.array(sub_df.apply(flat_list, axis=0).transpose())
            #values=np.array(sub_df.transpose())
            #print values
            embed_dict[row_name] = dict(zip(colnames, values))
        #convert dict to df
        embed_df = pd.DataFrame(embed_dict).transpose()
        #print embed_df
        return embed_df
            
#input is only two df: self.df and df2
    def embed_dfs(self, df2):
        if self.df.empty:
            embed_df = df2
        else:
            embed_df = self.embed_dflist([self.df, df2])
        return embed_df
            
#permute data frame by columns, and return a new one
#it might be slow and occupy too much memory if df is huge
    def permute_col(self, times=2, slice_dict=None):
        #add shuffled dict into embed_dict
        embed_dict = {}
        for i in range(times):
            #1: shuffle data frame
            shuffled_df = self.df.iloc[np.random.permutation(len(self.df))].copy()
            shuffled_df.index = self.df.index
             #print shuffled_df
            #2: convert permuted dfs into embeded dataframe
            shuffled_dict=shuffled_df.to_dict()
            embed_dict=myDict.basic(embed_dict).combine_dupdict2(shuffled_dict, i)
            
        #convert to dataframe
        if slice_dict is None:
            permute = embed_dict # col-name is key1, row-name is key2
        else:
            permute = {}
            embed_df = pd.DataFrame(embed_dict)
            for slice_name, row_indexs in slice_dict.iteritems():
                permute[slice_name] = {}
                #print slice_name
                #
                sub_df = embed_df.ix[row_indexs]
                for col_name, col in sub_df.iteritems():
                    permute[slice_name][col_name] = col #col is pd.Series
                    #print col
                    #print type(col)
                    #break
                
        #print permute
        return(permute)        


#interacte matrix
#A_B is kind of serires with A as names, B as values
    def interact_df(self, A_B, func, outfile='NA'):
        print('Get the interacted data frame ', outfile)
        #initiate
        inter_dict = {}
        inter_df = pd.DataFrame(index=self.df.index, columns=list(self.df) )
        #A_names = list(self.df.index)
        #remove the value with 'none'
        A_B = A_B[A_B!='none']
        for A_name, B_str in A_B.iteritems():
            B_names = str(B_str).split(',')
            #print B_names
            Bdf = self.df.ix[B_names]
            inter_dict[A_name] = list(Bdf.apply(func, axis=0))
            #print inter_df.ix[A_name]
        #get the whole data frame
        sdf = pd.DataFrame(inter_dict).transpose()
        sdf.columns = list(self.df)
        #replace values
        inter_df.ix[sdf.index] = sdf
        #print inter_df
        if outfile != 'NA':
            inter_df.to_csv(outfile, sep="\t", na_rep='nan')
        #print inter_df
        return inter_df           

#for specie alignment of virus
#the input is binary matrix from zscore-alignscore
    def unispecie(self, sim_threshold=0.8):
        #r=self.df.apply(sum, axis=0)
        #r.sort(axis=0,ascending=False)
        #print(r)
        #initiate ranked hits
        ranked_hits = self.df.apply(sum, axis=0)
        ranked_hits.sort_values(axis=0, ascending=False, inplace=True)
        #print ranked_hits[:10]

        #set conditions 
        tag = 0
        #initiate similar series
        sim_tag = pd.Series(float(0), index=list(self.df))
        while len(ranked_hits) > 0 and ranked_hits[0] > 0:
            #calculate similarity
            specie = list(ranked_hits.index)
            #print(len(specie))
            #get hits of the highest specie
            highrank = self.df[specie[0]]
            hits_name = tuple(highrank[highrank>=1].index)#peptide names
            
            #1: similarity
            sim = self.df[specie].apply(lambda x, y=highrank: sum(y*x)/sum(y), axis=0)
            #print 'Similarity:' 
            #print list(sim[sim>0])
            #get all dis-similiar specie against the highest specie
            sim_specie = list(sim[sim>=sim_threshold].index)
            dissim_specie = tuple(sim[sim<sim_threshold].index)
            
            #2: reset
            if len(dissim_specie) > 0:
                #print list(self.df.ix[hits_name,dissim_specie].apply(sum, axis=0))
                self.df.set_value(hits_name,dissim_specie, 0)
                #print list(self.df.ix[hits_name,dissim_specie].apply(sum, axis=0))
            if len(sim_specie) > 1:
                tag += 1
                #print list(sim_tag[sim_specie])
                sim_tag.ix[sim_specie] = float('.%03d' % tag) #tag indicate them as sim
                #if tag==1: print ranked_hits[sim_specie]
            
            #refresh the list specie
            specie = list( set(specie) - set(sim_specie) )#remove those similar species
            #refresh conditions
            ranked_hits = self.df[specie].apply(sum, axis=0)
            ranked_hits.sort_values(axis=0, ascending=False, inplace=True)
        #
        return self.df, sim_tag
########################



#for reassign hits and similar tags
    def binom_unispecie(self, prob_dir, input_num, p_threshold=0.001, hits_threshold=1):
        #Read probability files
        first_round_prob = pd.read_csv(prob_dir+"avg_total_probabilities.csv", index_col=0, header=None, squeeze=True)
        second_round_prob = pd.read_csv(prob_dir+"avg_unique_probabilities.csv", header=0, index_col=0)
        third_round_prob = pd.read_csv(prob_dir+"avg_shared_probabilities.csv", header=0, index_col=0)
        
        #Function for defining significance
        def is_sig(p, n, x):
            p_value = binom_test(p=p, n=n, x=x, alternative='greater')
            if p_value < p_threshold and x > hits_threshold:
                return True
            else:
                return False
        
        #Series of p-values for each virus
        p_series = pd.Series(index=list(self.df.columns))
        
        #Number of hits to each virus
        ranked_hits = self.df.apply(sum, axis=0)
        #Calculate p-values for each virus based on the number of total hits
        virus_pvalues_1 = pd.Series(index=list(self.df.columns))
        for i in virus_pvalues_1.index:
            virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=input_num, x=ranked_hits[i], alternative='greater')
        
        #Sort viruses by initial p-value
        virus_pvalues_1.sort_values(inplace=True, ascending=True)
        #Sort virus hits series by p-value
        ranked_hits = ranked_hits[list(virus_pvalues_1.index)]        
        specie=list(ranked_hits.index)
        n_rank = input_num
        
        while len(ranked_hits)>0 and ranked_hits.iloc[0]>0:
            #Comparing top hit to everything else (reassignments only)
            for i in specie[1:]:
                #Top hit virus (binary vector)
                highrank = self.df[specie[0]]
                high_num = sum(highrank)
                i_rank = self.df[i]
                i_num = sum(i_rank)
                #element-wise multiplication to get the vector of shared peptides
                shared_peps = np.multiply(i_rank, highrank)
                shared_num = sum(shared_peps)
                #Only do reassignment/sim tags if there is overlap between the two viruses
                if shared_num > 0:
                    #If only one passes the threshold, reassign peptides accordingly
                    if is_sig(second_round_prob.loc[i,specie[0]], n_rank-shared_num, high_num-shared_num) and not is_sig(second_round_prob.loc[specie[0],i], n_rank-shared_num, i_num-shared_num):
                        #Subtract shared peptides from the insignificant virus
                        self.df[i] = i_rank-shared_peps
                    elif not is_sig(second_round_prob.loc[i,specie[0]], n_rank-shared_num, high_num-shared_num) and is_sig(second_round_prob.loc[specie[0],i], n_rank-shared_num, i_num-shared_num):
                        #Same as above
                        self.df[specie[0]] = highrank-shared_peps
            
            #Add p-value to series after any potential reassignments
            top_hit = specie[0]
            p_series[top_hit] = binom_test(p=first_round_prob[top_hit], n=n_rank, x=sum(self.df[top_hit]), alternative='greater')
            
            #Figure out how many peptides were globally unique to highrank
            high_peptides = np.where(self.df[specie[0]] == 1)[0]
            high_unique = len(np.where(self.df.iloc[high_peptides,:].apply(sum, axis=1) == 1)[0])
            
            #Now remove the highest hit since it will not be involved in subsequent comparisons
            specie.remove(specie[0])
            #Re-rank virus hits by total binomial
            ranked_hits = self.df.apply(sum, axis=0)
            ranked_hits = ranked_hits[specie]
            virus_pvalues_1 = pd.Series(index=specie)
            #Adjust the n used for binomial tests
            n_rank -= high_unique
            
            #Sort viruses by p-value
            for i in specie:
                virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=n_rank, x=ranked_hits[i], alternative='greater')
            virus_pvalues_1.sort_values(inplace=True, ascending=True)
            specie = list(virus_pvalues_1.index)
            ranked_hits = ranked_hits[specie]
        
        #DOING SIM TAGS AFTER REASSIGNMENT '''
        ranked_hits = self.df.apply(sum, axis=0)
        virus_pvalues_1 = pd.Series(index=list(self.df.columns))
        for i in virus_pvalues_1.index:
            virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=input_num, x=ranked_hits[i], alternative='greater')
        virus_pvalues_1.sort_values(inplace=True, ascending=True)
        ranked_hits = ranked_hits[list(virus_pvalues_1.index)]
        
        #Sim tags for each virus, list of species to be examined
        specie=list(ranked_hits.index)
        sim_tag = pd.Series(float(0), index=specie)
        tag = 0
        
        n_rank = input_num
        
        while len(ranked_hits)>0 and ranked_hits.iloc[0]>0:
            for i in specie[1:]:
                highrank = self.df[specie[0]]
                high_num = sum(highrank)
                i_rank = self.df[i]
                i_num = sum(i_rank)
                shared_peps = np.multiply(i_rank, highrank)
                shared_num = sum(shared_peps)
                #Only do sim tags if there is overlap between the two viruses
                if shared_num > 0:
                    #If neither passes (using shared test, symmetric probability table) and neither can stand on their own:
                    if is_sig(third_round_prob.loc[i,specie[0]], n_rank-(high_num-shared_num)-(i_num-shared_num), shared_num) and not is_sig(second_round_prob.loc[i,specie[0]], n_rank-shared_num, high_num-shared_num) and not is_sig(second_round_prob.loc[specie[0],i], n_rank-shared_num, i_num-shared_num):
                        #If neither of them already has a sim tag
                        if np.array_equal(sim_tag[[specie[0],i]], [0,0]):
                            tag += 1
                            #Add a tag to the list for both viruses
                            sim_tag[[specie[0],i]] = 0.001*tag
                        #If either already has a tag, assign that tag to the other
                        elif sim_tag[specie[0]] != 0 and sim_tag[i] == 0:
                            sim_tag[i] = sim_tag[specie[0]]
                        elif sim_tag[i] != 0 and sim_tag[specie[0]] == 0:
                            sim_tag[specie[0]] = sim_tag[i]
            
            #Figure out how many peptides were globally unique to highrank
            high_peptides = np.where(self.df[specie[0]] == 1)[0]
            high_unique = len(np.where(self.df.iloc[high_peptides,:].apply(sum, axis=1) == 1)[0])
            #Now remove the highest hit since it will not be involved in subsequent comparisons
            specie.remove(specie[0])
            #Re-rank virus hits by total binomial
            ranked_hits = self.df.apply(sum, axis=0)
            ranked_hits = ranked_hits[specie]
            virus_pvalues_1 = pd.Series(index=specie)
            n_rank -= high_unique
            for i in specie:
                virus_pvalues_1[i] = binom_test(p=first_round_prob[i], n=n_rank, x=ranked_hits[i], alternative='greater')
            #Sort viruses by p-value
            virus_pvalues_1.sort_values(inplace=True, ascending=True)
            #Sort virus hits series by p-value
            specie = list(virus_pvalues_1.index)
            ranked_hits = ranked_hits[specie]
        
        #Generate adjusted p-values for the p-value output Series (using the R package)
        #stats = importr('stats')
        #p_adjusted = stats.p_adjust(FloatVector(p_series.values), method='BH')​
        
        print("\tBinomial testing for unique species is Done!")
        return self.df, sim_tag, p_series#, p_adjust_series

#
#end