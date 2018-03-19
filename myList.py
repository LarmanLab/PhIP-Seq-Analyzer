# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 10:31:52 2015

@author: yuan
"""
import collections
import random
#import networkx as nx
import timeit
import numpy as np
import pandas as pd

class basic:
    def __init__(self, array):
        self.list = array
        self.list_len=self.list.__len__()
        self.out_list=[]
        #print(self.list)
    
    #return tuple of min and max with a given list
    def min_max(self):
        #remove missing data
        numeric_list=[x for x in self.list if not x in [np.nan, -np.inf, np.inf] ]
        min_max=(min(numeric_list), max(numeric_list))
        #print(min_max)
        return min_max
    
#sort mixed list
    def sort_list(self):
        #mixed list is numeric and string mixed together
        str_list = []
        for element in self.list:
            if type(element) is str:
                str_list.append(element)
            else:
                self.out_list.append(element)
        #sort string list
        str_list.sort()
        #sort numeric list
        self.out_list.sort()
        #combine them
        self.out_list.extend(str_list)
        #print(self.out_list)
        return self.out_list

#remove duplicates, method 1
    def remove_duplicate1(self):
        duplicate = []
        while self.list.__len__() > 0:
            x = self.list.pop()
            if x in duplicate:
                pass
            else:
                if x in self.out_list:
                    self.out_list.remove(x)
                    duplicate.append(x)
                else:
                    self.out_list.append(x)
        #print(self.out_list) #unique
        #print(duplicate)
        return self.out_list

#remove duplicates, method2
    def remove_duplicate(self):
        x = 0  #index of unique boundary
        y = self.list_len # index of duplicate boundary
        while x < y:
            A = self.list[x]
            #print A, self.list
            dup = [x]
            for i in range(x+1, y):
                if self.list[i] == A:
                    dup.append(i)
            #
            if dup.__len__() == 1:
                x += 1
            else:
                #swithch index if there are duplicates
                for d in dup[::-1]:
                    y -= 1
                    self.list[d], self.list[y] = self.list[y], self.list[d]
            #print self.list, '\n'
        #print 'unique:', self.list[:x]
        #print 'duplicates:', self.list[y:]
        return self.list[:x]
        
#shrink a huge list
#select values from a sorted list with equal intervals
    def interval_list(self, points=1e3):
        if self.list_len > points:
            interval = int(self.list_len/points)
            flag = 0
            while flag < (self.list_len-1):
                self.out_list.append(self.list[flag])
                flag = flag + interval
            #always reserve the last element 
            self.out_list.append(self.list[-1])
        else:
            self.out_list = list(self.list)
        #
        return self.out_list #shrinked
    
#with a given list( sorted or not)
#return unjunct elements against another list pool
    def un_neighbours(self, in_list, return_type='hits_num'):
        #read intersected elements
        tmp_dict = dict([(x,0) for x in self.list])
        for e in in_list:
            tmp_dict[e] = 1
        #
        i = 0
        while i < self.list_len:
            a = self.list[i]
            if tmp_dict[a] == 1:            
                self.out_list.append(a)
                i += 2
            else:
                i += 1
        #
        #print unjunct_list
        if return_type == 'hits_num':
            return self.out_list.__len__()
        else:
            return self.out_list
            
#counts number of items of a compound list
#one on multiple items seperated by comma
    def elements_frequency(self, sep1=None, sep2=None):
        self.list = [str(x) for x in self.list]
        if sep1 is None:
            counts_dict = self.elements_frequency0()
        else:
            if sep2 is None:
                counts_dict = self.elements_frequency1(sep1)
            else:
                counts_dict = self.elements_frequency2(sep1, sep2)
        return counts_dict
#directly counts frequency such as [1,2,3,4,2,5]
    def elements_frequency0(self):
        items = [str(x) for x in self.list]
        #counts frequency of elements 
        counter_obj = collections.Counter(items)
        counts_dict = dict(zip(counter_obj.keys(),counter_obj.values()))
        #print counts_dict
        return counts_dict
    #Example: ['GO:0016021,GO:0019031,GO:0019050', 'GO:0003677,GO:0003887,GO:0006260']
    def elements_frequency1(self, sep):
        items = []
        for item in self.list:
            items1 = item.split(sep)
            items = items+items1
        #counts frequency of elements 
        counts_dict = array(items).elements_frequency0()
        return counts_dict
    #Example: ['PS00854,LPTDL;PS00441,LAL;PS01167,ARR, 'PS00778,ITAA;PS01324,GLCA;PS01221,PPGH;PS00063,TAA']
    def elements_frequency2(self, sep1, sep2):
        items = []
        for item in self.list:
            items1 = item.split(sep1)
            items2 = [r.split(sep2)[0] for r in items1]
            items = items+items2
        #counts frequency of elements 
        counts_dict = array(items).elements_frequency0()
        return counts_dict
                      
#export list to a text file seperated by return
    def list_to_file(self, out_file):
        out_obj = open(out_file, 'wt')
        for key in self.sort_list():
            out_obj.write(str(key)+'\n')
        out_obj.close()
        print('write a list to ', out_file)

#convert nested list to flat list
    def flat_list(self):
        for i in self.list:
            if isinstance(i,list): #judge the list type
                self.out_list += i
            else:
                self.out_list.append(i)
        return self.out_list  

#permute list, and return a permuted nested dict
    def permute_list(self, times=2, sampling_num=1):
        sampling_df = pd.DataFrame()
        for i in range(1,times+1):
            s = list(self.list)
            np.random.shuffle(s)
            sampling_df[i] = s[:sampling_num]
        #print(shuffled_df.ix[:10])
        #print(sampling_df)
        return sampling_df
            
#permute series, and return a permuted nested dict
    def permute_Series(self, times=2, slice_dict=None):
        self.list=pd.Series(self.list)
        #print(self.list)
        #print len(self.list)
        #initiate df
        shuffled_df = pd.DataFrame()
        for i in range(1, times+1):
            s = self.list.copy()
            np.random.shuffle(s)
            shuffled_df[i] = s
        #print shuffled_df.ix[:10]
        
        #output is dataframe or dict
        if slice_dict is None:
            permute = shuffled_df # not colnames, row-names is names of self.series
        else:
            permute = {}
            for slice_name, row_indexs in slice_dict.items():
                permute[slice_name] = shuffled_df.ix[row_indexs].copy()
                #print list(permute[slice_name].index)
                #print np.array(permute[slice_name])
        #print(permute)
        return(permute)            

#reserve the peptide with the highest value, and remove other overlapping peptides
#input: descending-ordered pd.Series
    def remove_overlap(self, overlap_dict):
        hits = pd.Series(self.list)
        debug = pd.DataFrame({'zscore':hits, 'overlap':hits.index})
        #remove overlapped hits
        indep_hits = pd.Series()
        while hits.__len__() > 0:
            #print len(hits)
            #unique hits
            high_name = hits.index[0]
            indep_hits[high_name] = hits[high_name]
            if high_name in overlap_dict:
                overlap_names = set(overlap_dict[high_name])
                other_names = set(hits.ix[1:].index)
                #intersection
                overlap_hits = list(overlap_names & other_names)
                #print len(overlap_hits)
                #remove overlapped lower hits
                hits = hits.drop([high_name]+overlap_hits)
                debug.set_value(tuple(overlap_hits), 'overlap', high_name)
            else:
                hits = hits.drop([high_name])
        #print indep_hits
        return indep_hits, debug
        
#input is pd.series. output is a name string of hit names
    def names_string(self, cutoff=0.001):
        hits = pd.Series(self.list)
        hits = hits[hits>=cutoff]
        names_str = ';'.join(list(hits.index) )
        return names_str

#convert duplicates in numeric list to unique values
#used for read counts only
#by adding a little decimal digits
    def jitter_duplicates(self):
        #count number of identical digits
        counter_obj = collections.Counter(self.list)
        freq = dict(zip(counter_obj.keys(),counter_obj.values()))
        #print(freq)
        
        #convert duplicate to unique values
        dict_unique = {}
        for key,value in freq.items():
            key = float(key)
            if value == 1:
                dict_unique[key] = [key]
            else:
                #get little decimal digits
                dicimal = [key+i/1e6 for i in range(value)]
                random.shuffle(dicimal)
                dict_unique[key] = dicimal
        #export list keep the same order
        for x in self.list:
            element = dict_unique[x].pop()
            self.out_list.append(element)
        #print(self.out_list)
        return self.out_list

#remove overlapped peptides based network topologic structure
    def gen_ind_hits(self, overlap_dict):
        #hits series
        hits_series = pd.Series(self.list)
        #Start time, for timing purposes
        start_time = timeit.default_timer()
        #Take only key-value pairs from overlap_dict for the sample hits
        hits_dict = hits_series.to_dict()
        peptide_hits = [str(i) for i in hits_dict.keys()]
        #only consider those overlapped hits
        sub_dict = dict([(i, overlap_dict[i]) for i in peptide_hits if i in overlap_dict])
        #Generated a subgraph of relationships between all sample hits
        G = nx.Graph(sub_dict)
        #remove non-hit nodes
        for i in G.nodes():
            if i not in peptide_hits:
                G.remove_node(i)
        #Add peptides that have no connections to others
        for i in peptide_hits:
            if i not in G.nodes():
                G.add_node(i)
        #Add z-scores as attributes to the graph
        num_hits = len(hits_series)
        zscore = hits_series.to_dict()
        nx.set_node_attributes(G,'Z-Score', zscore)
        zscore = nx.get_node_attributes(G,'Z-Score')
        
        #Reducing graph to max_degree 2
        if len(G.edges()) != 0:
            #dictinary: node~num connections
            degree=G.degree(nbunch=G.nodes())
            degrees = [i for i in degree.values()]
            max_degree = np.max(degrees)
            while max_degree > 2:
                degree=G.degree(nbunch=G.nodes())
                vertices = [i for i in degree.keys()]
                degrees = [i for i in degree.values()]
                vertices = np.array(vertices); degrees = np.array(degrees);
                max_degree = np.max(degrees)
                #Remove nodes of highest degree (look for lowest z-score if tie)
                if max_degree > 2:
                    max_vertex_indices = np.where(degrees==max_degree)
                    max_vertices = vertices[max_vertex_indices]
                    max_degree_scores = [zscore[i] for i in max_vertices]
                    min_z = min(max_degree_scores)
                    for i in max_vertices:
                        if zscore[i] == min_z:
                            G.remove_node(i)
                            break #so that multiple nodes are not removed
                
        #Eliminates one vertex from each cycle (lowest z-score) to convert them into paths
        len_cycles = 1
        #While loop to make sure that there are no cycles left
        while len_cycles != 0:
            cycles = []
            for i in G.nodes():
                try:
                    #nx.find_cycle() return edges making the cycle
                    #[(A,B),(A,C),(B,c)]
                    cycles.append(nx.find_cycle(G, source=i))
                except:
                    pass
            len_cycles = len(cycles)
            if len_cycles != 0:
                cycles = [np.array(i) for i in cycles]
                cycles = [np.ndarray.flatten(i) for i in cycles]
                cycles = [np.unique(i) for i in cycles]
                cycles = pd.DataFrame(cycles).drop_duplicates().values.tolist()
                for i in range(len(cycles)):
                    cycles[i] = [str(j) for j in cycles[i] if str(j) in G.nodes()]
                node_zscores = nx.get_node_attributes(G, 'Z-Score')
                for i in cycles: # i is list of unique hit nodes representing a cycle
                    cycle_scores=[node_zscores[k] for k in i]
                    min_score = min(cycle_scores)
                    for j in i:
                        if node_zscores[j] == min_score:
                            G.remove_node(j)
                            break #otherwise it will eliminate two nodes in one cycle which both have the same z-score
        
        #Code for deleting vertices from paths based on even or odd length
        node_zscores = nx.get_node_attributes(G,'Z-Score')
        degree=G.degree(nbunch=G.nodes())
        if len(G.edges()) != 0:
            degrees = [i for i in degree.values()]
            max_degree=max(degrees)
            while(max_degree > 0):
                components = list(nx.connected_component_subgraphs(G))
                for i in components:
                    #even paths
                    if len(i.nodes())%2 == 0:
                        path_scores = [node_zscores[k] for k in i]
                        min_score = min(path_scores)
                        for j in i.nodes():
                            if node_zscores[j] == min_score:
                                G.remove_node(j)
                                break #same as above
                    #odd paths
                    elif len(i.nodes())%2 == 1 and len(i.nodes()) != 1:
                        endpoints=[]
                        for j in i.nodes():
                            if degree[j] == 1:
                                endpoints.append(j)
                        if len(endpoints)==2:
                            path = nx.shortest_path(G, source=endpoints[0], target=endpoints[1])
                            middle_indices = np.arange(1,len(path)-1,2)
                            for i in middle_indices:
                                G.remove_node(path[i])
                degree=G.degree(nbunch=G.nodes())
                degrees = [i for i in degree.values()]
                max_degree = max(degrees)
        
        num_nodes = len(G.nodes())
        proportion = np.divide(float(num_nodes), float(num_hits))
        print("Proportion of peptides kept: " + str(proportion))
        
        #Creating pd.Series object from the nodes of the graph
        ind_hits_dict = dict([(i, node_zscores[i]) for i in G.nodes()])
        ind_hits_series = pd.Series(ind_hits_dict)
        #print("Number of peptides kept: " + str(len(ind_hits_series)))  
        
        end_time = timeit.default_timer()
        sample_time = end_time - start_time
        print("Time it took to remove overlaps: {}".format(sample_time))
        
        #print("Number of edges in G: " + str(len(G.edges()))) (should be 0, always)
        
        return ind_hits_series

#end