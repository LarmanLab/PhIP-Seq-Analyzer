# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 15:08:55 2015

@author: yuan
"""

#standard modules
import re

#personal modules
#import myDict
#import myIO
#import myList

#class: sequence

####################################################
class sequence:
    def __init__(self, seq, end3_seq=None):
        #capitalize 
        self.seq = str(seq).upper()
        #remove non alpha characters
        self.seq = re.sub("[^A-Z]", "", self.seq)
        self.seq_len = self.seq.__len__()
        if end3_seq is not None:
            self.adapter3 = end3_seq
            self.len_adapter3 = self.adapter3.__len__()
        #Isoleucine:I; Leucine:L, Valine:V, Phenylalanine: F
        #Methionine:M 	, Cysteine:C, Alanine:A, Glycine:G
        #Proline:P, Threonine:T, Serine:	S, Tyrosine:Y
        #ptophan:W, Glutamine:Q, Asparagine:N, Histidine:H
        #Glutamic acid:E, Aspartic acid:D, Lysine:K, Arginine:R
        #Stop codons 	Stop 	TAA, TAG, TGA 
        self.DNA_codons = {'ATT':'I', 'ATC':'I', 'ATA':'I',\
        'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'TTA':'L', 'TTG':'L',\
        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',\
        'TTT':'F', 'TTC':'F', 'ATG':'M', 'TGT':'C', 'TGC':'C',\
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',\
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',\
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',\
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',\
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',\
        'TAT':'Y', 'TAC':'Y', 'TGG':'W', 'CAA':'Q', 'CAG':'Q',\
        'AAT':'N', 'AAC':'N', 'CAT':'H', 'CAC':'H',\
        'GAA':'E', 'GAG':'E', 'GAT':'D', 'GAC':'D', 'AAA':'K', 'AAG':'K',\
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',\
        'TAA':'.', 'TAG':'.', 'TGA':'.'}

#
    def format_DNA(self):
        #capitalize 
        DNA = self.seq.upper()
        #remove characters except A/T/G/C/N
        DNA = re.sub("[^A-Z]", "", DNA)
        #replace with N except A/T/C/G
        DNA = re.sub("[^A|T|C|G]", "N", DNA)
        return DNA
        
#filter DNA sequence
    def filter_DNA(self, infile):
        #read txt file
        in_obj = open(infile, 'rt')
        DNA_seq = [ line.rstrip("\n") for line in in_obj]
        DNA_seq = ''.join(DNA_seq)

        #remove numeric character and white space
        DNA_seq = re.sub("[0-9]|\s", "", DNA_seq)
        #replace not A/T/C/G with N
        DNA_seq = re.sub("[A|T|G|C]", "N", DNA_seq)
        DNA_seq = DNA_seq.upper()
        return DNA_seq

#filter protein sequence
    def filter_protein(self, infile):
        #read txt file
        in_obj = open(infile, 'rt')
        pep_seq = [line.rstrip("\n") for line in in_obj]
        pep_seq = ''.join(DNA_seq)
        
        #remove numeric character and white space
        pep_seq = re.sub("[0-9]|\s","", pep_seq)
        pep_seq = pep_seq.upper()
        return pep_seq

#convert DNA string to Decimal integer for memory saving
    def numDNA(self):
        int_dna = self.seq.upper()
        adict = {'A':'1', 'T':'3', 'G':'2', 'C':'4', 'N':'5'}
        rx = re.compile('|'.join(list(map(re.escape, adict))))
        int_dna = rx.sub(lambda x: adict[x.group(0)], int_dna)
        int_dna = int(int_dna)
        return int_dna
        
#convert numeric DNA to DNA string
    def reverse_numDNA(self):
        adict = {'1':'A', '3':'T', '2':'G', '4':'C', '5':'N'}
        rx = re.compile('|'.join(list(map(re.escape, adict))))
        dna_str = rx.sub(lambda x: adict[x.group(0)], self.seq)
        return dna_str

#reversed complemented DNA    
    def revcom_DNA(self):
        #print 'Input DNA sequence:', self.seq
        DNA = self.format_DNA()
        #print 'format DNA sequence:', DNA
        #
        if DNA == '':
            print('Error: No DNA sequence input!')
        else:
            #reverse
            rev_seq = DNA[::-1]
            #print rev_seq
            #compliment
            rdict = {'A':'T','T':'A','G':'C','C':'G'}
            robj = re.compile('|'.join(rdict.keys()))
            revcom_DNA = robj.sub(lambda m: rdict[m.group(0)], rev_seq)
            #print rev_com
            
        return revcom_DNA
        
    #translate DNA
    def translate_DNA(self, de=0):
        #print 'Input DNA sequence:', self.seq
        DNA = self.format_DNA()
        
        #slice DNA sequence
        #de=0,1,2
        n = de
        aa = []
        while n < len(DNA):
            coden = DNA[n:n+3]
            if coden in self.DNA_codons.keys():
                aa.append(self.DNA_codons[coden])
            else:
                aa.append('X')
            n += 3
        aa = ''.join(aa)
        return aa
        
#GC content
    def GC_perc(self, digits=1):
        base_len = float(len(self.seq))
        G = self.seq.count('G')
        C = self.seq.count('C')
        #
        GC_ratio = round((G+C)*100/base_len, digits)
        #print GC_ratio
        GC_str = str(GC_ratio) + '%'
        return GC_ratio, GC_str

#calculate kmer
    def kmer(self, kmer_size=4): 
        kmer_counting = {}
        for start in range(self.seq_len-(kmer_size-1)): 
            #print start
            kmer = self.seq[start:start+kmer_size] 
            if kmer in kmer_counting:
                kmer_counting[kmer] += 1
            else:
                kmer_counting[kmer] = 1
            #print(kmer) 
        return kmer_counting

    #positions of enzyme sites
    def site_pos(self, enzyme, start=None, end=None):
        sites = []
        if start is None:
            start = 0
        if end is None or start >= end:
            end = self.seq_len - 1
        #read enzyme sites into list
        for m in re.finditer(enzyme, self.seq):
            pos = m.start() + 1
            if start <= pos <= end:
                sites.append(pos)
        return sites, start, end
    
    def fragmentation(self, enzyme):
        sites,start, end = self.site_pos(enzyme)
        #
        fragments = []
        sites_num = len(sites) - 1
        if sites[0] > start + 1:
            first_site = sites[0] - 1
            fragments.append(self.seq[start:first_site])
        for i in range(sites_num):
            site_start = sites[i] - 1
            site_end = sites[i+1] - 1
            fragments.append(self.seq[site_start:site_end])
        if sites[-1] < end:
            last_site = sites[-1] - 1
            fragments.append(self.seq[last_site:end])
        return sites,fragments

    def all_fragments(self, enzyme):
        #get enzyme sits and min fragments
        sites, fragments = self.fragmentation(enzyme)
            
        #all possible fragments including partial digestion 
        all_fragments = {}
        for in_sites in range(1,len(sites)):
            start_index = 0
            end_index = in_sites
            while end_index <= len(sites):
                frag = "".join(fragments[start_index:end_index])
                all_fragments[frag] = len(frag)
                #print '===', in_sites, start_index, end_index, frag                
                start_index += 1
                end_index += 1
        return all_fragments
        
#trim 3-end adaptor
    def trim_3end(self, seq):
        #print seq
        flag = 0
        #all adapter sequence
        index = re.search(self.adapter3, seq)
        if index is not None:
            trimmed_seq = seq[0:index.start()]
            flag = 1
        else:
            #partial adapter at the 3 end of seq
            for i in range(6, self.len_adapter3+1)[::-1]:
                part_adapter_seq = self.adapter3[0:i]
                #print part_adapter_seq
                index = re.search(r''+self.adapter3+r'\b', seq)
                if index is not None:
                    trimmed_seq = seq[0:index.start()]
                    #print '=', trimmed_seq
                    flag = 1
                    break
        #
        if flag == 1:
            seq = trimmed_seq
        #print(seq)
        return seq  
        
#amino acid sequence        
    def residue_percentage(self, query_aa):
        aa_counts = {}
        aa_perc = {}
        for aa in query_aa:
            aa.upper()
            target = re.findall(aa, self.seq)
            aa_counts[aa] = 0 if target == [] else len(target)
            aa_perc[aa] = aa_counts[aa]*100/float(self.seq_len)
            aa_perc[aa] = round(aa_perc[aa], 2)
        return aa_counts, aa_perc

#aa stretch
    def seq_walking(self, stretch_len, step=1):
        if stretch_len < 2 or stretch_len > len(self.seq):
            stretch_len = 3
        #
        aa_arr = []
        for i in range(0, len(self.seq)-stretch_len+1, step):
            aa_arr.append(self.seq[i:i+stretch_len])
        #
        return (aa_arr)


      


#end