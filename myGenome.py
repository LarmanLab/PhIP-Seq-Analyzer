# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 13:59:50 2015

@author: yuan
"""


#standard modules
import gzip
import itertools
import os
import re
import string
import numpy as np

#personal modules
import myDict
import myIO
import myList

#
#class genome
    
##############################################
#class 
class genome:
    def __init__(self, biofile=None, sep=None):
        self.biofile = biofile
        #seperate character
        if sep is None:
            self.sep = ',' if myIO.file_os(self.biofile).name_suffix() == 'csv' else "\t"
        else:
            self.sep = sep
        self.record_num = 0

#ouput is read-only file handle
    def readonly_handle(self, infile):
        if infile.endswith('.gz'):
            in_obj = gzip.open(infile, 'rt')
        else:
            in_obj = open(infile, 'rt')
        return in_obj
    
#fasta title:
#Ensembl example: MT dna:chromosome chromosome:GRCh38:MT:1:16569:1 REF
    def acc_number(self, line):
        line = line.rstrip()
        line = re.sub(r"^>", "", line)
        if len(re.findall(r'\|', line)) > 2:
            items = line.split('|')
            acc = items[1]
        elif len(re.findall(r' ', line)) > 2:
            items = line.split()
            acc = items[0]
        else:
            acc = line
        #print acc
        return acc
        
    #return is the dict
    def read_fa(self):
        #store accession: sequences
        ref_dict = {}
        ref_list = []
        in_obj=self.readonly_handle(self.biofile)
        #read fa file
        for line in in_obj:
            line = line.rstrip()
            if line.find(">") == 0:
                acc = self.acc_number(line)
                ref_list.append(acc)
                ref_dict[acc] = []
                self.record_num += 1
            else:
                ref_dict[acc].append(line)
        in_obj.close()
        #arrage ref
        for fa_id in ref_list:
            #print fa_id, ref_dict[fa_id][:3], "\n"
            ref_dict[fa_id] = ''.join(ref_dict[fa_id])
            #print fa_id, len(re.findall('^N*', ref_dict[fa_id])[0]), len(ref_dict[fa_id]),'\n'
        #print ref_list
        #print ref_dict
        return ref_dict, ref_list
 
#return is names of ref in fasta
    def fa_displayid(self):
        #store accession: sequences
        ref_id = []
        in_obj=self.readonly_handle(self.biofile)
        for line in in_obj:
            line = line.rstrip(string.whitespace)
            if line.find(">") == 0:
                ref_id.append(line[1:])
        in_obj.close()
        #print
        return ref_id
        
    #return is the first sequence
    def fa_first(self):
        num = 0
        seq = []
        in_obj=self.readonly_handle(self.biofile)
        for line in in_obj:
            line = line.rstrip()
            if line.find(">") == 0:
                num += 1
                if num > 1:
                    break
            else:
                seq.append(line)
        in_obj.close()
        #print
        seq = ''.join(seq)
        #print seq
        return seq

    def combine_fa(self, fa_dict):
        chrs = myList.basic(fa_dict.keys()).sort_list()
        out_obj = open(self.biofile, 'wt')
        for chr_name in chrs:
            fa_file = fa_dict[chr_name]
            seq=genome(fa_file).read_fa_first()
            out_obj.write('>{}\n{}\n'.format(chr_name, seq))
            print('\t{}:{}'.format(chr_name, fa_file))
        out_obj.close()
        print('Combine fa files into ', self.biofile)
        
    #read fasta record one by one with functional process
    def fa_loop(self, FUN=None, *args):
        #store accession: sequences
        ref = {}
        with open(self.biofile, 'rt') as in_obj:
            for L1,L2 in zip(*[in_obj]*2):
                self.record_num += 1
                #first line
                L1 = L1.rstrip()
                acc = self.acc_number(L1)
                #second line
                L2 = L2.rstrip()
                #ref information
                ref[acc] = L2 if FUN is None else FUN(L2, *args)
                #print ref[acc]
        in_obj.close()
        return ref
        
    #read fastq record one by one with functional processing
    def fq_loop(self, FUN=None, *args):
        #store accession: sequences
        reads = {}
        readsQ = {}
        with open(self.biofile, 'rt') as in_obj:
            for L1, L2, L3, L4 in zip(*[in_obj]*4):
                self.record_num += 1
                #first line
                L1 = L1.rstrip()
                items = L1.split()[0]
                items = items.split(':')
                coordinate = items[-2]+'-'+items[-1]
                #print items, coordinate
                
                #second line: sequence
                seq = L2.rstrip()
                #print seq
                reads[coordinate] = seq if FUN is None else FUN(seq, *args)
                    
                #fourth line:  Q value
                L4 = L4.rstrip()
                Q = [ord(x) for x in L4] #ASCII to decimal
                #print Q
                readsQ[coordinate] = Q
        in_obj.close()
        return reads, readsQ
    
    #sequencing direction: R1 or R2
    def R1R2(self):
        #print(self.biofile)
        first_line = myIO.file_os(self.biofile).first_line()
        direction = 'R2' if re.search(r'\/2$', first_line) else 'R1'
        #print direction
        return direction
        
#split fastq file determined by Illumina analyzer
#based on barcode
    def demultiplex_fq(self, par):
        #our directory
        out_dir = myIO.dir_os(par['dir_raw_data']).create_dir()
        #sequencing direction: R1 or R2
        direction = self.R1R2()
        #read relationship between barcode vs sample from sample_file
        barcode_sample = myIO.file_os(par['barcode_file'], '\t').to_dict()
        #barcode_sample={ mySequence.sequence(k).revcom_DNA():v for k,v in barcode_sample.items()}
        barcode_sample['unassigned'] = 'unassigned'
        #print barcode_sample
        #open file handles based on barcode_sample
        file_handle = {}
        barcode_file = {}
        known_dict = {}
        un_dict = {}
        for barcode, sample_name in barcode_sample.items():
            fq_file = '{}{}_{}.fq'.format(out_dir, sample_name, direction)
            file_handle[barcode] = open(fq_file, 'wt')
            barcode_file[barcode] = fq_file
            known_dict[barcode] = {'sample_name':sample_name,'read_counts':0}
        ###

        #file handle
        #with open(self.biofile, 'rt') as F1, open(index_file, 'rt') as F2:
        F1=self.readonly_handle(self.biofile)
        F2=self.readonly_handle(par['index_file'])
        n = 0 #total number of reads
        m = 0 # total number assigned reads
        stdout_format = '|{:^15}|{:^15}|{:^15}|{:^15}|'
        dash_line=stdout_format.format('-'*15, '-'*15, '-'*15, '-'*15)
        print(dash_line)
        print(stdout_format.format('Raw reads','Assigned reads', 'Percentage', 'Read Length'))
        print(stdout_format.format('millions','millions', '%', 'nt'))
        print(dash_line)
        with F1, F2:
            #read 4 lines at a time per file
            for L1,La, L2,Lb, L3,Lc, L4,Ld in itertools.zip_longest(*[F1,F2]*4):
                barcode = Lb.rstrip()
                #assign record based on barcode
                if barcode in file_handle and len(L2) >= par['seq_min']:
                    L_name = re.sub(r'\/', '#'+barcode+'/', L1)
                    #print L_name, La
                    #trim reads from 5 end
                    if par['seq_start'] > 0:
                        L2 = L2[par['seq_start']:]
                        L4 = L4[par['seq_start']:]
                    #trim the longer reads from 3-end
                    if par['seq_end'] != 0:
                        L2 = L2.rstrip()
                        L4 = L4.rstrip()
                        L2 = L2[:par['seq_end']]+"\n"
                        L4 = L4[:par['seq_end']]+"\n"
                    #output file handle
                    file_handle[barcode].writelines([L_name,L2,L3,L4])
                    #counting
                    known_dict[barcode]['read_counts'] += 1
                    m += 1
                else:
                    #output file handle
                    file_handle['unassigned'].writelines([L1,L2,L3,L4])
                    un_dict[barcode] = un_dict[barcode] + 1 if barcode in un_dict else 1
                    known_dict['unassigned']['read_counts'] += 1
                n += 1
                #export when 
                if m >= 1e6 and m%1e6 == 0: #million
                    print(stdout_format.format(n/1e6, m/1e6, round(m*100/n, 2), len(L2)-1 ))
                #if n==3e6: break
            else:
                print(dash_line)
                print(stdout_format.format(n/1e6, m/1e6, m*100/n, '---'))
                print(dash_line)
        #calculate percentage
        for bc in known_dict.keys():
            RC = float(known_dict[bc]['read_counts'])
            known_dict[bc]['percentage_%'] = round(RC*100/n, 2)
        #close file handle
        for b, F in file_handle.items():
            #close file handle
            F.close()
            #delete empty file
            if os.stat(barcode_file[b]).st_size == 0:
                os.remove(barcode_file[b])
        #export statistics
        myDict.basic(known_dict).dict2_to_file(out_dir+'known.log', '\t')
        myDict.basic(un_dict).dict_to_file(out_dir+'unknown.log', '\t')
        #no return
        
#cbind two fastq into one fastq
    def cbind_fq(self, infile2, outfile):
        print("Combine index fastq files {} and {} into {}\n".format(self.biofile, infile2, outfile))
        #get file handles of the two fastq files, and the output file
        F1=self.readonly_handle(self.biofile)
        F2=self.readonly_handle(infile2)
        out_obj = open(outfile, 'wt')
        with F1, F2:
            #read 4 lines at a time per file
            for L1, La, L2,Lb, L3,Lc, L4,Ld in zip(*[F1,F2]*4):
                #sequence
                L2 = L2.rstrip()
                Lb = Lb.rstrip()
                L_seq=L2+Lb+"\n"
                #revcom_L2=mySequence.sequence(L2).revcom_DNA()
                #revcom_Lb=mySequence.sequence(Lb).revcom_DNA()
                #L_seq=revcom_L2+revcom_Lb+"\n"
                #Q-scores
                Ld = Ld.rstrip()
                L_Q = L4.rstrip()+Ld[::-1]+"\n"
                #export to the output file                
                out_obj.writelines([L1,L_seq,L3,L_Q])
        
#trim fastq into a new fastq
    def trim_fq(self, outdir, seq_start=0, seq_end=0):
        file_name = myIO.file_os(self.biofile).file_name()
        outfile = outdir + re.sub('\.gz$', '', file_name)
        print("Trim fastq files {}, and save new file {}\n".format(self.biofile, outfile))
        #get file handles of the two fastq files, and the output file
        F1 = self.readonly_handle(self.biofile)
        out_obj = open(outfile, 'wt')
        with F1:
            #read 4 lines at a time per file
            for L1, L2, L3, L4 in zip(*[F1]*4):
                if seq_start > 0:
                    L2 = L2[seq_start:]
                    L4 = L4[seq_start:]
                #trim the longer reads from 3-end
                if seq_end != 0:
                    L2 = L2.rstrip()
                    L4 = L4.rstrip()
                    L2 = L2[:seq_end] + "\n"
                    L4 = L4[:seq_end] + "\n" 
                #export to the output file                
                out_obj.writelines([L1,L2,L3,L4])
                
#read gtf file, get features and export to a text file
#column can be ['chr','start','end','transcript_id','gene_id', 'gene_name', 'gene_biotype']
    def export_annot(self, feature_name, columns, out_file):
        #
        print("get annotation of", feature_name)
        annot = self.read_gtf(feature_name)
        #
        print("export to", out_file)
        out_obj = open(out_file, 'wt')
        #header line
        header = list(columns)
        header.insert(0,'id')
        header = '\t'.join(header)
        out_obj.write(header+'\n')
        #
        for annot_dict in annot.keys():
            arr = [annot_dict]
            for col in columns:
                #print col
                arr.append(annot[annot_dict][col])
            out_obj.write('\t'.join(arr)+'\n')
        out_obj.close()
    
#annotation.txt
#extract annotation information based input name
    def normalize_by_factor(self, annot_file, col_index, col_factor='pro_len'):
        #read RC file
        #RC_df=pd.read_table(self.biofile, sep=self.sep, index_col=True)
        RC_df = myIO.file_os(self.biofile, sep=self.sep).to_df(header=True, rowname=True)
        #print(RC_df.shape)

        #read annotation file
        #annot_df=pd.read_table(annot_file, sep="\t", index_col=True)
        annot_df = myIO.file_os(annot_file, sep="\t").to_df(header=True, rowname=True)
        #print list(annot_df)
        sub_annot = annot_df[[col_index, col_factor]].drop_duplicates([col_index], take_last=True)
        sub_annot.index = list(sub_annot[col_index])
        sub_annot = sub_annot.ix[:,1:] #remove the column with col_index
        #print sub_annot
        
        #sort annot_df by row names of RC_df
        sub_annot = sub_annot.ix[list(RC_df.index)]
        #for missing proteins, pro_len equal ave-pro_len
        ave_pro_len = np.mean(sub_annot[col_factor])
        pro_len_df = sub_annot.fillna(ave_pro_len)
        #print(pro_len_df.shape)
        
        #normalization by aa length of proteins
        RC_df.insert(0, col_factor, list(pro_len_df[col_factor]) )
        normRC_df = RC_df.apply(lambda x: x[1:]/x[0], axis=1)
        #scaling normalization by million reads
        def norm_func(x):
            sum_x = np.sum(x)
            norm_x = x*10e6/sum_x if sum_x >0 else x
            norm_x = np.round(norm_x)
            norm_x = norm_x.astype(int)
            return norm_x
        normRC_df = normRC_df.apply(norm_func, axis=0)
        #print normRC_df
        
        return normRC_df
        
#read gtf/gff3 file into a dictionary
#feature may be one of 'gene', 'exon', 'transcript', 'chromosome'
    def read_ncbi_gff(self, feature_key='seqid', feature_value='chromosome'):
        feature_annot = {}
        #store accession: sequences
        in_obj = self.readonly_handle(self.biofile)
        for line in in_obj:
            if not line.startswith('#'):
                items = line.rstrip(';\n').split("\t")
                #return one recorder and save into a dictionary
                annot = {}
                annot['seqid'], annot['source'], annot['feature']=items[0], items[1], items[2]
                annot['start'], annot['end'], annot['strand']=items[3], items[4], items[6]
                #print annot['source']'
                #treat dbxref
                attr = items[8]
                s_obj = re.search(r"Dbxref=(.*?);", attr)
                if s_obj:
                    dbxref = re.sub(r"HGNC=HGNC=", "HGNC=HGNC:", s_obj.groups()[0])
                    dbxref = re.sub(r"MGI=MGI=","MGI=MGI:", dbxref)
                    dbxref = re.sub(r":","=", dbxref)
                    dbxref = re.sub(r",",";", dbxref)
                    #print dbxref
                    attr = attr[:s_obj.span()[0]] + dbxref + ';' + attr[s_obj.span()[1]:]
                #
                attr = attr.split(";")
                #print attr
                for a in attr:
                    try:
                        name, value = a.split("=")
                        annot[name] = value
                    except ValueError:
                        print(attr)
                    #print name, annot[name], "\n"
                print(annot)
                #update feature_annot
                if feature_key and feature_value in annot:
                    feature_annot[annot[feature_key]] = annot[feature_value]
                    print(annot[feature_key], annot[feature_value])
                    self.record_num += 1
                #break
        in_obj.close()
        return feature_annot
        
#identical ids between fasta and gtf in order    
#self.biofile is .fa file
    def match_ncbi_fa(self, gtf_file):
        #read matched information from gtf_file
        seqid_chr = genome(gtf_file).read_ncbi_gff('seqid','chromosome')
        #out_gtf_file
        file_head = myIO.file_os(self.biofile).file_prefix()
        out_gtf = file_head + '.' + myIO.file_os(gtf_file).name_suffix()
        out_obj = open(out_gtf, 'wt')
        #get order of chr
        ref_arr = self.fa_displayid()

        #export matched gtf
        print('Match the first column of {} with {}, => {}'.format(gtf_file, self.biofile, out_gtf))
        for chr_id in ref_arr:
            chromosome = re.sub('chr', '', chr_id,flags=re.IGNORECASE)
            n = 0
            #read gtf            
            in_obj = self.readonly_handle(self.biofile)
            for line in in_obj:
                if not line.startswith('#'):
                    items = line.split("\t")
                    seqid = items[0]
                    if seqid in seqid_chr and seqid_chr[seqid] == chromosome:
                        items[0] = chr_id
                        myline = "\t".join(items)
                        out_obj.write(myline)
                        n += 1
            in_obj.close()
            print('{}:{}'.format(chr_id, n))
        out_obj.close()     

#read gtf/gff3 file into a dictionary
#feature may be one of 'gene', 'exon', 'transcript', 'chromosome'
    def read_ensembl_gff(self, feature_key='seqid', feature_value='gene_id'):
        feature_annot = {}
        #store accession: sequences
        in_obj = self.readonly_handle(self.biofile)
        for line in in_obj:
            if not line.startswith('#'):
                line = line.rstrip(';\n')
                items = line.split("\t")
                #return one recorder and save into a dictionary
                annot = { 'seqid':items[0], 'source':items[1], 'feature':items[2],\
                         'start':items[3], 'end':items[4], 'strand':items[6] }
                #print annot['source']'
                attr = re.sub(r"\"", "", items[8])
                attr = attr.split("; ")
                #print attr
                for a in attr:
                    try:
                        name, value = a.split(" ")
                        annot[name] = value
                    except ValueError:
                        print(attr)
                    #print name, annot[name], "\n"
                #print annot
                #update feature_annot
                if feature_key and feature_value in annot:
                    feature_annot[annot[feature_key]] = annot[feature_value]
                    print(annot[feature_key], annot[feature_value])
                    self.record_num += 1
                #break
        in_obj.close()
        return feature_annot  


#identical ids between fasta and gtf in order    
#self.biofile is .fa file
    def match_ensembl_fa(self, gtf_file):
        #get order of chr
        ref_arr = self.fa_displayid()
        #print ref_arr
        
        #export matched gtf
        file_head = myIO.file_os(self.biofile).file_prefix()
        out_gtf = file_head + '.' + myIO.file_os(gtf_file).name_suffix()
        out_obj = open(out_gtf, 'wt')#out_gtf_file
        print('Match the first column of {} with {}, => {}'.format(gtf_file, self.biofile, out_gtf))
        for chr_id in ref_arr:
            n = 0
            #read gtf
            in_obj=self.readonly_handle(self.gtf_file)
            for line in in_obj:
                if not line.startswith('#'):
                    items = line.split("\t")
                    seqid = items[0]
                    if chr_id == seqid:
                        out_obj.write(line)
                        n += 1
            in_obj.close()
            print('{}:{}'.format(chr_id,n))
        out_obj.close()  
        
#get idmapping from UniProt *.dat
#simple mode:always uniprot_acc ~ *
    def uniprot_idmapping_simple1(self, value='Ensembl'):
        id_dict = {}
        #read *.dat file
        in_obj=self.readonly_handle(self.biofile)
        for line in in_obj:
            line = line.rstrip()
            UniProtKB_AC,ID_type,ID=line.split("\t")
            if ID_type == value:
                if UniProtKB_AC in id_dict:
                    id_dict[UniProtKB_AC].append(ID)
                else:
                    id_dict[UniProtKB_AC]=[ID]
        in_obj.close()
        return id_dict

#simple mode:always * ~ uniprot_acc
    def uniprot_idmapping_simple2(self, value='Ensembl'):
        id_dict = {}
        #read *.dat file
        in_obj=self.readonly_handle(self.biofile)
        for line in in_obj:
            line = line.rstrip()
            UniProtKB_AC, ID_type, ID=line.split("\t")
            if ID_type == value:
                if ID in id_dict:
                    id_dict[ID].append(UniProtKB_AC)
                else:
                    id_dict[ID] = [UniProtKB_AC]
        in_obj.close()
        return id_dict
                    
#normal mode:key and value can be one of ID_type
    def uniprot_idmapping(self, type1='KEGG', type2='Ensembl'):
        #retrieve seperately
        key_dict = self.uniprot_idmapping_simple1(type1)
        value_dict = self.uniprot_idmapping_simple1(type2)
        #combine
        id_dict = {}
        for uniprot_acc, key_list in key_dict.items():
            if uniprot_acc in value_dict:
                for key_element in key_list:
                    id_dict[key_element] = value_dict[uniprot_acc]
        return id_dict
            
            
#end