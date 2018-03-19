# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 17:41:09 2015

@author: yuan
"""

#methods related to sequence alignment
import gzip
import os
import re
import subprocess
#personal modules
import myDict
import myGenome
import myIO
import myList
#class: alignment

#
class alignment:

    def __init__(self, parameters=None):
        #parameters is a dictionary
        self.par=parameters
        #print self.par['tool_alinger'] #namely bowtie1
        #print self.par['dir_aligner'] #namely /home/yuan/eRNA/bowtie1/
        #print self.par['file_ref_fa'] #namely /home/yuan/eRNA/ref_seq/human.fa
        self.init_aligner_par()
        
    def init_aligner_par(self):
        if self.par['tool_aligner'] == 'bowtie1':
            self.par['bowtie_aligner'] = self.par['dir_aligner'] + 'bowtie'
            self.par['bowtie_builder'] = self.par['dir_aligner'] + 'bowtie-build'
        elif self.par['tool_aligner'] == 'bowtie2':
            self.par['bowtie_aligner'] = self.par['dir_aligner'] + 'bowtie2'
            self.par['bowtie_builder'] = self.par['dir_aligner'] + 'bowtie2-build'
        #print self.par['bowtie_aligner'], self.par['bowtie_builder']

        #bowtie index
        self.par['bowtie_index_name'] = myIO.file_os(self.par['file_ref_fa']).file_prefix()
        self.par['bowtie_index'] = self.par['dir_aligner'] + self.par['bowtie_index_name']
        #print self.par['bowtie_index'], self.par['bowtie_index_name']
        
    def check_bowtie1(self, tail='.ebwt'):
        #print 'check bowtie1 index of %s' % (self.par['file_ref_fa'])
        #file tail of the bowtie index files
        index_tails = ['.1', '.2','.3','.4','.rev.1','.rev.2']
        index_tails = [ o+tail for o in index_tails]
        #	
        flag=1
        for tail in index_tails:
            index_file = ''.join([self.par['genome_index'],tail])
            if not os.path.isfile(index_file):
                flag = 0
                break
        return flag

    def check_bowtie2(self):
        tail = '.bt2'
        flag = self.check_bowtie1(tail)
        return flag

    def build_bowtie_index(self):
        #check bowtie index
        flag = 1
        if self.par['tool_aligner'] == 'bowtie1':
            flag = self.check_bowtie1()
        elif self.par['tool_aligner'] == 'bowtie2':
            flag = self.check_bowtie2()
        #build bowtie index    
        output = 'Bowtie index of %s exists' % self.par['file_ref_fa']
        if flag == 0:
            print('Build bowtie index {} based on {}.'.format(self.par['genome_index'], self.par['file_ref_fa']))
            shell_command = '{} {} {}'.format(self.par['bowtie_builder'], self.par['file_ref_fa'], self.par['genome_index'])
            print("@@@@@@@@@@@@", shell_command)
            output = subprocess.Popen(shell_command, stdout=subprocess.PIPE, shell=True).stdout.read()
            #output = output.split("\n")
        print(output)
        return output
    
    def bowtie1_alignment(self):
        print('Sequence alignment of ', self.par['sample_name'])
        shell_command = "{} {} -q \"{}\" -S \"{}\"".format(self.par['aligner_options'], \
            self.par['genome_index'], self.par['sample_raw_files'], self.par['sample_sam_file'])
        print("@@@@@@@@@@@@", shell_command)
        #print shell_command
        output = subprocess.Popen(shell_command, stdout=subprocess.PIPE, shell=True).stdout.read()
        #output=output.split("\n")
        
        #compress sam file
        #-f force overwrite of output file and compress links
        shell_command="gzip -f \"{}\"".format(self.par['sample_sam_file'])
        #print "@@@@@@@@@@@@", shell_command
        output = subprocess.Popen(shell_command, stdout=subprocess.PIPE, shell=True).stdout.read()

#read sam file and count reads
#Note: counting by ref    
    def count_reads(self):
        #key is ref name, value is reads string sep by comma, the first is ref seq
        unique_seq = dict((a, []) for a in self.par['ref_dict'].keys()) 
        #unique and multiple counts in dict
        unique = {} #key is ref name, value is counts
        multiple = {} # key is query name, value is the list of refs
        num = {}# counts statistics
        saturation = {0:{1:0, 5:0, 10:0, 'max':0 }} # count number for saturation analysis
        last_index = 0
        
        print('\tread sam file: {}.gz'.format(self.par['sample_sam_file']))
        IN = gzip.open(self.par['sample_sam_file']+'.gz', 'rt')
        UN = gzip.open(self.par['sample_dir']+self.par['sample_name']+'_unknown.fa.gz', 'wt')
        maxRC = 0
        for line in IN:
            #print(line)
            #counts
            num['raw_reads_num'] = num.setdefault('raw_reads_num',0)+1
            #analyze sam line
            info = self.analyze_SAM(line)
            qname, ref= info['qname'], info['ref']
            #unique alignment
            if info['aligned'] == '1':
                unique[ref] = unique.setdefault(ref,0) + 1
                if unique[ref] > maxRC: maxRC = unique[ref]
                #counting of saturation
                if unique[ref] in [1,5,10]:
                    last_counts = saturation[last_index].copy()# copy() is essential!!!!!
                    last_counts[unique[ref]] += 1
                    last_counts['max'] = maxRC#the maximum RC at the time of raw reads we get
                    saturation[num['raw_reads_num']] = last_counts
                    #print num['raw_reads_num'], last_index, saturation[num['raw_reads_num']]
                    last_index = num['raw_reads_num']
                #export aligned sequences of reads
                unique_seq[ref].append(info['seq'])
                num['unique_aligned_reads_num'] = num.setdefault('unique_aligned_reads_num',0)+1
            #multiple alignment
            elif info['aligned'] == '3':
                multiple[qname] = multiple[qname] + [ref] if qname in multiple else [ref]
                num['multialigned_reads_num'] = num.setdefault('multialigned_reads_num',0)+1
            #unalignment
            else:
                UN.write('>'+qname+'\n'+info['seq']+'\n')
                num['unaligned_reads_num'] = num.setdefault('unaligned_reads_num',0) + 1
        IN.close()
        UN.close()
        #counting of saturation
        if num['raw_reads_num'] > last_index:
            saturation[num['raw_reads_num']] = saturation[last_index].copy()
        #for key in sorted(saturation.keys()):
        #    print key, saturation[key]
        
        #upate num statistics
        myIO.file_os(self.par['sample_log'], '=').line_add(num)
        
        print('\tcombine RCs from unique and multiple alignments of ', self.par['sample_name'])
        #reversed multiple
        #print multiple
        rev_multiple = myDict.basic(multiple).counting_reversed_dict()
        #print unique
        RC_dict = self.multiple_counts(unique, rev_multiple)
        #export
        print('\tSave read counts into ', self.par['sample_RC_file'])
        myDict.basic(RC_dict).dict2_to_file(self.par['sample_RC_file'], pattern='\t')
        myDict.basic(saturation).dict2_to_file(self.par['sample_saturation_file'], pattern='\t')
        #
        seq_counts = {}
        for ref, reads_list in unique_seq.items():
            key=ref+'\t'+self.par['ref_dict'][ref]+'\t'+str(len(reads_list))
            if len(reads_list)>0:
                freq_dict = myList.basic(reads_list).elements_frequency0()
                seq_counts[key] = ';'.join(str(a)+':'+str(b) for a,b in freq_dict.items())
            else:
                seq_counts[key] = 'NA'
        myDict.basic(seq_counts).dict_to_file(self.par['sample_dir']+'unique_aligned_reads.txt', pattern='\t')
        
#initiate RC_dict        
    def init_RCdict(self):
        RC_dict = {}
        #get all ref names from the refseq file
        ref_names = myGenome.genome(self.par['file_ref_fa']).fa_displayid()
        for ref in ref_names:
            RC_dict[ref] = {'lowRC':0, 'midRC':0, 'highRC':0}
        #
        return RC_dict

#            
    def multiple_counts(self, unique, multiple):
        #init RC dict
        RC_dict = self.init_RCdict()
        #
        for ref, uRC in unique.items():
            RC_dict[ref]['lowRC'] += uRC
            RC_dict[ref]['midRC'] += uRC
            RC_dict[ref]['highRC'] += uRC

        #
        uniq_total = 0
        for ref_name_str, mRC in multiple.items():
            #print ref_name_str,mRC
            ref_names = ref_name_str.split(',')
            #assign multiple counts based on unique ratio
            mcounts = {}
            #uni_counts = 0
            for ref in ref_names:
                mcounts[ref] = unique[ref] if ref in unique else 1
                uniq_total += mcounts[ref]
            #assign formula
            mcounts = dict([(ref, float(RC)*mRC/uniq_total) for ref, RC in mcounts.items()])
            #update RC_dict
            for ref in ref_names:
                RC_dict[ref]['midRC'] += mcounts[ref]
                RC_dict[ref]['highRC'] += multiple[ref]
              
        return RC_dict

#extract info from one record line in a sam file
    def analyze_SAM(self, sam_line):
        #loop of the sam line
        items = sam_line.split("\t")
        sam_info = {'qname':items[0],'FLAG':items[1],'ref':items[2], \
                  'offset':items[3],'mapQ':items[4],'CIGAR':items[5], \
                  'mate_ref':items[6],'mate_offset':items[7], \
                  'fragment_size':items[8],'seq':items[9], 'Qual':items[10] }
        #optional fields
        for item in items[11:]:
            (tag, mid, value) = item.split(':')
            sam_info[tag] = value
            #print key,value
        
        #analyze sam info
        #coordinate from fastq file
        qname_list = sam_info['qname'].split(':')   
        if len(qname_list) > 3 and qname_list[-2].isdigit() and qname_list[-1].isdigit():
            sam_info['coordinate'] = '-'.join(qname_list[-2:])
        else:
            sam_info['coordinate'] = 'NA'
        #print qname_list, sam_info['coordinate']
        
        #analyze FLAG
        flag_info = self.analyze_FLAG(sam_info['FLAG'])
        sam_info.update(flag_info)
        
        #analyze CIGAR
        cigar_info = self.analyze_CIGAR(sam_info['CIGAR'], sam_info['offset'])
        sam_info.update(cigar_info)
        
        return sam_info

    def analyze_FLAG(self, FLAG):
        binary_str = format(int(FLAG), 'b')
        #fixed length 12 characters
        binary_str = ''.join(['0']*(12-len(binary_str))) + binary_str
        #print sam_info['FLAG'],binary_str
        #
        flags = {}
        #paired-end, single-end
        flags['paired'] = binary_str[-1]
        #aligned, unaligned, multiple-aligned
        if binary_str[-3] == '1':
            flags['aligned'] = '0' #unaligned
        elif binary_str[-9] == '1':
            flags['aligned'] = '3' #multiple-alinged
        else:
            flags['aligned'] = '1'
        #reversed complementary
        flags['revcom'] = binary_str[-5]
        #R1 or R2
        if binary_str[-8] == 1:
            flags['end'] = 'R2' #unaligned
        else:
            flags['end'] = 'R1'
        #print FLAG, flags['aligned']
        return flags
        
    def analyze_CIGAR(self, CIGAR, offset):
        cigar_info = {}
        exons = [int(offset)]
        insertions = []
        deletions = []
        #analyze CIGAR
        cigar_num = re.split(r'[A-Z]',CIGAR)[:-1]
        cigar_chr = re.findall(r'[A-Z]', CIGAR)
        for a, b in zip(cigar_num, cigar_chr):
            a = int(a)
            if b == 'M':
                if len(exons)%2 == 1: #matches
                    exons.append(exons[-1]+a-1)
                else:
                    exons[-1] += a
            elif b == 'N': # intron
                exons.append(exons[-1]+a+1)
            elif b == 'D': #deletion
                exons[-1] += a
                deletions.append(exons[-1]+1)
                deletions.append(exons[-1]+a)
            elif b == 'I': #insertion
                insertions.append(exons[-1])
        #print offset, exons
        #
        #for start, end in zip(exons[::2],exons[1::2]):
        #    print start, end
        cigar_info['exons'] = exons
        cigar_info['insertions'] = 'NA' if insertions == [] else insertions 
        cigar_info['deletions'] = 'NA' if deletions == [] else deletions 
        #print cigar_info
        return cigar_info
    
#combine counting files into a matrix
    def combine_countfiles(self, args_tuple):
        #row_names should be None or list type
        infile_tail, RC_level, out_file, row_names = args_tuple
        #
        counting_dict2 = {}
        for sample_name in self.par['sample_names']:
            #get read counts of a given sample
            counting_file = '{}{}/{}{}'.format(self.par['dir_result'], sample_name, sample_name, infile_tail)
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





