 # -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 13:08:40 2015

@author: yuan
"""
#standard modules
import collections
import gzip
import pandas as pd
import re
import os
import urllib
import urllib2
#personal modules
import myDict
import myList
import myGenome
import myIO
import mySystem


################################################################################
#class web
#class ensembl
#class NCBI
#class annotation
################################################################################
class web:
    def __init__(self, url):
        self.url = url
        #self.url = url if re.search(r'/$', url) else url+'/'
        
    #capture html page
    def get_html(self):
        print('retrieve the linkage:', self.url)
        try:
            url_obj = urllib2.urlopen(self.url)
        except:
            #print self.url
            return None
        else:
            html = url_obj.read()
            lines = html.split('\n')
            #for l in lines:
            #    print l
            return lines
            
    #download file from HTTP  or FTP      
    def download_file(self, local_file):
        print("downloading", self.url)
        testfile = urllib.URLopener()
        try:
            testfile.retrieve(self.url, local_file)
        except:
            print("Error: no", self.url)
            return None
    
    #return dictionary of relative name and absolute name    
    def ls_html(self):    
        #capture html page
        lines = self.get_html()
        #get names (directory or filename)
        dir_names = {}
        file_names = {}
        for l in lines:
            #print l
            items = l.split()
            if len(items) > 0:
                name = items[-1]
                abs_name = self.url+name+'/'
                #file return None
                if web(abs_name).get_html() is None:
                    file_names[self.url+name] = name
                else:
                    dir_names[abs_name] = name
                    #print items[-1]
        #print dir_names,
        #print file_names
        return dir_names, file_names
        
    def recursive_ls(self):
        #get dir and file of top level
        dir_names,file_names = self.ls_html()
        dirs = dir_names.keys()
        #print dirs
        #
        while len(dirs)>0:
            sub_url = dirs.pop()
            subdir_names, subfile_names = web(sub_url).ls_html()
            print(sub_url, subdir_names)
            #
            dir_names.update(subdir_names)
            file_names.update(subfile_names)
            dirs.extend(subdir_names.keys())
        #for i in dir_names:
        #    print i
        return dir_names, file_names
            
            
################################################################################
class annotation:
    #specie: human, mouse
    def __init__(self, specie='human', local_dir='/home/yuan/Downloads/'):
        self.specie = specie
        self.local_dir = local_dir
        
    #default export entrez refseq protein id vs. GO id    
    def GO_annot_human(self, by='RefSeq_NT'):
        #download id-associated file
        id_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/'
        local_id_file = self.local_dir + 'HUMAN_9606_idmapping.dat.gz'
        web_id_file = id_url + 'HUMAN_9606_idmapping.dat.gz'
        #download
        web(web_id_file).download_file(local_id_file)
        #get association between uniprotKB and refseq_nt
        id_assoc = {}
        with gzip.open(local_id_file,'rb') as f:
            for line in f:
                line = line.rstrip("\n")
                (uniprot, annot_type, annot_id) = line.split("\t")
                #print line
                if annot_type == by: 
                    annot_id = annot_id[:-2]
                    #one annot_id and one uniprot, one uniprot might be with multiple annot_ids
                    if uniprot in id_assoc:
                        id_assoc[uniprot].append(annot_id)
                        #print '====', uniprot, id_assoc[uniprot]
                    else:
                        id_assoc[uniprot] = [annot_id]
        #print id_assoc
        
        #download GO-associated file
        go_url = 'http://geneontology.org/gene-associations/'
        local_go_file = self.local_dir + 'gene_association.goa_ref_human.gz'
        web_go_file = go_url + 'gene_association.goa_ref_human.gz'
        web(web_go_file).download_file(local_go_file)
        #get association between uniprotKB and refseq_nt
        go_assoc=collections.defaultdict(dict)
        with gzip.open(local_go_file,'rb') as f:
            for line in f:
                if line.find('!') != 0:
                    line = line.rstrip("\n")
                    items = line.split("\t")
                    uniprot = items[1]
                    pro_symbl = items[2]
                    go_id = items[4]
                    pro_name = items[9]
                    #print go_id, pro_name
                    if uniprot in id_assoc:
                        for annot_id in id_assoc[uniprot]:
                            if annot_id in go_assoc:
                                go_assoc[annot_id]['go_id'] += ',' + go_id
                            else:
                                go_assoc[annot_id]['go_id'] = go_id
                                go_assoc[annot_id]['uniprot_id'] = uniprot
                                go_assoc[annot_id]['uniprot_symbl'] = pro_symbl
                                go_assoc[annot_id]['uniprot_name'] = pro_name
        #export to file
        go_assoc = pd.DataFrame(go_assoc).transpose()
        outfile = self.local_dir + self.specie + '_GO_annotation.txt'
        go_assoc.to_csv(outfile, sep="\t",index_label=by)

    #default export entrez refseq protein id vs. GO id    
    def GO_annot_mouse(self,by='MGI'):
        #download id-associated file
        id_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/'
        local_id_file = self.local_dir + 'MOUSE_10090_idmapping.dat.gz'
        web_id_file = id_url + 'MOUSE_10090_idmapping.dat.gz'
        #download    
        #web(web_id_file).download_file(local_id_file)
        #get association between uniprotKB and refseq_nt
        id_assoc = {}
        with gzip.open(local_id_file,'rb') as f:
            for line in f:
                line = line.rstrip("\n")
                (uniprot, annot_type, annot_id) = line.split("\t")
                #print line
                if annot_type == by: 
                    #one annot_id and one uniprot, one uniprot might be with multiple annot_ids
                    if uniprot in id_assoc:
                        id_assoc[uniprot].append(annot_id)
                        #print('====', uniprot, id_assoc[uniprot])
                    else:
                        id_assoc[uniprot] = [annot_id]
        #print id_assoc
        go_url = 'http://geneontology.org/gene-associations/'
        local_go_file = self.local_dir + 'gene_association.mgi.gz'
        #web_go_file=go_url + 'gene_association.mgi.gz'
        #download GO-associated file
        web(web_go_file).download_file(local_go_file)
        #get association between uniprotKB and refseq_nt
        go_assoc = collections.defaultdict(dict)
        with gzip.open(local_go_file,'rb') as f:
            for line in f:
                if line.find('!') != 0:
                    line = line.rstrip("\n")
                    items = line.split("\t")
                    uniprot = items[1]
                    pro_symbl = items[2]
                    go_id = items[4]
                    pro_name = items[9]
                    #print go_id, pro_name
                    if uniprot in id_assoc:
                        for annot_id in id_assoc[uniprot]:
                            if annot_id in go_assoc:
                                go_assoc[annot_id]['go_id'] += ',' + go_id
                            else:
                                go_assoc[annot_id]['go_id'] = go_id
                                go_assoc[annot_id]['uniprot_id'] = uniprot
                                go_assoc[annot_id]['uniprot_symbl'] = pro_symbl
                                go_assoc[annot_id]['uniprot_name'] = pro_name
        #export to file
        go_assoc = pd.DataFrame(go_assoc).transpose()
        outfile = self.local_dir + self.specie + '_GO_annotation.txt'
        go_assoc.to_csv(outfile, sep="\t",index_label=by)        
        
        
################################################################################
#download the newest genome sequences from Ensembl
class ensembl:
    def __init__(self, specie, out_dir):
        self.specie = specie
        self.out_dir = myIO.dir_os(out_dir+specie).create_dir()
        #initiate url list
        self.url_list()
                
    def url_list(self):
        species = {'human':'homo_sapiens', 'mouse':'mus_musculus',\
                 'rat':'rattus_norvegicus', 'maize':'zea_mays'}
        #annot_type=['dna_fa','cds_fa','pep_fa','ncrna_fa','gff3', 'gtf']
        url = {} 
        #ensembl: human, rat, mouse
        ensembl_specie = ['human', 'rat', 'mouse']
        #ensembl plants
        ensemblgenome_specie = ['maize']
        if self.specie in ensembl_specie:
            url['pub'] = 'ftp://ftp.ensembl.org/pub/'
            #fasta format
            fa_format = url['pub'] + 'current_fasta/' + species[self.specie]+'/'
            #gff3 format
            url['gff'] = url['pub'] + 'current_gff3/' + species[self.specie] + '/'
            #gtf format
            url['gtf'] = url['pub'] + 'current_gtf/' + species[self.specie] + '/'
        elif self.specie in ensemblgenome_specie:
            url['pub'] = 'ftp://ftp.ensemblgenomes.org/pub/current/plants/'
            #fasta format
            fa_format = url['pub'] + 'fasta/' + species[self.specie]+'/'
            #gff3 format
            url['gff'] = url['pub'] + 'gff3/' + species[self.specie] + '/'
            #gtf format
            url['gtf'] = url['pub'] + 'gtf/' + species[self.specie] + '/'
        #
        #newest release version
        #self.ver=self.newest_release(url['pub'])
        url['dna_fa'] = fa_format + 'dna/'
        url['cdna_fa'] = fa_format + 'cdna/'
        url['cds_fa'] = fa_format + 'cds/'
        url['ncrna_fa'] = fa_format + 'ncrna/'
        url['protein'] = fa_format + 'pep/'
        self.url = url
        
    #release version of genome annotation
    def newest_release(self, url):
        lines = web(url).get_html()
        #print html
        new_ver = 1
        for line in lines:
            #print line
            m = re.search(r'release-[0-9]*', line)
            if m:
                ver = line[(m.start()+8):m.end()]
                if new_ver < int(ver):
                    new_ver = int(ver)
        new_ver = 'release-' + str(new_ver)
        print('Release:', new_ver)
        return new_ver
    
    #get DNA choromosome file names
    def dna_files(self, lines):
        chr_files = {}
        for line in lines: 
            if line.find('dna.chromosome') > 0:
                if not line.find('CHR_') > 0:
                    items = line.split()
                    file_name = items[-1]
                    chromosome = re.sub(r'^.*chromosome\.|\.fa.gz', '', file_name)
                    chr_files[chromosome] = file_name 
        #print dict
        #myDict.basic(chr_files).print_dict()
        return chr_files

    def single_file(self, lines):
        chr_files = {}
        for line in lines: 
            if re.search(r'\.gz', line):
                if line.find('abinitio') > 0 or line.find('chromosom') > 0:
                    pass
                else:
                    items = line.split()
                    file_name = items[-1]
                    chr_files['all'] = file_name
                    print(file_name)
        return chr_files
       
    #download genome files        
    def download_annot(self,genome_type):
        #print genome_type
        url = self.url[genome_type]
        #get html
        lines = web(url).get_html()
        #get the list of files
        chr_files = self.single_file(lines)
        print(chr_files)
        #download and decompress genome files
        local_chr_files = {}
        for key in chr_files.keys():
            gz_file = myIO.file_os(url+chr_files[key]).download(self.out_dir)
            #decompress file
            ungz_file = myIO.file_os(gz_file).decompress_gz()
            local_chr_files[key] = ungz_file
        return local_chr_files
    
    #download genome files        
    def download_dna(self):
        url = self.url['dna_fa']
        #get genome files
        #get html
        lines = web(url).get_html()
        chr_files = self.dna_files(lines)
        
        #download and decompress genome files
        local_chr_files = {}
        for key in chr_files.keys():
            self.ver = re.sub(r"\.chromosome.*", '', chr_files[key])
            gz_file = myIO.file_os(url+chr_files[key]).download(self.out_dir)
            #decompress file
            #ungz_file=myIO.file_os(gz_file).decompress_gz()
            local_chr_files[key] = gz_file
        #combine fa files
        out_file = self.out_dir+self.ver+'.fa'
        #print out_file
        myGenome.genome(out_file).combine_fa(local_chr_files)
        return local_chr_files, out_file
        
    #usually this method is the top choice        
    def batch_download(self,annot_type):
        #annot_type=['dna_fa','cds_fa','pep_fa','ncrna_fa','gff3', 'gtf']
        for a in annot_type:
            if a == 'dna_fa':
                files_turple = self.download_dna_fa()
                genome_fa = files_turple[1]
            elif a == 'gff3' or a == 'gtf':
                files_dict = self.download_annot_file(a)
                for f in files_dict:
                    self.match_fa_gtf(genome_fa, files_dict[f])
            else:
                self.download_chr_file(a)

################################################################################
#download the newest genome sequences from NCBI
class NCBI:
    def __init__(self, specie, out_dir):
        self.specie = specie
        self.out_dir = myIO.dir_os(out_dir+specie).create_dir()
        #initiate url list
        self.url_list()

#
    def url_list(self):
        species = {'human':'Homo_sapiens', 'mouse':'Mus_musculus',\
                 'rat':'Rattus_norvegicus', 'maize':'Zea_mays'}
        #annot_type=['dna_fa','cds_fa','pep_fa','ncrna_fa','gff3', 'gtf']
        self.url = {'pub':'ftp://ftp.ncbi.nlm.nih.gov/genomes/'}
        if self.specie in species.keys():
            head=self.url['pub'] + species[self.specie] + '/
            #genome DNA in fasta format
            self.url['dna_fa'] = head + 'Assembled_chromosomes/seq/'
            #annotations in gff3 format
            self.url['gff'] = head + 'GFF/'
            #protein amino acid in fasta format
            self.url['protein'] = head + 'protein/'
            #RNA
            self.url['RNA'] = head + 'RNA/'

#get DNA choromosome file names
    def dna_files(self, lines):
        chr_files = {}
        for line in lines: 
            if line.find('_ref_')>0 and line.find('.fa.')>0 and line.find('_chr')>0:
                items = line.split()
                file_name = items[-1]
                chromosome = re.findall(r"chr(\w*)", file_name)[0]
                chr_files[chromosome] = file_name 
                #print '=', file_name
        #print dict
        #myDict.basic(chr_files).print_dict()
        return chr_files
            
#download genome DNA
    def download_dna(self):
        #get html
        lines = web(self.url['dna_fa']).get_html()
        chr_files = self.dna_files(lines)
        
        #download and decompress genome files
        local_chr_files = {}
        for key in chr_files.keys():
            #release version
            self.ver = re.sub(r"_chr.*", '', chr_files[key])
            url = self.url['dna_fa']+chr_files[key]
            gz_file = myIO.file_os(url).download(self.out_dir)
            #decompress file
            #ungz_file=myIO.file_os(gz_file).decompress_gz()
            local_chr_files[key] = gz_file
        #combine fa files
        out_file = ''.join([self.out_dir, self.ver,'_dna.fa'])
        #print out_file
        myGenome.genome(out_file).combine_fa(local_chr_files)
        return local_chr_files, out_file

#download genome files        
    def download_annot(self,genome_type):
        #print genome_type
        url = self.url[genome_type]
        #get html and the list of files
        url_dir, url_files = web(url).ls_html()
        #print url_files
        #download and decompress genome files
        local_chr_files = {}
        for file_url in url_files.keys():
            gz_file = myIO.file_os(file_url).download(self.out_dir)
            #decompress file
            ungz_file = myIO.file_os(gz_file).decompress_gz()
            local_chr_files[file_url]=ungz_file
        return local_chr_files

#####################3
class uniprot:
        def __init__(self,out_dir):
            self.out_dir = myIO.dir_os(out_dir).create_dir()
            self.url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/'
            
        
        def download_idmapping(self):
            #get web file list
            url_idmapping = self.url+'knowledgebase/idmapping/by_organism/'
            web_dir, web_files = web(url_idmapping).ls_html()
            #print web_files
            
            #select file
            file_names = filter(lambda x: '.dat.' in x, web_files.values())
            file_names.sort()
            file_name = mySystem.system().select_key(file_names, 'Select web file')
            #download idmapping dat file
            url_file = url_idmapping + file_name
            local_file = self.out_dir + file_name
            web(url_file).download_file(local_file)
            #decompress file
            ungz_file = myIO.file_os(local_file).decompress_gz()
            print('Save ', url_file, ' as ', ungz_file)
            return ungz_file

            
#end



