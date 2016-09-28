'''
Created on 02/06/14

@author: kibanez
'''

#!/usr/bin/python

import sys, re, shlex , os, string, urllib, time, math, random, subprocess, shutil

from multiprocessing import Process, Manager

import inspect

from itertools import izip_longest,groupby

from operator import itemgetter

import ConfigParser

import optparse

from os import path as osp

localModulesBase = osp.dirname(osp.realpath(__file__))


modulesRelDirs = ["modules/"]

for moduleRelDir in modulesRelDirs:
    sys.path.insert(0,osp.join(localModulesBase,moduleRelDir))
        

from mod_CoverageControl import c_CoverageControl

import mod_CNV 
import logging

import numpy


######################################################################

class OptionParser(optparse.OptionParser):

    def check_required (self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)
        
        if atrib is None:
#            self.error("%s option not supplied" % option)
            return False
        else:
            return True
            

######################################################################

def read_cfg_file(cfg_filename):
    
    fi = open(cfg_filename,'r')
    
    config = ConfigParser.ConfigParser()
    config.readfp(fi)
    
    hash_cfg = {}
        
    for field in config.options('INPUT'):
        hash_cfg[field] = config.get('INPUT',field)
    
    for field in config.options('REFERENCE'):
        hash_cfg[field] = config.get('REFERENCE',field)
    
    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT',field)
     
    for field in config.options('BDs'):
        hash_cfg[field] = config.get('BDs',field)
        
    for field in config.options('SOFTWARE'):
        hash_cfg[field] = config.get('SOFTWARE',field)
        
    fi.close()
    
    return hash_cfg

#######################################################################

def parse_coverage(f_cov):
    
    #Sample_id    Gene    Refseq    Exon    Chr    Exon Start    Exon End    Exon length    Avg_Coverage    % Coverage < 5    % Coverage 5-10    % Coverage 10-15    % Coverage 15-20    % Coverage < 20
    
    hash_table = {}
    
    fi = open(f_cov,'r')
    
    l_lines = map(lambda l: l.strip().split('\t') , fi.readlines())[1:]
    map(lambda (i,(sample,gene,transc,exon,chr,start,end,len,avg_cov,cov_5,cov_5_10,cov_10_15,cov_15_20,cov_less_20)): hash_table.setdefault((i,gene,transc,exon,chr,start,end),[]).append([sample,float(avg_cov.replace(',','.')),float(cov_less_20.replace(',','.'))]), enumerate(l_lines))
    
    fi.close()
    
    return hash_table

#######################################################################
# cleans the bed file, removing lines where different things apart from "chr start end" appear

def clean_bed(bed_file,alignment_path):
    bed_name = bed_file.split('/')
    bed_name = bed_name[len(bed_name) - 1]
    
    bedName, bedExtension = os.path.splitext(bed_name)
    
    new_analysis_bed = bedName + "_reheadered" + bedExtension
    
    new_analysis_bed_path = alignment_path + "/" + new_analysis_bed
    aux = alignment_path + "/" + "aux.bed"
    aux2 = alignment_path + "/" + "aux2.bed"  
    
    iOutFile = open(new_analysis_bed_path,"w")
    iOutFile2 = open(aux,"w")
    # we conserve all those lines with the following info: chr \t start \t end
    os.system("grep -e '^chr[[:upper:]]' -e '^chr[[:digit:]][[:blank:]]' -e '^chr[[:digit:]][[:digit:]][[:blank:]]' %s > %s" %(bed_file,aux))
    # we remove "chrM" lines if they appear 
    os.system("grep -v '^chrM' %s > %s" %(aux,aux2))
    os.system("sort -k1,1V %s > %s" %(aux2,new_analysis_bed_path))
    iOutFile.close()
    iOutFile2.close()

    os.remove(aux)
    os.remove(aux2)
    return new_analysis_bed_path


#######################################################################

def __create_reference_tables(refseq_filename,output_path):
    
    ref_annotation_bed  = os.path.join(output_path,os.path.basename(refseq_filename) + str(random.randrange(1,500)) + '.ref.bed')
    ref_annotation_file = os.path.join(output_path,os.path.basename(refseq_filename)+ str(random.randrange(1,500)) + '.ref.annot')
    
    """
        Columns in the table:
            hg19.refGene.name    
            hg19.refGene.chrom    
            hg19.refGene.strand    
            hg19.refGene.cdsStart    
            hg19.refGene.cdsEnd    
            hg19.refGene.exonCount    
            hg19.refGene.exonStarts    
            hg19.refGene.exonEnds    
            hg19.refGene.name2    
            hg19.refSeqStatus.status
            
        Columns in ref_annotation_file:
            Chr    
            txStart    
            txEnd    
            exonCount    
            exonStart    
            exonEnd    
            Gene    
            Refseq
    """
    
    f = open(refseq_filename,'r')
    l_lines = map(lambda x: x.strip().split('\t'), f.readlines())[1:]
    f.close()
    
    #NM_032291    chr1    +    67000041    67208778    25    66999824,67091529,67098752,67101626,67105459,67108492,67109226,67126195,67133212,67136677,67137626,67138963,67142686,67145360,67147551,67154830,67155872,67161116,67184976,67194946,67199430,67205017,67206340,67206954,67208755    67000051,67091593,67098777,67101698,67105516,67108547,67109402,67126207,67133224,67136702,67137678,67139049,67142779,67145435,67148052,67154958,67155999,67161176,67185088,67195102,67199563,67205220,67206405,67207119,67210768    SGIP1    Validated
    
    #select the principal isoform as the largest transcript
    hash_gene_2_isoform = {} 
    hash_table          = {}
    hash_index_gene     = {}
    
    for index,line in enumerate(l_lines):
        
        refseq_id = line[0]
        chrm      = line[1]
        strand    = line[2]
        cds_start = line[3]
        cds_end   = line[4]
        n_exon    = line[5]
        l_exon_s  = map(int,line[6].split(','))
        l_exon_e  = map(int,line[7].split(','))
        gene_name = line[8]
        
        if len(l_exon_s) <> len(l_exon_e):
            raise RuntimeError('coverage_statistics.__create_reference_tables: The number of initial exon coordinates does not agree with the number of end coordinates: %s\n' % (gene_name))
        if len(l_exon_s) <> int(n_exon):
            raise RuntimeError('coverage_statistics.__create_reference_tables: The number of exons does not agree with exon_count: %s: %s\n' % (gene_name,n_exon))
        
        cs = int(cds_start)
        ce = int(cds_end)
        
        if strand == '+':
            l_coding_exons = filter(lambda (i,s,e): float(min(e,ce)-max(s,cs))/(e-s)>0, zip(range(1,len(l_exon_s)+1),l_exon_s,l_exon_e))
        else:
            l_coding_exons = filter(lambda (i,s,e): float(min(e,ce)-max(s,cs))/(e-s)>0, zip(range(len(l_exon_s)+1,1,-1),l_exon_s,l_exon_e))
        
        l_exon_coord = map(lambda (i,s,e): (i,max(s,cs),min(e,ce)), l_coding_exons)
        
        hash_gene_2_isoform.setdefault(gene_name,[]).append((sum(map(lambda (i,s,e): (e-s) , l_exon_coord)),refseq_id))
        hash_table[refseq_id] = (chrm,l_exon_coord,gene_name)
        hash_index_gene.setdefault(gene_name,[]).append(index)
        
    l_genes = map(itemgetter(1),sorted(map(lambda gene: (sorted(hash_index_gene[gene])[0],gene), hash_index_gene.keys())))
    
    f_bed   = open(ref_annotation_bed,'w')
    f_annot = open(ref_annotation_file,'w')
    
    f_annot.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("Chr","txStart","txEnd","exonCount","exonStart","exonEnd","Gene","Refseq"))
    for gene_name in l_genes:
        
        length_trans,refseq_id = sorted(hash_gene_2_isoform[gene_name],reverse=True)[0]
        
        if length_trans == 0:
            continue
        
        (chrm,l_exon_coord,gene_name) = hash_table[refseq_id]
    
        f_bed.write('\n'.join(map(lambda (i,s,e): "%s\t%d\t%d\t%s\t%s" % (chrm,s,e,gene_name,refseq_id), l_exon_coord))+'\n')
        f_annot.write('\n'.join(map(lambda (i,s,e): "%s\t-\t-\t%d\t%d\t%d\t%s\t%s" % (chrm,i,s,e,gene_name,refseq_id), l_exon_coord))+'\n') 
    
    f_bed.close()
    f_annot.close()    

    
    return ref_annotation_bed,ref_annotation_file


#######################################################################

def run(argv=None):
    
    if argv is None: argv = sys.argv    
   
    parser = OptionParser(add_help_option=True,description="The script performs CNV estimation within the regions of interest")
    
    parser.add_option("--cfg",default=None,help="Config file with the complete information of the target regions and paths of the files needed for the calling",dest="f_cfg")

                    
    # Se leen las opciones aportadas por el usuario
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)
    
    if not parser.check_required("--cfg"):
        raise IOError('The cfg file does not exist')
        
               
    try:
        
        if options.f_cfg <> None:
            
            cfg_file = options.f_cfg        
          
            if not os.path.exists(cfg_file):
                raise IOError('The file %s does not exist' % (cfg_file))
            
            hash_cfg = read_cfg_file(cfg_file)
           

            cnv_output_path = hash_cfg.get('cnv_path','')
            controlCoverage_path = hash_cfg.get('coverage_path','')
            alignment_path  = hash_cfg.get('alignment_path','')
            ref_fasta = hash_cfg.get('ref_fasta','')
            fasta_cnv_path = hash_cfg.get('ref_fasta_cnv','')
            gatk_path = hash_cfg.get('gatk_path','')
            analysis_bed = hash_cfg.get('analysis_bed','')
            l_samples  = hash_cfg.get("sample_names",'').split(',')            
            l_gender = hash_cfg.get("sample_gender",'').split(',')
            window_length = hash_cfg.get("window_length",'')
            annotation_file = hash_cfg.get("annotation_file",'')
            
            
            if not os.path.exists(alignment_path):
                raise IOError('The path does not exist. %s' % (alignment_path))
           
            if not os.path.isfile(ref_fasta):
                raise IOError('The file does not exist. %s' % (ref_fasta))
            
            if not os.path.isfile(fasta_cnv_path):
                raise IOError('The file does not exist. %s' % (fasta_cnv_path))
            
            if not os.path.exists(controlCoverage_path):
                os.mkdir(controlCoverage_path)
            
            if not os.path.exists(cnv_output_path):
                os.mkdir(cnv_output_path)

            if not os.path.isfile(gatk_path):
                raise IOError('The file does not exist. %s' % (gatk_path))
            
            if not os.path.isfile(analysis_bed):
                raise IOError('The file does not exist. %s' % (analysis_bed))

            
            if not os.path.isfile(annotation_file):
                raise IOError("annotation_file not exist. %s" % (annotation_file))


                
            #Configure logger
            formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(logging.INFO)
            logger = logging.getLogger("preprocess")
            logger.setLevel(logging.INFO)
            logger.addHandler(console)
            
            l_bams = []
            for bam_f in l_samples:
                abs_path = os.path.join(alignment_path,bam_f)
                if not os.path.exists(abs_path):
                    raise IOError("The bam file does not exist. Check if the introduced path is correct: %s" %(abs_path))
                else:
                    l_bams.append(abs_path)
                
            logger.info("CNV estimation will be done in the following files: %s \n" %(l_bams))
            
            
            logger.info("Human genome bed generation, for the annotation process")
            refseq_filename = "/ingemm/ref/DBs/genes_tables/RefSeq_table.txt"
            (ref_annotation_bed,ref_annotation_file) = __create_reference_tables(refseq_filename,controlCoverage_path)
            
            
    
            logger.info("Bed generation")
            new_analysis_bed_path = clean_bed(analysis_bed,alignment_path)
            
            # before calling for CNV, it is necessary to create GATK coverage files , normalize them and do the estimation CNV 
            
            # 1 - GATK coverage average calling            
            cov_control = c_CoverageControl(l_bams,logger)
            cov_control.set_bed_analysis(new_analysis_bed_path)            
            
            # all the resulting files are generated in coverage/files directory
            path_files = os.path.join(controlCoverage_path,"files") 
            if not os.path.exists(path_files):
                    os.mkdir(path_files)            
            
            logger.info("Starting Coverage Control...")
            cov_output = cov_control.perform_coverage_control_gatk(path_files,gatk_path,ref_fasta)
                    
            # 2 - Normalization and 3 - CNV estimation calling (via mod_CNV)
            
            l_covFiles = []
            for i in cov_output:
                aux = i + ".sample_interval_summary"
                l_covFiles.append(aux)
    
            
            aux1 = filter(lambda x:'H' in x, l_gender)
            aux2 = filter(lambda x:'M' in x, l_gender)
            aux3 = filter(lambda x:'X' in x, l_gender)
             
            if (len(aux1) + len(aux2) + len(aux3) == len(l_gender)):
                
                logger.info("CNV estimation starts")
                cnv_estimation = mod_CNV.mod_CNV(l_covFiles,ref_fasta,fasta_cnv_path,gatk_path,cnv_output_path,l_gender,annotation_file,logger)
                cnv_estimation.set_bed_analysis(new_analysis_bed_path)
                cnv_estimation.perform_cnv()
             
            else:
                 
                raise IOError('The gender list must have only M (xx) or H (xy) or X (unknown) characters. Please review the cfg file.')
            
            
            logger.info("CNV estimation done! ")    
        
    except:
        print >> sys.stderr , '\n%s\t%s' % (sys.exc_info()[0],sys.exc_info()[1])
        sys.exit(2)

############################################################################333

if __name__=='__main__':
    
    run()


