'''
Created on 21/11/2013

@author: kibanez
'''

# /usr/bin/python

import shlex, os, subprocess, sys

#Absolute Path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))


class c_CoverageControl:
    
    def __init__(self,l_bam,logger,threshold=15,**kwargs):
        
        self.bam_list = l_bam       # list of bam files of a same run

        self.threshold    =  threshold
        
        self.bed_analysis = None
        
        self.logger = logger

                
              
    def set_bed_design(self,bed_design):
        
        self.bed_design   =  bed_design
        

    def set_bed_analysis(self,bed_analysis):
        
        filename_bed_analysis = bed_analysis
            
        if os.path.exists(filename_bed_analysis):
            
            self.bed_analysis = bed_analysis 
                  
        else:
            
            raise IOError('mod_CovControl.c_coverage.set_bed_analysis: The bed file of analysis has not been found: %s' % (filename_bed_analysis))
        

    def depth_of_coverage(self,gatk_path,fasta_file,bed,outFile,l_bams):
        
        try:
            args = ['java','-jar',gatk_path,'-T','DepthOfCoverage','-omitBaseOutput','-omitLocusTable','-R',fasta_file,'-I',l_bams,'-L',bed,'-o',outFile]   
            subprocess.Popen(args).wait()
        
        except:
            raise RuntimeError('mod_ClinicNV.mod_ClinicNV.depth_of_coverage: Error when launching gatk DepthOfCoverage \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
        
        
    def perform_coverage_control_gatk(self,path,gatk_path,fasta_file):
        l_cov = []
        for i in self.bam_list:
            try:
                bam_f = i.split('/')
                name_f = bam_f[len(bam_f) - 1]
                
                fileName = name_f +  "_GATKcoverage"                
                fileName = os.path.join(path,fileName)
                l_cov.append(fileName)
                args = ['java','-jar',gatk_path,'-T','DepthOfCoverage','-omitBaseOutput','-omitLocusTable','-R',fasta_file,'-I',i,'-L',self.bed_analysis,'-o',fileName]  
                self.logger.info("GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R %s -I %s -L %s -o %s" %(fasta_file,i,self.bed_analysis,fileName))
                subprocess.Popen(args).wait()
                
            except:
                raise RuntimeError('mod_CoverageControl.perform_coverage_control_gatk: Error when launching GATK \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
            
        sys.stdout.write("End of performing the coverage of each sample within the panel\n")
        return l_cov