'''
Created on 06/11/2013

@author: kibanez
'''

# /usr/bin/python

import shlex, os, shutil, subprocess, sys 
import inspect
from multiprocessing import Process, Manager

class mod_CNV:    
                      
    def __init__(self,list_files,ref_fasta,fasta_cnv_path,gatk_path,output_path,l_gender,annotation_file,logger):
        
        self.list_files = list_files       # list of GATK coverage files

        self.bed_analysis = None
        
        self.fasta_cnv_path = fasta_cnv_path
        
        self.gatk_path = gatk_path
        
        self.output_path = output_path
        
        self.l_gender = l_gender
        
        self.ref_fasta = ref_fasta
        
        self.annotation_file = annotation_file
        
        self.logger = logger
        
                
              
    def set_bed_analysis(self,bed_analysis):
        
        filename_bed_analysis = bed_analysis
            
        if os.path.exists(filename_bed_analysis):
            self.bed_analysis = bed_analysis       # not the bed file, its path
        else:
            raise IOError('mod_CopyNumberVariation.set_bed_analysis: The bed file of analysis has not been found: %s' % (filename_bed_analysis))

    # Creates an output folder (specified in cfg file) and the temporal scratch subfolder as well        
    def checkOutputFolder(self,outF):
        self.logger.info("Creating Output Folder...")
     
        if outF[len(outF)-1] == "/":
            outF = outF[:len(outF)-1]
    
        try:            
            os.mkdir(outF)
                    
        except:            
            self.logger.info("The folder '%s already exists'" %outF)
        
        individual_folder = os.path.join(self.output_path,"individual_results")
         
        try:            
            os.mkdir(individual_folder)
    
        except:
            if not os.path.isdir(outF):
                raise IOError('The folder %s does not exist and cannot be created. Maybe privileges problem? ' % (outF))
            
        return individual_folder         
       
    def perform_cnv(self):
       
        scratch_path = self.checkOutputFolder(self.output_path)        
  
        try:
            gc_content_template = scratch_path + "/gc_content_template.txt"
            # ACHTUNG !! el outBed ahora mismo no es el de las windows 150 ventanas (!!)
            args = ['java','-jar',self.gatk_path,'-T','GCContentByInterval','-R',self.fasta_cnv_path,'-L',self.bed_analysis,'-o',gc_content_template]   
            subprocess.Popen(args).wait()
         
        except:
            raise RuntimeError('mod_CNV: GCContentByInterval sucks')

        project_path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
        scriptPath   = os.path.join(project_path[:project_path.find('src')],'src/modules/scripts')
           
        self.logger.info("Normalizing test GATK coverage outputs ...")
        
        cases_files = []
        cases_files_rel = []  
        rScriptName = os.path.join(scriptPath,"normalized_DoC2.R")
        for index,i in enumerate(self.list_files):
            args = shlex.split("Rscript %s %s %s %s" %(rScriptName,i,gc_content_template,self.l_gender[index]))
            subprocess.call(args)
            aux = i + "_normalized"
            aux = aux + "_" + self.l_gender[index]
            cases_files.append(aux)
            aux2 = aux.split('/')
            aux2 = aux2[len(aux2)-1]           
            cases_files_rel.append(aux2)
        
        rScriptName = os.path.join(scriptPath,"cnv_estimation_individual2.R")
        samples = ";".join(cases_files)
        
        gender = ";".join(self.l_gender)
        args = shlex.split("Rscript %s %s %s %s %s %s" %(rScriptName, samples, self.output_path,self.bed_analysis,gender,self.annotation_file))
        subprocess.call(args)        