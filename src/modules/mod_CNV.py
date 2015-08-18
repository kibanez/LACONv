'''
Created on 06/11/2013

@author: kibanez
'''

# /usr/bin/python

import shlex, os, shutil, subprocess, sys 
import inspect
from multiprocessing import Process, Manager

class mod_CNV:    
    
    def __init__(self,list_files,ref_fasta,fasta_cnv_path,gatk_path,output_path,l_gender,panel_id,annotation_file,control_cnv,window=150,pipeline=1,all_controls = 0):
        
        self.list_files = list_files       # list of GATK coverage files

        if (window != ""):
            self.window = window        # length of window (default = 150)
        else:
            self.window = 150
            
        if (pipeline != ""):
            self.pipeline = pipeline    # pipeline = 0, when mod_CNV is called from tools/CNV_analysis --> it changes script_path parameter (!!)
        else:
            self.pipeline = 1     
        
        self.bed_analysis = None
        
        self.fasta_cnv_path = fasta_cnv_path
        
        self.gatk_path = gatk_path
        
        self.output_path = output_path
        
        self.l_gender = l_gender
        
        self.ref_fasta = ref_fasta
        
        self.panel_id = panel_id
        
        self.annotation_file = annotation_file
        
        self.control_cnv = control_cnv
        
        self.all_controls = all_controls
        
        #self.fastaIndex_cnv_path = fastaIndex_cnv_path
        
        #self.control_path = control_path
        
        #self.splitBam_path = splitBam
        
        #self.splitBam_template = splitBam_template
        
        
        
        
                
              
    def set_bed_analysis(self,bed_analysis):
        
        filename_bed_analysis = bed_analysis
            
        if os.path.exists(filename_bed_analysis):
            self.bed_analysis = bed_analysis       # not the bed file, its path
        else:
            raise IOError('mod_CopyNumberVariation.set_bed_analysis: The bed file of analysis has not been found: %s' % (filename_bed_analysis))

    # Creates an output folder (specified in cfg file) and the temporal scratch subfolder as well        
    def checkOutputFolder(self,outF):
        print "Creating Output Folder :",
     
        if outF[len(outF)-1] == "/":
            outF = outF[:len(outF)-1]
    
        try:
            os.mkdir(outF)        
        except:
            print "The folder '%s already exists'" %outF
        
        try:
            os.mkdir(outF + "/scratch")
    
        except:
            if not os.path.isdir(outF):
                raise IOError('The folder %s does not exist and cannot be created. Maybe privilegies problem? ' % (outF))        
    
        # control and test subfolders within scratch
        
        try:
            os.mkdir(outF + "/scratch/control")
            
        except: 
            if not os.path.isdir(outF + "/scratch"):
                raise IOError('The folder %s does not exist and CONTROL folder cannot be created. Maybe privilegies problem? ' % (outF + "/scratch"))

        try:
            os.mkdir(outF + "/scratch/test")
            
        except: 
            if not os.path.isdir(outF + "/scratch"):
                raise IOError('The folder %s does not exist and TEST folder cannot be created. Maybe privilegies problem? ' % (outF + "/scratch"))

        
        return outF
    
  
    def depth_of_coverage(self,gatk_path,fasta_file,bed,outFile,l_bams):
        try:
            args = ['java','-jar',gatk_path,'-T','DepthOfCoverage','-omitBaseOutput','-omitLocusTable','-R',fasta_file,'-I',l_bams,'-L',bed,'-o',outFile]   
            subprocess.Popen(args).wait()
        except:
            raise RuntimeError('mod_CNV.mod_CNV.depth_of_coverage: Error when launching gatk DepthOfCoverage \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))


        # java -jar /opt/gatk/GenomeAnalysisTK.jar -T DepthOfCoverage -omitLocusTable -R /mnt/bioinfo01-ingemm/ref/Ref/hg19/fasta/ucsc.hg19.fasta -I R63-4_S2_align.realign.recal.bam -L /mnt/bioinfo01-ingemm/scratch/endoscreen_viejo/jc_target_genes.bed -o R63-4

    def depth_of_coverage_individual(self,gatk_path,fasta_file,bed,bam_file,scratch_path):
        try:
            bam_name = bam_file.split('/')[len(bam_file.split('/'))]
            
            outFile = "coverage"
            args = ['java','-jar',gatk_path,'-T','DepthOfCoverage','-omitBaseOutput','-omitLocusTable','-R',fasta_file,'-I',bam_file,'-L',bed,'-o',outFile]   
            subprocess.Popen(args).wait()
        except:
            raise RuntimeError('mod_CNV.mod_CNV.depth_of_coverage: Error when launching gatk DepthOfCoverage \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
        
        
    def bam_template(self,bedF,bamF,fileName,output_path,type):
        
        args = shlex.split("samtools view -b -L %s %s" %(bedF,bamF))
        
        if (type == "control"):
        
            bam_template_output = os.path.join(output_path + '/scratch/control',fileName)
        
        elif (type == "test"):
        
            bam_template_output = os.path.join(output_path + '/scratch/test',fileName)

        iOutFile = open(bam_template_output,"w")
        output = subprocess.Popen(args,stdout=iOutFile).wait()
        iOutFile.close()
        
    def sam_to_bam(self,samfile, bamfile, index):
        args = shlex.split("samtools view -S -bt %s %s -o %s" %(index, samfile ,bamfile))
        try:
            subprocess.call(args)
        except:
            raise RuntimeError('mod_CNV.mod_CNV.sam_to_bam: Error when launching samtools view -S -bt %s %s -o %s \n %s\t%s' % (index, samfile ,bamfile, sys.exc_info()[0],sys.exc_info()[1]))
   
        
    def bam_to_sam(self,bamfile,samfile):
        args = shlex.split("samtools view -h -o %s %s" %(samfile,bamfile))

        try:
            subprocess.Popen(args).wait()
            return samfile
        except:
            raise RuntimeError('mod_CNV.mod_CNV.bam_to_sam: Error when launching samtools view -h -o %s %s \n %s\t%s' % (samfile ,bamfile, sys.exc_info()[0],sys.exc_info()[1]))
            
        try:    
            os.system("rm %s" %(bamfile))
        except:
            raise RuntimeError('mod_CNV.mod_CNV.bam_to_sam: Error when launching rm %s \n %s\t%s' % (bamfile, sys.exc_info()[0],sys.exc_info()[1]))            
    
    # removes chrM, ungl_chr, etc.
    def reheader_test(self,samfile,output_path,type,samName):
        try:
            iOutFile = open(samName,"w")
            os.system("grep -E -v '^@SQ[[:blank:]]SN[[:punct:]]chrM' %s | grep -E -v '^@SQ[[:blank:]]SN[[:punct:]]chrUn[[:punct:]]' |grep -E -v '^@SQ[[:blank:]]SN[[:punct:]]chr[[:digit:]][[:digit:]][[:punct:]]' |grep -E -v '^@SQ[[:blank:]]SN[[:punct:]]chr[[:digit:]][[:punct:]]' | grep -E -v '^@RG[[:blank:]]ID[[:punct:]]ERR'> %s" %(samfile,samName))
            iOutFile.close()
        except:
            raise RuntimeError('mod_CNV.mod_CNV.reheader_test: Error when launching  grep -E -v with %s and %s \n %s\t%s' % (samfile,samName,sys.exc_info()[0],sys.exc_info()[1]))

        if (type == "test"):
            try:
                os.system("rm %s" %(samfile))
            except:
                raise RuntimeError('mod_CNV.mod_CNV.reheader_test: Error when launching  rm %s \n %s\t%s' % (samfile,sys.exc_info()[0],sys.exc_info()[1]))                
       
        return samName

    def create_indexBam(self,bamF,folder):
        bam_file = os.path.join(folder,bamF)
        args = shlex.split("samtools index %s" %(bam_file))
        try:
            subprocess.call(args)
        except:
            raise RuntimeError('mod_CNV.mod_CNV.create_indexBam: Error when launching  samtools index %s \n %s\t%s' % (bam_file,sys.exc_info()[0],sys.exc_info()[1]))
    
    
    def removeTempFolder(self,tempFolderPath):    
        shutil.rmtree(tempFolderPath)
#        print "Temp Folder Removed"

    #def splitBam(bamFile,splitBam_jarPath):
    #    args = shlex.split("java -jar", comments, posix)
           
    def splitBam(self,bamF,bamName,ref_fasta,splitBam_path,template,output_path):
        '''java -jar jvarkit-read-only/dist/splitbam.jar -p /home/azureuser/dataset/DistrofiaEjemplo/split_NA12878/__CHROM__.bam -R /home/azureuser/dataset/ref/hg19/fasta/chr1toY.fasta 
        -G /home/azureuser/dataset/CNVapps/CNAnorm/chr_1_S1 /home/azureuser/dataset/DistrofiaEjemplo/NA12878_chr_CFTRtemplate.bam '''
        
        '''
         new splitbam release: pierre changed also the parameters calling and the fasta thing:
         java -jar /home/kibanez/Descargas/jvarkit-master/dist/splitbam.jar OFP=__CHROM__.bam R=/mnt/bioinfo01-ingemm/ref/Ref/hg19/fasta/chr1toY.fasta GP=/mnt/bioinfo01-ingemm/ref/DBs/Templates/template_splitBam I=../../../1-CR2/1-CR2_S1_align.realign.recal.bam'''
        
        # we create a folder with the name and inside it ==> __CHROM__.bam
        
        fileName, fileExtension = os.path.splitext(bamName)
        folder = "split_" + fileName
        #split_nameBam(without .bam)/__CHROM__.bam        
        outF = output_path + "/scratch/test/"        
        outF = outF + folder + "/__CHROM__.bam"
         
        
        #args = ['java','-jar',splitBam_path,'-p',outF,'-R',fasta_cnv_path,'-G',template,bamF]
        
                
    # new release
        args = ['java','-jar',splitBam_path,'OFP=',outF,'R=',ref_fasta,'GP=',template,'I=',bamF]
        
        try:
            output    = subprocess.Popen(args).wait()
        except:
            raise RuntimeError('mod_CNV.mod_CNV.splitBam: Error when launching  splitBam \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
              
        #iOutFile = open(outF, "w")
        #output    = subprocess.Popen(args, stdout = iOutFile).wait()
        #iOutFile.close()
        
        
    # args=(fasta_cnv_path,bam_control,i (el nombre),output_path
    # not used in the pipeline (only in the preparation of it)        
    def reduce_bams(self,fasta_file,input_Bam,output_Bam,output_path,gatk_path):
        # ACHTUNG !!! problema importate: GATK al ser un java, usa como temporary folder /tmp
        # En azure no tenemos casi espacio => hay qeu moverlo a /home/azureuser/dataset/tmp
        # export _JAVA_OPTIONS=-Djava.io.tmpdir=/home/azureuser/dataset/tmp
        
        outBam = output_Bam + "_reduced.bam"
        outputBam = os.path.join(output_path + '/scratch/control',outBam)

        args = ['java','-jar',gatk_path,'-T','ReduceReads','-R',fasta_file,'-I',input_Bam,'-o',outputBam]
        try:   
            subprocess.Popen(args).wait()
        except:
            raise RuntimeError('mod_CNV.mod_CNV.reduce_bams: Error when launching  gatk reduceReads \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
    
    
    def perform_cnv(self):
       
        self.checkOutputFolder(self.output_path)
        scratch_path = self.output_path + "/scratch"
  
        try:
            gc_content_template = scratch_path + "/gc_content_template.txt"
            # ACHTUNG !! el outBed ahora mismo no es el de las windows 150 ventanas (!!)
            args = ['java','-jar',self.gatk_path,'-T','GCContentByInterval','-R',self.fasta_cnv_path,'-L',self.bed_analysis,'-o',gc_content_template]   
            output = subprocess.Popen(args).wait()
         
        except:
            raise RuntimeError('mod_CNV: GCContentByInterval sucks')

        project_path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))
        scriptPath   = os.path.join(project_path[:project_path.find('src')],'src/modules/scripts')
           
        print "Normalizing test GATK coverage outputs..."
        
        # If there are samples with the same panel's bed => we are also going to use them as controls (more...better)       
        
        # Look if the controls_folder (panel_id) exists and take the elements inside it
        controls_folder = self.control_cnv + self.panel_id
        
         
        cases_files = []
        cases_files_rel = []  
        rScriptName = os.path.join(scriptPath,"normalized_DoC2.R")
        for index,i in enumerate(self.list_files):
            args = shlex.split("Rscript %s %s %s %s" %(rScriptName,i,gc_content_template,self.l_gender[index]))
            subprocess.call(args)
            aux = i + "_normalized"
            print "%s normalized" %(i)
            aux = aux + "_" + self.l_gender[index]
            cases_files.append(aux)
            aux2 = aux.split('/')
            aux2 = aux2[len(aux2)-1]           
            cases_files_rel.append(aux2)
            try:
                # Copy the new file to  /mnt/bioinfo01-ingemm/ref/DBs/Data/ + panel_id
                if (os.path.isfile(controls_folder + "/" + aux2 )):
                    print 'The file %s already exists. The system will not replace it. Check if it is ok'  %(aux2)
                else:
                    # ACHTUNG !! 01/04/2014: mirar el avg coverage total => dependiendo del panel, se filtran aquellas que no tienen un cov
                    
                    os.system("cp %s %s" %(aux,controls_folder))
            except:
                raise RuntimeError('mod_CNV.perform_cnv: Error when copying the normalized coverage file \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
        
        
        l_controls = []
        
        if (os.path.isdir(controls_folder)):
            l_controls = os.listdir(controls_folder) 
        else:
            print "It does not exist the folder corresponding to the panel_id. Check if there are not more controls for the panel_id introduced."
        
        # ACHTUNG!! only take those files that are different from the samples being currently processed
        l_controls2 = []
        for index,item in enumerate(l_controls):
            # ACTHUNG !! hay que quedarse con el nombre relativo , no absoluto !!!
            if item not in cases_files_rel:                
                l_controls2.append(item)
        # Add the absolute path
        l_controls = []
        if (len(l_controls2) > 0):
            for i in l_controls2:
                aux = os.path.join(controls_folder,i)
                l_controls.append(aux)
        
        rScriptName = os.path.join(scriptPath,"cnv_estimation_individual2.R")
        samples = ";".join(cases_files)
        if (len(l_controls) != 0):
            if self.all_controls == 1:
                controls = ";".join(l_controls)
            else:
                controls = 0
        else:
            controls = 0
        
        gender = ";".join(self.l_gender)
        args = shlex.split("Rscript %s %s %s %s %s %s %s" %(rScriptName, samples, self.output_path,self.bed_analysis,gender,controls,self.annotation_file))
        subprocess.call(args)        