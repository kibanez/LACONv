'''
Created on 21/11/2013

@author: kibanez
'''

# /usr/bin/python

import shlex, os, shutil, subprocess, sys
import shutil
import pandas as pd
import numpy as np
from multiprocessing import Process, Manager
from pybedtools import BedTool
import warnings

#Absolute Path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))


class c_CoverageControl:
    
    def __init__(self,l_bam,l_samples,l_ids,threshold=15,**kwargs):
        
        self.bam_list = l_bam       # list of bam files of a same run (same bed, same panel)
        
        self.l_samples= l_samples
        
        self.l_ids = l_ids

        self.threshold    =  threshold
        
        self.bed_analysis = None
        #self.bed_design   =  kwargs.get('design_bed','')
        #self.bed_analysis =  kwargs.get('analysis_bed','')   # to be calculated from the bam file

                
              
    def set_bed_design(self,bed_design):
        
        self.bed_design   =  bed_design
        

    def set_bed_analysis(self,bed_analysis):
        
        filename_bed_analysis = bed_analysis
            
        if os.path.exists(filename_bed_analysis):
            #self.feat_analysis = BedTool(filename_bed_analysis)
            self.bed_analysis = bed_analysis       # not the bed file, its path
        else:
            raise IOError('mod_CovControl.c_coverage.set_bed_analysis: The bed file of analysis has not been found: %s' % (filename_bed_analysis))
        
#        hash_feat_analysis_tmp  = {}
#        self.hash_feat_analysis = {}
        
#        map(lambda f: hash_feat_analysis_tmp.setdefault(f.chrom,[]).append(f) , self.feat_analysis)
        
#        for chrm in hash_feat_analysis_tmp.keys():
#            self.hash_feat_analysis[chrm] = BedTool(hash_feat_analysis_tmp[chrm])


    # Creates an output folder (specified in cfg file) and the temporal scratch subfolder as well        
    def checkOutputFolder(self,outF):
        print "Creating Output Folder :",
     
        if outF[len(outF)-1] == "/":
            outF = outF[:len(outF)-1]
    
        try:
            os.mkdir(outF)        
        except:
            print "cannot create folder '%s'" %outF
            print "because it already exists"
            sys.exit(1)
        
        try:
            os.mkdir(outF + "/scratch")
    
        except:
            print "[ERROR: CANNOT CREATE THE scratch SUBFOLDER]"
            sys.exit(1)        
    
        return outF


    def depth_of_coverage(self,gatk_path,fasta_file,bed,outFile,l_bams):
        try:
            args = ['java','-jar',gatk_path,'-T','DepthOfCoverage','-omitBaseOutput','-omitLocusTable','-R',fasta_file,'-I',l_bams,'-L',bed,'-o',outFile]   
            subprocess.Popen(args).wait()
        except:
            raise RuntimeError('mod_ClinicNV.mod_ClinicNV.depth_of_coverage: Error when launching gatk DepthOfCoverage \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
        
        
    def coverage_parallel(self,bamFile,outFile):
        args = shlex.split("bedtools genomecov -ibam %s -bg " %(bamFile))
        
        print "bedtools genomecov -ibam %s -bg" %(bamFile)
        
        outCov = open(outFile,"w")
        subprocess.Popen(args, stdout = outCov).wait() 
        outCov.close()
    
    
    # Oct 2014: We add BAM QUALITY PARAMETERS in order to compute # reads on (1) ORIGINAL bam (after sam creation) (2) Dedupped bam (3) >Q20 bam and (4) %ON_TARGET bam. 
    # #total reads, #mapped reads and #unmapped reads and %     
    def bam_quality_parameters(self,alignment_path,sample,id,bed_file,samtools_path,composed_path=True):
        
        qual_info = []
        
        raw_bam_name = id + "_" + sample
        
        if composed_path:
            bam_path = os.path.join(alignment_path,id)
        else:
            bam_path = alignment_path
        
        # ORIGINAL bam
        orig_bam = os.path.join(bam_path,raw_bam_name + "_align_original.bam")
        
        if not os.path.exists(orig_bam):
            #raise IOError("The original bam does not exist. %s" % (orig_bam))
            warnings.warn("mod_CoverageControl.bam_quality_parameters: The original BAM file does not exist")
            
            orig_bam = os.path.join(bam_path,raw_bam_name+"_align.realign.recal.bam")
                    
        # number of reads
        args = [samtools_path,'view','-c',orig_bam]
        
        try:
            reads_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)

        
        qual_info.append(reads_orig.rstrip())
        
        # number of mapped reads
        args = [samtools_path,'view','-c','-F 4',orig_bam]
        
        try:
            mapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        
        qual_info.append(mapped_orig.rstrip())
        
        
        # number of unmapped reads
        args = [samtools_path,'view','-c','-f 4',orig_bam]
        
        try:
            unmapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        
        qual_info.append(unmapped_orig.rstrip())
            
        # DEDUPPED bam
        dedupped_bam = os.path.join(bam_path,raw_bam_name + "_align.realign.recal.bam")
        if not os.path.exists(dedupped_bam):
            raise IOError("The dedupped bam (delete duplicated reads) does not exist. %s" % (dedupped_bam))
        
        # number of reads
        args = [samtools_path,'view','-c',dedupped_bam]
        
        try:
            reads_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)

        
        qual_info.append(reads_orig.rstrip())
        
        # number of mapped reads
        args = [samtools_path,'view','-c','-F 4',dedupped_bam]
        
        try:
            mapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        
        qual_info.append(mapped_orig.rstrip())
        
        # number of unmapped reads
        args = [samtools_path,'view','-c','-f 4',dedupped_bam]
        
        try:
            unmapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        
        qual_info.append(unmapped_orig.rstrip())
        
        # >Q20 bam
        q20_bam = os.path.join(bam_path,raw_bam_name + "_align.realign.recal_Q20.bam")
        if not os.path.exists(q20_bam):
            raise IOError("The bam containing reads with mapping Q>20 does not exist. %s" % (q20_bam))
        
        # number of reads
        args = [samtools_path,'view','-c',q20_bam]
        
        try:
            reads_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)

        
        qual_info.append(reads_orig.rstrip())
        
        # number of mapped reads
        args = [samtools_path,'view','-c','-F 4',q20_bam]
        
        try:
            mapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        
        qual_info.append(mapped_orig.rstrip())
        
        
        # number of unmapped reads
        args = [samtools_path,'view','-c','-f 4',q20_bam]
        
        try:
            unmapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        
        qual_info.append(unmapped_orig.rstrip())
        
        
        # ON_TARGET bam 
        onTarget_bam = os.path.join(bam_path,raw_bam_name + "_align.realign.recal_onTarget.bam")
        tmpBam = os.path.join("/scratch", raw_bam_name + "_align.realign.recal_onTarget.bam")        
        # creamos el bam que contiene los reads WITHIN el bed que usamos en el analisis
        # intersectBed -abam 2-8908_S2_align_sorted_fixed_sorted_dedup_header.bam -b DE_Roche_v3.bed > intersect.bam
        
        try:
            '''
               intersectBed -abam any.bam -b any.bed \| samtools view -h - \| samtools view -bS - > tmp.bam"               
            '''
            #args = shlex.split("intersectBed -wa -abam %s -b %s" %(dedupped_bam,bed_file))    
            #args = shlex.split("intersectBed -wa -abam %s -b %s | samtools view -h - | samtools view -bS - > %s" %(dedupped_bam,bed_file,onTarget_bam))
            
            #out_bed = open(onTarget_bam,"w")
            #subprocess.Popen(args,stdout = out_bed).wait()            
            #out_bed.close()
            
            BedTool(dedupped_bam).intersect(bed_file,wa=True,output=tmpBam)
            shutil.copyfile(tmpBam, onTarget_bam)
            
            if os.path.exists(tmpBam):
                os.remove(tmpBam)
            
        except:         
            raise RuntimeError('Error when launching bedtools \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))

        if not os.path.exists(onTarget_bam):
            raise IOError("The bam containing reads within the targets (bed file) does not exist. %s" % (onTarget_bam))
        
        # number of reads
        args = [samtools_path,'view','-c',onTarget_bam]
        
        try:
            reads_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)

        
        qual_info.append(reads_orig.rstrip())
        
        # number of mapped reads
        args = [samtools_path,'view','-c','-F 4',onTarget_bam]
        
        try:
            mapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        qual_info.append(mapped_orig.rstrip())
        
        # number of unmapped reads
        args = [samtools_path,'view','-c','-f 4',onTarget_bam]
        
        try:
            unmapped_orig = subprocess.check_output(args)
        except subprocess.CalledProcessError, e:
            msg = 'mod_CoverageControl.bam_quality_parameters: Error in samtools view -c  :%s' % str(e)
            print msg
            raise RuntimeError(msg)
        
        
        qual_info.append(unmapped_orig.rstrip())
        return qual_info        
    
  # bedtools genomecov -ibam *.bam -bg  
  # NOTE: when using higher performance CPU: paralelyze !
    def perform_coverage_control_bed(self,coverage_path,samtools_path,alignment_path,l_sample_ids,composed_path=True):

        l_proc = []
        l_cov  = []
        
        for index,i in enumerate(self.bam_list):
            
            fileName, fileExtension = os.path.splitext(os.path.basename(i))
            bam_Q20 = fileName + "_Q20.bam"
            if composed_path:
                bam_Q20 = os.path.join(alignment_path + "/" + l_sample_ids[index],bam_Q20)
            else:
                bam_Q20 = os.path.join(alignment_path,bam_Q20)
                
            fileName = fileName + "_coverage_Q20"
            fileName = os.path.join(coverage_path,fileName)
            
            l_cov.append(fileName)          
            
            ################################################################################################################################
            # Oct 2014: Improvement within coverage analysis: After reads are dedupped (delete duplicated reads) we only count reads which have >Q20 (mapping quality)
            #################################################################################################################################                
            out_bam = open(bam_Q20,"w")
            samtools_sal    = subprocess.Popen([samtools_path,'view','-b','-q','20',i],stdin=subprocess.PIPE, stdout=out_bam, stderr=subprocess.PIPE,close_fds=True,bufsize=-1)
            (trash,logdata) = samtools_sal.communicate()
            out_bam.close()
                    
            if logdata <> "":
                if logdata.lower().find("error") <> -1:
                    raise RuntimeError("c_CoverageControl.perform_coverage_control_bed: Error in samtools view: %s" % (logdata))
                
            #################################################################################################################################
                
            
            #args = shlex.split("bedtools genomecov -ibam %s -bga" %(i))
            args = shlex.split("bedtools genomecov -ibam %s -bga" %(bam_Q20))   # ahora le metemos el bam filtrado con los reads con >Q20 (y el bam ya esta previamente dedupped)
            
            sys.stdout.write("bedtools genomecov -ibam %s -bga\n" %(bam_Q20))
            sys.stdout.flush()
                
                
            outCov = open(fileName,"w")
            bedtools_sal    = subprocess.Popen(args, stdin=subprocess.PIPE,stdout=outCov,stderr=subprocess.PIPE,close_fds=True,bufsize=-1)
            (trash,logdata) = bedtools_sal.communicate()   
            outCov.close()
            
            if logdata <> "":
                if logdata.lower().find("error") <> -1:
                    raise RuntimeError("c_CoverageControl.perform_coverage_control_bed: Error in bedtools genomecov: %s" % (logdata))
                

        sys.stdout.write("End of performing the coverage of each sample within the panel\n")
        sys.stdout.flush()
            
        return l_cov
    
    
    def perform_coverage_control_gatk(self,path,gatk_path,fasta_file):
        l_cov = []
        for index,i in enumerate(self.bam_list):
            try:
                bam_f = i.split('/')
                name_f = bam_f[len(bam_f) - 1]
                #fileName, fileExtension = os.path.splitext(name_f)
                #fileName = fileName + "_GATKcoverage"
                
                fileName = self.l_samples[index] + "_" + self.l_ids[index] +  "_GATKcoverage"                
                fileName = os.path.join(path,fileName)
                l_cov.append(fileName)
                args = ['java','-jar',gatk_path,'-T','DepthOfCoverage','-omitBaseOutput','-omitLocusTable','-R',fasta_file,'-I',i,'-L',self.bed_analysis,'-o',fileName]
                print "GATK -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R %s -I %s -L %s -o %s" %(fasta_file,i,self.bed_analysis,fileName)   
                subprocess.Popen(args).wait()
                
            except:
                raise RuntimeError('mod_CoverageControl.perform_coverage_control_gatk: Error when launching GATK \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))
            
        sys.stdout.write("End of performing the coverage of each sample within the panel\n")
        return l_cov

    
    def coverage_statistics(self,l_cov,analysis_bed,alignment_path):
        bed = pd.read_csv(analysis_bed,sep="\t",header=None)
        bed_chr = bed[0]            
        bed_start = bed[1]
        bed_end = bed[2]
        print bed.head()
        if (len(bed.columns) > 3):
            print len(bed.columns)
            bed_gene = bed[4]            
        else:
            bed_gene = bed[0]   # ChrX:start-end
            
             
             
        df_final = pd.DataFrame()
        
        print l_cov
        
        for i,file in enumerate(l_cov):
            if not os.path.isfile(file):
                raise IOError('The file does not exist. %s' %(file))
        
            data = pd.read_csv(file,sep="\t",header=None)
            print "Read the fucking big coverage file"
            # index <- which(((analyzed_bed$V2[j] <= data$V2) & (analyzed_bed$V3[j] >= data$V3) & (analyzed_bed$V1[j] == data$V1)) | 
            #((analyzed_bed$V2[j] >= data$V2) & (analyzed_bed$V3[j] <= data$V3) & (analyzed_bed$V1[j] == data$V1)) | 
            # ((analyzed_bed$V2[j] > data$V2) & (analyzed_bed$V2[j] <= data$V3) & (analyzed_bed$V1[j] == data$V1)) | 
            # ((analyzed_bed$V3[j] < data$V3) & (analyzed_bed$V3[j] >= data$V2) & (analyzed_bed$V1[j] == data$V1)))
            data_chr = data[0]
            data_start = data[1]
            data_end = data[2]
            data_cov = data[3]
            
            cov_avg = []
            for j in xrange(0,len(bed_chr)):
                index1 = set(data_start[data_start >= bed_start[j]].index) & set(data_end[data_end <= bed_end[j]].index) & set(data_chr[data_chr == bed_chr[0]].index)
                index2 = set(data_start[data_start <= bed_start[j]].index) & set(data_end[data_end >= bed_end[j]].index) & set(data_chr[data_chr == bed_chr[0]].index)
                index3 = set(data_start[data_start < bed_start[j]].index) & set(data_end[data_end >= bed_start[j]].index) & set(data_chr[data_chr == bed_chr[0]].index)
                index4 = set(data_end[data_end > bed_end[j]].index) & set(data_start[data_start <= bed_end[j]].index) & set(data_chr[data_chr == bed_chr[0]].index)
                index = set(index1) | set(index2) | set(index3) | set(index4) 
                if len(index) == 0:
                    mediana = -1
                else:
                    coverages = data_cov[index]
                    mediana =  np.median(coverages)
                print j
#                print mediana    
                cov_avg.append(mediana)
            
            if (i == 0):
                df = pd.concat([bed_chr,bed_start,bed_end,bed_gene],axis=1)    # Chr, Start, End, Gene (In exome analysis, there is more info.)
                df_averages = pd.DataFrame(cov_avg)
            else:
                df_averages = pd.concat([df_averages,pd.DataFrame(cov_avg)],axis=1)
                
            
            print "Control coverage done ...%s" %(file)
            
               
        l_names = []
        for i in l_cov:
            aux = i.split('/')
            l_names.append(aux[len(aux)-1])
         
            
        
        # Once all of them are computed => we compute the overall average of the samples
        all_cov = []
        for i in xrange(0,len(df_averages[0])):
            aux = np.median(df_averages.ix[i])
            all_cov.append(aux)
        
        cov_index = "all_coverage"
        l_names.insert(0, cov_index)
        #df_averages.columns = l_names
        
        # Add the overall coverage to the global dataframe
        df_final = pd.concat([df,pd.DataFrame(all_cov),df_averages],axis=1)
        df_final.columns = ['Chr','Start','End','Gene'] + l_names
               
        no_capture = (len(df_final[df_final[cov_index] == -1].index) / len(df_final[cov_index])) * 100
        less15 = (len(df_final[df_final[cov_index] < 15].index) / len(df_final[cov_index])) * 100
        # Write output
        #pd.save_xls(df_final,alignment_path)
        
        
        new_path = os.path.join(alignment_path,"CoverageControlPYTHON.xls")
        #df_final.to_excel(new_path,sep="\t")
        df_final.to_excel(new_path)
        print "Coverage Control generated"
        print "Percentage of no captured regions: %s" %(str(no_capture))
        print "Percentage of regions with coverage less than 15: %s" %(str(less15))
        
        no_capture = df_final[df_final[cov_index] == -1].index
        df_no_capture = pd.concat([pd.DataFrame(df_final[no_capture])])
        
        new_path2 = os.path.join(alignment_path,"CoverageControl_NoCapturedPYTHON.xls")
        df_no_capture.to_excel(new_path2)
                
        less15 = df_final[df_final[cov_index] < 15].index
        df_less15 = pd.concat([pd.DataFrame(df_final[less15])])
        new_path3 = os.path.join(alignment_path,"CoverageControl_less15PYTHON.xls")
        df_less15.to_excel(new_path3)        
        
        
        
    def perform_coverage_control(self,gatk_path,ref_fasta):
        # save in test.list the list of the bam files to process with GATK DoC
        
        test_list = os.getcwd() + "/test.list"
        f = open(test_list,"w")
        for index,item in enumerate(self.bam_list):
                f.write(item + "\n")
        f.close()
       
                
        l_cov = []     
        covTest = os.getcwd() + "/coverage_control"

#        depth_of_coverage(gatk_path,ref_fasta,self.feat_analysis,covTest,self.bam_list)



        try:
            args = ['java','-jar',gatk_path,'-T','DepthOfCoverage','-omitBaseOutput','-omitLocusTable','-R',ref_fasta,'-I',test_list,'-L',self.feat_analysis,'-o',str(covTest)]   
            subprocess.Popen(args).wait()
        except:
            raise RuntimeError('mod_CoverageControl.perform_coverage_control: Error when launching gatk DepthOfCoverage \n %s\t%s' % (sys.exc_info()[0],sys.exc_info()[1]))

        return covTest
