ReadBig <- function( File, sep="\t", header=TRUE, comment.char = "", skip = 0){

   tab5rows <- read.table( File, header = header, nrows = 5, stringsAsFactors=FALSE, skip=skip, sep=sep)
   classes <- sapply( tab5rows, class )
   Nrows <- system( paste("wc -l", File, sep=" "), intern=TRUE )
   Nrows <- as.numeric( unlist( strsplit( Nrows, split=" " ) )[ 1 ] )
   read.table( File, header = header,
   			  sep = sep,
              colClasses = classes,
              nrows = Nrows,
              comment.char = comment.char,
              skip = skip,
              , stringsAsFactors=FALSE)
}


# we separate data between XX and XY rows, in order to compare them
# SOLO AQUELLOS QUE COMPARTEN XX y XY (el resto, al comparar específicamente los sexuales)
remove_sexual <- function(df.sexual,avg_values,gender_list){
 xy = which(gender_list == "H")
 xx = which(gender_list == "M")
 
 # mirar cuantos hay de xy
 male = avg_values[xy]
 female = avg_values[xx]

 for (i in 1:2){
	if (i == 1) {
		sexo = male
	}else{
		sexo = female
	}
	m_aux2 = as.matrix(df.sexual[,sexo])
	m_aux2 = apply(m_aux2,1,as.numeric)
	m_aux2 = t(m_aux2)
	median_avg_row2 <- apply(m_aux2,1,median)
	if (i == 1){
		less_30xy = which(median_avg_row2 < 30)
	}else{
		less_30xx = which(median_avg_row2 < 30)
	}
 }

 # the best thing would be to separate those targets <30 in xx and those <30 in xy
 identicos = intersect(less_30xy,less_30xx)
 return(identicos)
}

 


sexual_CNV <- function (df,gender_list,sample,sample_name){
 xy = which(gender_list == "H")
 xx = which(gender_list == "M")

 test_index_norm_avg <- seq(from = 7, to = length(colnames(df)), by = 4) 

 if (sample %in% xx){
	group = xx
 }else{
	group = xy
 }

# ACHTUNG !!! ANTES DE NADA, BORRAMOS AQUELLAS ROWS QUE TENGAN COVERAGE < 30!!!
 average_indexes = seq(from = 6, to = length(colnames(df)), by = 4)
 m_aux = as.matrix(df[,average_indexes[group]])
 m_aux = apply(m_aux,1,as.numeric)
 m_aux = t(m_aux)
# median_avg_row <- apply(m_aux,1,median)
# less_30 = which(median_avg_row < 30)
# if (length(less_30) > 0){
#	 df = df[-less_30,]
# }
# aux = paste("SexualChromosomes",sample_name,sep="_")
# aux = paste(aux,"AverageCoverage_less30.tsv",sep="_")
# output_fileName = paste(output_path,aux,sep="/")
# tryCatch(write.table(file=output_fileName,out_less_30,quote=F,row.names=F,col.names=T,sep="\t"),error = function(e) print(paste(e$message,output_fileName,sep=" ")))



 if (length(group) > 1){  # at least there are 2 samples to compare...not enough but at least
	 c_values <- c() 
	 
	 for (i in 1:length(group)){
		# mean_reads contains the mean reads of each sample: sum(all reads in all the intervals) / num.intervals
		index = test_index_norm_avg[group]
		mean_reads <- sum(as.numeric(as.character(df[index][,i]))) / length(df$Target)
		c_sample = as.numeric(as.character(df[index][,i])) / mean_reads
		c_values = cbind(c_values,c_sample)
	 }


	 # Comparation of c across all the samples PER EACH LINE/ROW/TARGET
	 final_c <- c()
	 mean_rows =apply(c_values,1,mean)

	 for (i in 1:length(c_values[1,])){
		 column_partial = c_values[,i] / mean_rows
		 final_c = cbind(final_c,column_partial)
	 }

	 # Look positions with < 0.69 (deletions) or > 1.5
	 del_pos <- which( final_c < 0.69, arr.ind=T )
	 dup_pos <- which( final_c > 1.5, arr.ind=T )

         # We are going to see the info across all the samples
	 # ACHTUNG!! In sexual chromosomes, indexes are not sequential as in not sexual ==> 1,2,3....1,5,9, for instance
	 index_sexual = pmatch(sample,group)
         row_sample <- as.integer(which(del_pos[,2] == index_sexual))
         row_sample2 <- as.integer(which(dup_pos[,2] == index_sexual))

         if ((length(row_sample2) > 0) | (length(row_sample) > 0)){
		if (class(del_pos[row_sample,]) == "integer"){
			lista <- del_pos[row_sample,][1]
		}else{
			lista <- del_pos[row_sample,][,1]
		}

		if (class(dup_pos[row_sample2,]) == "integer"){
			lista2 <- dup_pos[row_sample2,][1]
		}else{
			lista2 <- dup_pos[row_sample2,][,1]
		}

		todo <- c(lista,lista2)

		aux_df <- data.frame()
		Chr <- as.integer(as.character(df$Chr[todo]))
		Start <- as.integer(as.character(df$Start[todo]))
		End <- as.integer(as.character(df$End[todo]))
		numBases <- df$End[todo] - df$Start[todo]

		
		value_sample <- final_c[todo,index_sexual]
	
		samples <- c() 			# ratio of each sample, except the one looking for

		if (length(todo) > 1){
			for (j in 1:length(todo)){
#				samples <- c(samples,paste(as.character(final_c[todo,][j,]),collapse=','))
				all_samples = as.numeric(as.character(final_c[todo,][j,]))
				all_samples = round(all_samples,digits=2)
				value_sample[j] = round(value_sample[j],digits = 2)
				to_remove =  pmatch(value_sample[j],all_samples)
				all_samples = all_samples[-to_remove]
				samples <- c(samples,paste(all_samples,collapse=','))				
			}
		}else{
			samples <- c(samples,paste(as.character(final_c[todo,]),collapse=','))
 		}

		gain_loss <- c(rep("loss",length(lista)),rep("gain",length(lista2)))
		type <- rep("CNV",sum(length(lista),length(lista2)))
		sample <- rep(sample_name,sum(length(lista),length(lista2)))
	
		aux_df <- cbind(Chr,Start,End,sample,type,numBases,gain_loss,value_sample,samples)
		aux_df <- as.data.frame(aux_df)

		return(aux_df)
	}# if length(row_sample)	
 } # if length(group) > 1
}# function



separate <- function(x){
	  aux <- strsplit(x,":")[[1]][1]
	  aux <- strsplit(aux,"chr")[[1]][2]	
}

separate_start <- function(x){
	aux <- strsplit(strsplit(x,":")[[1]][2],"-")[[1]][1]
}

separate_end <- function(x){
	aux <- strsplit(strsplit(x,":")[[1]][2],"-")[[1]][2]
}


separate_start_end <- function(x){
	aux = strsplit(x,":")[[1]][2]
}

split_cov_outcome <- function(list_file){
df <- data.frame()
index_remove <- c()

for (j in 1:length(list_file)){ 
 if (!file.exists(list_file[j])){
	write("split_cov_outcome R function: the coverage file is not given.", stderr())
 }

 tryCatch(
 {
#	    test <- read.table(file,stringsAsFactors=F,header=T)
	    test <- ReadBig(list_file[j],header=FALSE)
	    names(test) <- test[1,]
	    test <- test[-1,]

 },
	    error=function(e) {
            print(paste(e$message,list_file[j],sep=" "))
            return(NA)
 },
	
	    warning=function(cond) {
            message(cond)
 },
	    finally={
#	    test <- read.table(file,stringsAsFactors=F,header=T)
	    test <- ReadBig(list_file[j],header=FALSE)
	    #print("Read normalized file")
	    names(test) <- test[1,]
	    test <- test[-1,]

 }
 )
 
 if (j == 1){
        paseparar <- test$Target

	chr <- sapply(paseparar,separate)
	chr = gsub("Y","24",chr)
	chr = gsub("X","23",chr)
	chr = as.numeric(chr)

	start.end = sapply(paseparar,separate_start_end)
	start = sapply(paseparar,separate_start)
	end = sapply(paseparar,separate_end)

	test <- as.data.frame(append(test, list(Chr = chr, Start = start, End = end), after = 1))
	test$Chr <- as.numeric(as.character(test$Chr))
	test$Start <- as.numeric(as.character(test$Start))
	test$End <- as.numeric(as.character(test$End))
	if (all(is.na(test$gc_avg[-1])) == TRUE){
		index_remove <- c(index_remove,j)
		df = data.frame()
	}else{
		df <- data.frame(test)
	}
 }else{
	if (all(is.na(test$gc_avg[-1])) == TRUE){
		index_remove <- c(index_remove,j)
	}else if (length(df$Target) == length(test[,1])){
		df <- cbind(df,test[,c(2:5)])
	}else{
		stop("ERROR!! The bed file between coverage files is different !!")
	}
 }# if (j==1) 1st sample
} # for j
 result = list(df,index_remove)
 return (result)
} # function


######################################################################################################################

# Parameters Parsing (from Command Line)
options <- commandArgs(trailingOnly = T)

case_files <- as.character(options[1])
#control_file <- as.character(options[2])
output_path <- as.character(options[2])
bed_file <- as.character(options[3])
gender <- as.character(options[4])
controls <- as.character(options[5])
annotation.file <- as.character(options[6])

options(java.parameters = "-Xmx2000m")

#cnv_estimation_individual <- function(case_file,control_file,output_path,sample_names){
#install.packages("exactRankTests")
library(exactRankTests)
library("xlsx")

if (!file.exists(bed_file)){
	write("cnv_estimation_individual R function: bed_file is not given. The bed file is required", stderr())
}

if (controls != 0){
	control_files = strsplit(controls,";")[[1]]
}else{
	control_files = c()
}

case_files <- strsplit(case_files,";")[[1]]

 if (length(control_files) > 0){
	all_files <- c(case_files,control_files)	
	result <- split_cov_outcome(all_files)
	df.test = result[[1]]
	index_remove = result[[2]]
 }else{
	result <- split_cov_outcome(case_files)
	df.test = result[[1]]
	index_remove = result[[2]]
 } 
 

 sample_names <- c()

 gender_list <- strsplit(gender,";")[[1]]
 # We have to add the gender info for the control files (!!!)
 if (length(control_files) > 0){
	for (i in 1:length(control_files)){
		aux = strsplit(control_files[i],"_")[[1]]
		gender_list <- c(gender_list,aux[length(aux)])
	}
 }


 # if index_remove > 0 => remove the file from gender_list
 if (length(index_remove) > 0){
	 gender_list = gender_list[-index_remove]
  	 case_files = case_files[-index_remove]
 }


 # we separate the comparation between samples considering sexual chromosomes
 # ACHTUNG ! ALSO before filtering by coverage (!!) 
 chr_list <- sapply(as.character(df.test$Target),separate)
 sexual_chr = which((chr_list == "X") | (chr_list == "Y"))

 if (length(sexual_chr) > 0){
	 df.sexual <- df.test[sexual_chr,]
	 df.test <- df.test[-sexual_chr,]
 }


 # Filter those regions with avg coverage < 30 ==> when there is not enough coverage in ALL THE SAMPLES we cannot say anything 
 avg_values <- seq(from = 6, to = length(colnames(df.test)), by = 4)

 # chr1...chr22 rows
 m_aux = as.matrix(df.test[,avg_values])
 m_aux = apply(m_aux,1,as.numeric)
 m_aux = t(m_aux)
 median_avg_row <- apply(m_aux,1,median)
 less_30 = which(median_avg_row < 30)

 # sexual rows
 if (length(sexual_chr) > 0){
	 index_sexual_remove = remove_sexual(df.sexual,avg_values,gender_list)
 }

 if (length(sexual_chr) > 0){
	 out_less_30 = rbind(df.test[less_30,],df.sexual[index_sexual_remove,])
	 output_fileName = paste(output_path,"AverageCoverage_less30.tsv",sep="/")
	 tryCatch(write.table(file=output_fileName,out_less_30,quote=F,row.names=F,col.names=T,sep="\t"),error = function(e) print(paste(e$message,output_fileName,sep=" ")))
 }else{
	 out_less_30 = df.test[less_30,]
	 output_fileName = paste(output_path,"AverageCoverage_less30.tsv",sep="/")
	 tryCatch(write.table(file=output_fileName,out_less_30,quote=F,row.names=F,col.names=T,sep="\t"),error = function(e) print(paste(e$message,output_fileName,sep=" ")))
 }

 if (length(sexual_chr) > 0){ 
	 if (length(index_sexual_remove) > 0){
		df.sexual = df.sexual[-index_sexual_remove,]
	 }
 }
 df.test <- df.test[-less_30,]


 # ACHTUNG !! what if the sample is X ??	

 tryCatch(
 {
	    orig_bed <- ReadBig(bed_file,header=FALSE)
 },
	    error=function(e) {
            print(paste(e$message,bed_file,sep=" "))
            return(NA)
 },
	
	    warning=function(cond) {
            message(cond)
 },
	    finally={
	    orig_bed <- ReadBig(bed_file,header=FALSE)
 }
 )
 

 orig_bed$V1 = gsub("chrX","chr23",orig_bed$V1)
 orig_bed$V1 = gsub("chrY","chr24",orig_bed$V1)

 test_index_norm_avg <- seq(from = 7, to = length(colnames(df.test)), by = 4) 
 test_values_norm_avg <- cbind(df.test[,test_index_norm_avg])
# control_index_norm_avg <- seq(from = 8, to = length(colnames(df.control)), by = 4) 
# control_values_norm_avg <- cbind(df.control[,control_index_norm_avg])
 


 for (i in 1:length(case_files)){
        aux <- strsplit(case_files[i],"/")[[1]]
	aux2 <- gsub("_GATKcoverage.sample_interval_summary_normalized","",aux[length(aux)])
	sample_names <- c(sample_names,aux2)
 }



 #########################################################################################################################
 # Algoritmo tipo MLPA, aCGH
# wb <- createWorkbook()
 wb_significant <- createWorkbook()
 wb_filtered <- createWorkbook()

 # c_values contains for each sample IN each interval: = (num reads in interval ) / mean_reads
 c_values <- c() 
 
 for (i in 1:length(test_index_norm_avg)){
        # mean_reads contains the mean reads of each sample: sum(all reads in all the intervals) / num.intervals
	mean_reads <- sum(as.numeric(as.character(df.test[test_index_norm_avg][,i]))) / length(df.test$Target)
	c_sample = as.numeric(as.character(df.test[test_index_norm_avg][,i])) / mean_reads
	c_values = cbind(c_values,c_sample)
 }

 #c_values <- as.data.frame(c_values)
 #names(c_values) <- sample_names
 

 # Comparation of c across all the samples PER EACH LINE/ROW/TARGET
 final_c <- c()
 mean_rows =apply(c_values,1,mean)

 for (i in 1:length(c_values[1,])){
	 column_partial = c_values[,i] / mean_rows
	 final_c = cbind(final_c,column_partial)
 }

 # We save all the ratios in a separate file
 rownames(final_c) = rownames(df.test)
 colnames(final_c) = case_files
 all_ratios = as.data.frame(final_c)

 # Look positions with < 0.69 (deletions) or > 1.5
 del_pos <- which( final_c < 0.69, arr.ind=T )
 dup_pos <- which( final_c > 1.5, arr.ind=T )

 for (i in 1:length(sample_names)){
      output_fileName <- paste(sample_names[i],"_all_ratios.mut",sep="")
      output_fileName <- paste(output_path,output_fileName,sep="/")
      tryCatch(write.table(file=output_fileName,all_ratios,quote=F,row.names=T,col.names=T,sep="\t"),error = function(e) print(paste(e$message,output_fileName,sep=" ")))	
      log <- paste("Writing the output ",output_fileName,sep="")
      print(log)

      # We are going to see the info across all the samples
      row_sample <- as.integer(which(del_pos[,2] == i))
      row_sample2 <- as.integer(which(dup_pos[,2] == i))
      aux_df <- data.frame()
      if ((length(row_sample2) > 0) | (length(row_sample) > 0)){
	if (class(del_pos[row_sample,]) == "integer"){
		lista <- del_pos[row_sample,][1]
        }else{
		lista <- del_pos[row_sample,][,1]
	}

	if (class(dup_pos[row_sample2,]) == "integer"){
		lista2 <- dup_pos[row_sample2,][1]
	}else{
		lista2 <- dup_pos[row_sample2,][,1]
	}

	todo <- c(lista,lista2)

	
	Chr <- as.integer(as.character(df.test$Chr[todo]))
	Start <- as.integer(as.character(df.test$Start[todo]))
	End <- as.integer(as.character(df.test$End[todo]))
	numBases <- df.test$End[todo] - df.test$Start[todo]

			
	value_sample <- final_c[todo,i]

	other_samples <- c()
	samples <- c() 			# ratio of each sample, except the one looking for

	if (length(todo) > 1){
		for (j in 1:length(todo)){
#			samples <- c(samples,paste(as.character(final_c[todo,][j,]),collapse=','))
			all_samples = as.numeric(as.character(final_c[todo,][j,]))
			to_remove =  pmatch(value_sample[j],all_samples)
			all_samples = all_samples[-to_remove]
			all_samples = round(all_samples,digits=2)
			samples <- c(samples,paste(all_samples,collapse=','))
			other_samples <- rbind(other_samples,all_samples)
		}
	}else{
		samples <- c(samples,paste(as.character(final_c[todo,]),collapse=','))
 	}

	gain_loss <- c(rep("loss",length(lista)),rep("gain",length(lista2)))
	type <- rep("CNV",sum(length(lista),length(lista2)))
	sample <- rep(sample_names[i],sum(length(lista),length(lista2)))

	aux_df <- cbind(Chr,Start,End,sample,type,numBases,gain_loss,value_sample,samples)
	aux_df <- as.data.frame(aux_df)
	} # if (length(row_sample2)

	# JUNTAR CON LAS SALIDAS DE LOS CHR. SEXUALES ??
	if (length(sexual_chr) > 0){
		df.sexual2 = sexual_CNV(df.sexual,gender_list,i,sample_names[i])


		if (length(df.sexual2) != 0){   # if we do have results for sexual chr

			aux_df <- rbind(aux_df,df.sexual2)
			names(aux_df) <- c("Chr","Start","End","Sample","Type","numBases","gain_loss","Ratio (<0.69 loss,1 neutral,>1.5 dup)","Ratio others")

			aux_df <- aux_df[order(as.numeric(as.character(aux_df$Chr)),aux_df$Start),]

		}else{
			names(aux_df) <- c("Chr","Start","End","Sample","Type","numBases","gain_loss","Ratio (<0.69 loss,1 neutral,>1.5 dup)","Ratio others")
			aux_df <- aux_df[order(as.numeric(as.character(aux_df$Chr)),aux_df$Start),]
		}
	}else{
		names(aux_df) <- c("Chr","Start","End","Sample","Type","numBases","gain_loss","Ratio (<0.69 loss,1 neutral,>1.5 dup)","Ratio others")
		aux_df <- aux_df[order(as.numeric(as.character(aux_df$Chr)),aux_df$Start),]
	}
	
	# CNV annotation section
	if (length(aux_df) == 0){
		print("There is no CNV estimation within the bed file for")
		print(paste(sample_names[i],".mut",sep=""))    
	}else{

	# Filter the significant ones respect the other samples
	# perm.test => returns the significance level (pvalue) of a value in a distribution
	# in cases tipo: 0.5  VS. c(2.08,1.02,0.96,0.98,1,0.91,0.66,0.89) ==> pvalue is 0.11 ==> es porque tenemos un valor con un max (descompensacion por el coverage de la muestra)
	significant_filtered <- c()
	for (j in 1:length(aux_df$Chr)){
		gain_loss = aux_df$gain_loss[j]
		s_value = as.numeric(as.character(aux_df[j,8]))
		others_value = as.numeric(strsplit(as.character(aux_df[j,9]),",")[[1]])

		# normalization : std deviation = 1, mean = 0
		norm2 = c(s_value,others_value)	
		norm2 = scale(norm2,center=T,scale=T)	
		s_value2 = norm2[1]		
		others2 = norm2[-1]
		
		maximo = max(others2)
		minimo = min(others2)
		if (as.character(gain_loss) == "loss"){
#			pvalue = perm.test(s_value,others_value,alternative="less",exact=FALSE)$p.value
			pvalue = perm.test(s_value2,others2,alternative="less",exact=FALSE)$p.value
			if (pvalue > 0.05){
				index_maximo = which(others2 == maximo)
				others2 = others2[-index_maximo]
				pvalue = perm.test(s_value2,others2,alternative="less",exact=FALSE)$p.value
			}
		}else{
#			pvalue = perm.test(s_value,others_value,alternative="greater",exact=FALSE)$p.value
			pvalue = perm.test(s_value2,others2,alternative="greater",exact=FALSE)$p.value
			if (pvalue > 0.05){
				index_minimo = which(others2 == minimo)
				others2 = others2[-index_minimo]
				pvalue = perm.test(s_value2,others2,alternative="greater",exact=FALSE)$p.value
			}		
		}

		if ((pvalue != "NaN") & (pvalue <= 0.05)){
			significant_filtered <- c(significant_filtered,j)
		}
	}

	# JULIO 2014: antes veiamos que los targets del cnv estaban en el bed ==> lo estan, porque GATK llama con ese bed (antes: habia un lio de bed's)
	# AHORA: como en el coverage_control script R, anotaremos GENE, EXONS. pero ACHTUNG!! de manera diferente
	# porque ahora, podemos tener casos como en karen(DE), que tiene regiones/intervalos intrónicos


	tryCatch(
	     {
	       annot <- ReadBig(annotation.file,header=FALSE)
	     },

	     error=function(e) {
		 print(paste(e$message,annotation.file,sep=" "))
		 return(2)
	     },

	     warning=function(cond) {
		   message(cond)
	     },
	     finally={
		  annot <- ReadBig(annotation.file,header=FALSE)
		  colnames(annot) = annot[1,]
		  annot = annot[-1,]
	    }
	)

	# aux_df no debe ser muy grande => la recorremos y examinamos
	gene_info = c()
	exon_info = c()
	refseq_info = c()
	for (j in 1:length(aux_df$Chr)){
		if (as.character(aux_df$Chr[j]) == "23"){
			partial_chr = annot[which(annot$Chr == paste("chr","X",sep="")),]	
		}else if (as.character(aux_df$Chr[j]) == "24"){
			partial_chr = annot[which(annot$Chr == paste("chr","Y",sep="")),]	
		}else{
			partial_chr = annot[which(annot$Chr == paste("chr",as.numeric(as.character(aux_df$Chr[j])),sep="")),]
		}

		extreme_izdo = which(as.numeric(partial_chr$exonStart) >= as.numeric(as.character(aux_df$Start[j])))
		extreme_dcho = which(as.numeric(partial_chr$exonEnd) <= as.numeric(as.character(aux_df$End[j])))
		index = intersect(extreme_izdo,extreme_dcho)

		# Targets within the bed that include more than one exon
		if (length(index) > 0){
			if (length(index) > 1){
				for (l in 1:length(index)){
					if (l == 1){
						aux_gene = partial_chr$Gene[l]
						aux_exon = partial_chr$exonCount[l]
						aux_refseq = partial_chr$Refseq[l]
					}else{
						aux_gene = paste(aux_gene,partial_chr$Gene[l],sep=",")
						aux_exon = paste(aux_exon,partial_chr$exonCount[l],sep=",")
						aux_refseq = paste(aux_refseq,partial_chr$Refseq[l],sep=",")
					}
				}
				gene_info = c(gene_info,aux_gene)
				exon_info = c(exon_info,aux_exon)
				refseq_info = c(refseq_info,aux_refseq)
			}else{
				gene_info = c(gene_info,partial_chr$Gene[index])
				exon_info = c(exon_info,partial_chr$exonCount[index])
				refseq_info = c(refseq_info,partial_chr$Refseq[index])
			}

		# targets with its corresponding individual target in the bed
		}else{
			aux = "no_info"
			gene_info = c(gene_info,aux)
			exon_info = c(exon_info,aux)
			refseq_info = c(refseq_info,aux)						

		}
	}

		
	# aux_df => dataframe which contains all the "loss" and "gain" targets. ACHTUNG! not all are significant!!
	# aux_df_significant => dataframe which contains all the "loss" and "gain" targets BUT being significant (!!!) --> significant_filtered
        # Filter the resulted CNV > 1kbp (teoretically) => > 500bp + significantly
	
	 # TODOS DA IGUAL SIENDO SIGNIFICATIVO ENTRE LOS DEMAS
	 target <- paste(orig_bed$V2,orig_bed$V3,sep="-")
	 target <- paste(orig_bed$V1,target,sep=":")

	 new_chr = paste("chr",as.character(aux_df$Chr),sep="")
	 df_out <- cbind(new_chr,as.character(aux_df$Start),as.character(aux_df$End),gene_info,exon_info,refseq_info,as.character(aux_df$Sample),as.character(aux_df$Type),as.character(aux_df$numBases),as.character(aux_df$gain_loss),as.character(aux_df[,8]),as.character(aux_df[,9]))
	 df_out <- as.data.frame(df_out)
 	 names(df_out) <- c("Chr","Start","End","Gene","Exon","Refseq","Sample","Type","numBases","gain_loss","Ratio (<0.69 loss,1 neutral,>1.5 dup)","Ratio others")
	
	 df_out$Chr = gsub("chr23","chrX",df_out$Chr)
	 df_out$Chr = gsub("chr24","chrY",df_out$Chr)

	 rm(aux_df)
	 aux_df <- df_out
	 rm(df_out)

	 # LOS SIGNIFICATIVOS
	 #cuales <- which(significant_filtered %in% indexes_filtered)
         #significant_filtered_rows <- significant_filtered[cuales]    # <== esto lo haciamos despues de mirar que targets del bed estaban en el original


	 # Format aux_df in order to represent better te data
	 aux_df2 = aux_df
	 # 0. as.numeric()
	 # 1. Round
	 # 2. as.character
	 # 3. escritura
	 aux_df2[,11] = as.numeric(as.character(aux_df2[,11]))
	 aux_df2[,11] = round(aux_df2[,11],digits = 2)
	 aux_df2[,11] = as.character(aux_df2[,11])
	 aux_df2[,11] = gsub("\\.",",",aux_df2[,11])
	
#	 aux_df_significant <- aux_df2[significant_filtered_rows,]
	 aux_df_significant <- aux_df2[significant_filtered,]

	 output_fileName <- paste(sample_names[i],"_CNV_all.mut",sep="")
	 output_fileName <- paste(output_path,output_fileName,sep="/")
	 tryCatch(write.table(file=output_fileName,aux_df2,quote=F,row.names=F,col.names=T,sep="\t"),error = function(e) print(paste(e$message,output_fileName,sep=" ")))	
	 log <- paste("Writing the output ",output_fileName,sep="")
	 print(log)

	 output_fileName <- paste(sample_names[i],"_CNV_significant.mut",sep="")
	 output_fileName <- paste(output_path,output_fileName,sep="/")
	 tryCatch(write.table(file=output_fileName,aux_df_significant,quote=F,row.names=F,col.names=T,sep="\t"),error = function(e) print(paste(e$message,output_fileName,sep=" ")))	
	 log <- paste("Writing the output ",output_fileName,sep="")
	 print(log)

	# write into excel
	sheet_significant <- createSheet(wb_significant,sheetName=sample_names[i])
	# add the data to the new sheet
        addDataFrame(aux_df_significant, sheet_significant,row.names=F)	


	# Once we have all intervals "loss" and "gain" (aux_df) and those that are significant (aux_df_significant) ==> we want to take those major than 500bp being consecutive within the bed file
	# That is, 'gain' and 'loss'-es that are solapped

	# Iremos por chr, porque cada uno tiene un start-end ==> ordenandolo y mirando con aux_df_significant

	aux_df_significant$Chr = gsub("chrX","chr23",aux_df_significant$Chr)
        aux_df_significant$Chr = gsub("chrY","chr24",aux_df_significant$Chr)


	# We filter those SIGNIFICANT copy number variations > 500 BUT ONLY IF THERE ARE SIGNIFICANT ONES!!!

	if (dim(aux_df_significant)[1] > 0 ){
		unique_chr = unique(aux_df_significant$Chr)			# chr within the aux_df_significant info
		aux_df_filtered = c()
		for (z in 1:length(unique_chr)){
			# separar por chr
			bed_chr = orig_bed[which(orig_bed$V1 == unique_chr[z]),]
			aux_df_chr = aux_df_significant[which(as.character(aux_df_significant$Chr) == unique_chr[z]),]

			# ordenar por start cada chr por separado
			bed_chr_ordered = bed_chr[order(bed_chr$V2),]
			aux_df_chr_ordered = aux_df_chr[order(as.numeric(as.character(aux_df_chr$Start))),]

			# Ahora ver si aux_df_chr_ordered$Start siguen intervalos seguidos en el bed_chr_ordered$Start
			indices = findInterval(as.numeric(as.character(aux_df_chr_ordered$Start)),bed_chr_ordered$V2)
			breaks <- c(0, which(diff(indices) != 1), length(indices))
	     	        lista <- sapply(seq(length(breaks) - 1), function(x) indices[(breaks[x] + 1):breaks[x+1]])

			indice_salida <- 1
			filtered <- c()
	 
			if (length(lista) > 1){
				for (j in 1:length(lista)){
				     sum <- 0
				     if (length(lista[[j]]) > 1){
					for (l in 1:length(lista[[j]])){
					     sum <- sum + bed_chr_ordered[lista[[j]][l],]$V3 - bed_chr_ordered[lista[[j]][l],]$V2
					}
					gol <- c() # array with "gain" or "loss" 
					same_chr <- c()    # array with the "chr"
					indice_aux <- indice_salida
					for (l in 1:length(lista[[j]])){
						gol <- c(gol,as.character(aux_df$gain_loss[indice_aux]))
						same_chr <- c(same_chr,as.character(aux_df$Chr[indice_aux]))
						indice_aux <- indice_aux + 1
					}
			
				    	if ((sum >= 500) & (length(unique(gol)) == 1) & (length(unique(same_chr)) == 1)){    # if sum total >= 500 and all are "gain" or "loss" and in the same chr
				 	   aux <- indice_salida + length(lista[[j]]) -1
				 	   interval <- seq(from=indice_salida,to=aux,by=1)
				 	   filtered <- c(filtered,interval)
			    		}	
				    }else{
					# If the orig_bed file is small and there is only one entrance for a target (ha pasado una vez con OI-muestra4, de julio) ==>
					sum <- bed_chr_ordered[lista[[j]],]$V3 - bed_chr_ordered[lista[[j]],]$V2
					if (sum >= 500){
					      filtered <- c(filtered,indice_salida)
					}	
		     		   }
				 indice_salida <- indice_salida + length(lista[[j]]) 
				}# for
			}else{
				if (lista == 0){
					sum <- bed_chr_ordered[1,]$V3 - bed_chr_ordered[1,]$V2
					if (sum >= 500){
					      filtered <- c(filtered,indice_salida)
					}			
				}else{
					sum <- bed_chr_ordered[lista,]$V3 - bed_chr_ordered[lista,]$V2
					if (sum >= 500){
					      filtered <- c(filtered,indice_salida)
					}	
				}
			}
			if (length(filtered > 0)){
				aux_df_filtered = rbind(aux_df_filtered,aux_df_chr_ordered[filtered,])
			}

	   	   }# for unique_chr

		   aux_df_filtered$Chr = gsub("chr23","chrX",aux_df_filtered$Chr)
		   aux_df_filtered$Chr = gsub("chr24","chrY",aux_df_filtered$Chr)
	}else{
		aux_df_filtered = colnames(aux_df2)

	} # if dimension > 0


	   output_fileName <- paste(sample_names[i],"_CNV_filtered_significant.mut",sep="")
	   output_fileName <- paste(output_path,output_fileName,sep="/")
	   tryCatch(write.table(file=output_fileName,aux_df_filtered,quote=F,row.names=F,col.names=T,sep="\t"),error = function(e) print(paste(e$message,output_fileName,sep=" ")))
	   log <- paste("Writing the output ",output_fileName,sep="")
	   print(log)	

	   # write into excel	
	   sheet_filtered <- createSheet(wb_filtered, sheetName=sample_names[i])
	   # add the data to the new sheet
	   addDataFrame(aux_df_filtered, sheet_filtered,row.names=F)	

      } # else annotation section
 }# for i (sample_names)

# saving the excel
#saveWorkbook(wb, paste(output_path,"CNV_all.xlsx",sep="/"))
saveWorkbook(wb_filtered, paste(output_path,"CNV_filtered_significant.xlsx",sep="/"))
saveWorkbook(wb_significant, paste(output_path,"CNV_significant.xlsx",sep="/"))