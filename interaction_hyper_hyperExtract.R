
#my own interaction count script for Adam's data analysis. Specific for Jenner
#run in R_info directory
#setwd("/home/andrew/Documents/Emma_RNA/adam/publication/R")
#setwd("../../Adam")
#install.packages("stringr", lib="R_info/library", repos = "http://cran.us.r-project.org")
#install.packages("tidyverse", lib="R_info/library", repos = "http://cran.us.r-project.org")
library("stringr", lib.loc="R_info/library") #to allow string_trim later on

library("backports", lib.loc="R_info/library")
library("vctrs", lib.loc="R_info/library")
library("crayon", lib.loc="R_info/library")
library("pillar", lib.loc="R_info/library")
library("readr", lib.loc="R_info/library")
library("tibble", lib.loc="R_info/library")
library("tidyr", lib.loc="R_info/library")
args = commandArgs(trailingOnly=TRUE)
#install.packages("readr", lib = "/home/andrew/Documents/Emma_RNA/adam/publication/R/library", repos = "http://cran.us.r-project.org")
annotation <- as.data.frame(read_delim("R_info/complete_annotation.tsv",
                                       "\t", escape_double = FALSE, col_names = FALSE,
                                       trim_ws = TRUE))
filename=args[1]
input=paste0("featureCounts_output/",filename,"_chimera_Chimeric_Aligned_single_merged_.out.sam.featureCounts")
print(input)
feature_counts <- as.data.frame(read_delim(input,
                                           "\t", escape_double = FALSE, col_names = FALSE,
                                           trim_ws = TRUE))


# Features #1 and #10191 have the same name ("new_4215473_4215670"), because the feature spans from the
# 'end' of the genome to the 'start'. Feature #1 is removed from the feature list when counting 
# interactions, so that reads mapped to both features #1 and #10191 are counted for #10191.

features <- annotation[-1,13]

interaction_counts <- matrix(0,ncol = length(features), nrow = length(features))

#colnames(interaction_counts) <- features
#rownames(interaction_counts) <- features

print("start aggregation")
feature_counts <- as_tibble(feature_counts)
feature_counts=feature_counts[!duplicated(feature_counts),]
featurs_per_read <- as.data.frame(aggregate(feature_counts$X4~feature_counts$X1,FUN = toString))

output_file=paste0("featureCounts_output/",filename,"_chimera_Chimeric_Aligned_single_merged_aggregated.out.sam.featureCounts")
write.table(featurs_per_read, file=output_file,quote = FALSE,sep = '\t', row.names = FALSE, col.names = FALSE)

y <- (lapply(featurs_per_read$`feature_counts$X4`, function(x)
  unlist(strsplit(as.character(x),","))))
print("start for loop")
for (i in 1:length(y)){
  
  z <- y[i]
  z0=unlist(z)
  z0=unique(z0)
  print(z0)
  print(paste0(i,"/",length(y)))
  for(j in 1:length(z0)){
    a <- str_trim(z0[j]) 
    row_index=which(features==a)
    for(k in 1:length(z0)){
      #currently z is in a list with multiple elements. The unlist function turns those elements into vectors
      #str_trim are added make sure no unexpectected white space were included 
      # print(a)
      
      if ((i != j) | (length(z0)==1)){ #to include non-chimeric read count, while avoiding multiple counting
        b <- str_trim(z0[k])
        # print(b)
        #print(row_index)
        col_index=which(features==b)
        # print(col_index)
        interaction_counts[row_index,col_index] <- interaction_counts[row_index,col_index] + 1
      }
    }
  }
}
print("for loop completed")
output=paste0("interaction_counts/",filename,"_chimera_Aligned_Single_merged_Chimeric.tsv")
write.table(interaction_counts, file=output,quote = FALSE,sep = '\t', row.names = FALSE, col.names = FALSE)


#================================================================================================

name=filename
annotation <- as.data.frame(read_delim("R_info/complete_annotation.tsv",
                                       "\t", escape_double = FALSE, col_names = FALSE,
                                       trim_ws = TRUE))
features <- annotation[-1,13]
#name_vector=rep(NA,length(list_of_files))
#number=1
#for (files in list_of_files){
#names=unlist(strsplit(files,".", fixed = TRUE))[1]
#name_vector[number]=names
#print(name_vector)
# number=number+1
#}  

for (h in 1:1){
  print(h)
  input=output
  print(input)
  int_counts <- as.data.frame(read_delim(input,
                                         "\t", escape_double = FALSE, col_names = FALSE,
                                         trim_ws = TRUE))
  
  
  # int_counts=int_counts[,-1]
  #int_counts=int_counts[-1,]
  # rownames(int_counts)=features
  #  colnames(int_counts)=features
  print("Read interaction counts file")
  thr.pv = 0.05
  print(dim(int_counts)[1])
  print(dim(int_counts)[1])
  pv.feature.feature.int <- matrix(1,nrow=dim(int_counts)[1],ncol=dim(int_counts)[1])
  int_matrix=data.matrix(int_counts, rownames.force = NA)

  
  n0=sum(int_matrix)
  #for(i in 1: 100) { 
  for(i in 1: dim(int_matrix)[1]) {
    # print("start for loop")
    m = sum(int_matrix[i,])
    n = n0 - m
    # print("test begin")
    for(j in 1: dim(int_matrix)[1]) {
      k <- sum(int_matrix[,j])
      x <- int_matrix[i,j]
      pv.feature.feature.int[i,j] <- phyper(x-1, m, n, k,lower.tail = F)
      
    }
    # print("test completed")
  }
  #  rownames(pv.feature.feature.int) = rownames(int_counts)
  # colnames(pv.feature.feature.int) = colnames(int_counts)
  hyper_output=paste0("hypergeometric_output/",name,"_chimera_Aligned_Single_merged.tsv")
  print(output)
  write.table(pv.feature.feature.int, file=hyper_output,quote = FALSE,sep = '\t', row.names = FALSE, col.names = FALSE)
  hyper.padj = matrix(p.adjust(as.numeric(unlist(pv.feature.feature.int)),method = "BH"),nrow=nrow(pv.feature.feature.int),ncol=ncol(pv.feature.feature.int))
  hyper_padj_output=paste0("hypergeometric_output/",name,"_chimera_Aligned_Single_merged_padj.tsv")
  write.table(hyper.padj, file=hyper_padj_output,quote = FALSE,sep = '\t', row.names = FALSE, col.names = FALSE)
  print("done")
}
#================================================================================================

input=as.data.frame(read_delim(hyper_padj_output,
                               "\t", escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE))

hyper <- as.data.frame(sapply(input, as.numeric))
#add the gtf annotation as row names
#annotation <- as.data.frame(read_delim("R_info/complete_annotation.tsv",
                      #                 "\t", escape_double = FALSE, col_names = FALSE,
                    #                   trim_ws = TRUE))
#features <- annotation[-1,13]
features_property <- annotation[-1,3]
features_start <- annotation[-1,4]
features_end <- annotation[-1,5]
features_strand <- annotation[-1,7]

#colnames(hyper)=features
#rownames(hyper)=features
#hyper=hyper[-1,]
#pick out row and columns with p-values below 0.05
highlighted_index=which(hyper<=0.05,arr.ind = T)
#write into new files. Include their start and end position of each faatures
#create empty vector
partner1_vec=c()
partner2_vec=c()
partner1_property_vec=c()
partner2_property_vec=c()
partner1_start_vec=c()
partner2_start_vec=c()
partner1_end_vec=c()
partner2_end_vec=c()
partner1_strand_vec=c()
partner2_strand_vec=c()
for (i in 1:nrow(highlighted_index)){
  partner1=features[highlighted_index[i,1]]
  partner1_vec=c(partner1_vec,partner1)
  partner2=features[highlighted_index[i,2]]
  partner2_vec=c(partner2_vec,partner2)
  
  partner1_property=features_property[highlighted_index[i,1]]
  partner1_property_vec=c(partner1_property_vec,partner1_property)
  partner2_property=features_property[highlighted_index[i,2]]
  partner2_property_vec=c(partner2_property_vec,partner2_property)
  
  
  partner1_start=features_start[highlighted_index[i,1]]
  partner1_start_vec=c(partner1_start_vec,partner1_start)
  partner2_start=features_start[highlighted_index[i,2]]
  partner2_start_vec=c(partner2_start_vec,partner2_start)
  
  
  partner1_end=features_end[highlighted_index[i,1]]
  partner1_end_vec=c(partner1_end_vec,partner1_end)
  partner2_end=features_end[highlighted_index[i,2]]
  partner2_end_vec=c(partner2_end_vec,partner2_end)
  
  partner1_strand=features_strand[highlighted_index[i,1]]
  partner1_strand_vec=c(partner1_strand_vec,partner1_strand)
  partner2_strand=features_strand[highlighted_index[i,2]]
  partner2_strand_vec=c(partner2_strand_vec,partner2_strand)
}
#output_data=cbind(partner1_vec,partner2_vec,partner1_property_vec,  partner2_property_vec, partner1_start_vec, partner2_start_vec, partner1_end_vec,partner2_end_vec, partner1_strand_vec,partner2_strand_vec)
output_data=cbind(partner1_vec,partner1_property_vec,partner1_start_vec,partner1_end_vec,partner1_strand_vec,partner2_vec,  partner2_property_vec,  partner2_start_vec, partner2_end_vec, partner2_strand_vec)
output_name=paste0("hypergeometric_output/",filename,"_chimera_Aligned_Single_merged_padj_extracted.tsv")
write.table(output_data, file=output_name,quote = FALSE,sep = '\t', row.names = FALSE, col.names = TRUE)

