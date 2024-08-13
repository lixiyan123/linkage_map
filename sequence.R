install.packages("VariantAnnotation")
install.packages("Biostrings")
library(Biostrings)
library(VariantAnnotation)
library(openxlsx)

library(Biostrings)
library(VariantAnnotation)
library(openxlsx)

# ?????????????????????
genome_file <- "E:/NEWRNA/RNASEQ/hap/f153_v7.fasta"
genome <- readDNAStringSet(genome_file)
ref_genome <- genome[["16"]]
head(ref_genome)
# ??????VCF??????
vcf_file <- "E:/NEWRNA/RNASEQ/hap/102.vcf.vcf"
vcf <- readVcf(vcf_file, "f153v7")  

# ????????????
chromosome <- "16"
start_position <- 10396106 - 500
end_position <- 10396106 + 550
chromosome1 <- "contig16"

# ??????????????????VCF??????
subset_vcf <- vcf[seqnames(vcf) == chromosome1 &
                    start(vcf) >= start_position &
                    end(vcf) <= end_position]

# ?????????????????????
sample_names <- colnames(vcf)

print(sample_names)
# ??????????????????????????????????????????????????????
samples <- rowData(vcf)$ID
ref_sequence <- unlist(strsplit(as.character(subseq(genome[[chromosome]], start_position, end_position)), ""))
# ??? ref_sequence ???????????????????????????
sequence_string <- paste(ref_sequence, collapse = "")

# ??? sequence_string ????????? txt ?????????????????????????????????
write(sequence_string, file = "ref_sequence.txt")


getwd()


# ??????????????????????????? _sequence.fasta ??????????????????
file_names <- list.files(pattern = "_sequence.fasta")

# ????????????????????????????????????????????????
merged_content <- character()

# ???????????????????????????????????????????????????????????? merged_content ???
for (file in file_names) {
  # ???????????????????????????
  file_content <- readLines(file)
  
  # ???????????? "Aco013395_" ???????????????????????????
  prefixed_content <- paste0(file_content)
  
  # ?????????????????? merged_content ???
  merged_content <- c(merged_content, prefixed_content)
}

# ????????????????????????????????????????????????
writeLines(merged_content, "merged_sequences_with_prefix.fasta")

prefixed_content
file_content[5]

alt_allele <- alt(subset_vcf)
ref_allele <- ref(subset_vcf)
allele_counts <- sapply(alt_alleles, length)

width(alt_alleles)[[1]]
print( ref_alleles,max = 1000)
alt_allele <- alt(subset_vcf)[112]
allele[2]
head(alt_allele[])
head(subset_vcf)
first_alt_value <- strsplit(as.character(subset_vcf$ALT[1]), ",")[[1]][1]
first_alt_value <- alt(subset_vcf)[1]
print( ref_alleles[])
ref_alleles[2]
sample_sequence1[replace_position-2] 

length(alt_alleles[[2]])
sample_sequence1[replace_position] <as.character(alt_alleles[[1]])
# ??? alt_allele ???????????????????????????
write.table(ref_allele, "ref_allele_data.txt", sep = "\t", quote = FALSE)
for (i in 1:length(alt_alleles[[1]])) {
  sample_sequence1[replace_position] <- as.character(alt_alleles[[1]][i])
  replace_position <- replace_position + 1
}
str(first_alt_value)
sample_sequence[2]<- as.character(first_alt_value[[1]])
sample_sequence[position - start_position + 1] <- as.character(alt_allele[[1]])
allele <- subset_vcf[394, '1-1']
head(geno(subset_vcf)) 
head(geno(subset_vcf)$GT)
genotype_info <- geno(subset_vcf)$GT[sample_name]
# ?????? variant ?????????????????????????????????1 ??? 37 ????????????????????????sample_name ???????????????
# ?????????????????????variant???????????????????????????
genotype_info <- geno(subset_vcf)$GT
# ?????? geno(subset_vcf) ??? GT ??????????????????
str(geno(subset_vcf)$GT)

# ????????????????????? head() ???????????????????????????
head(geno(subset_vcf)$GT)
genotype_info <- geno(subset_vcf)$GT[variant, sample_name]
genotype_info <- geno(subset_vcf)$GT[1,'1-13']####?????????
sample_genotype <- as.character(subset_vcf[[1-13]])
length(start(subset_vcf))
subset_vcf[1-1]
start(subset_vcf)[1]
genotype_info <- geno(subset_vcf)$GT[variant, sample_name]
allele <- strsplit(geno(subset_vcf)$GT[variant, 'sample_name'], "|")[[1]]
ref_alt <- with(reftrack(vcf), data.frame(REF, ALT))
print(ref_genome)
cumulative_difference <- 0
length(ref_alleles)[[1]]
ref_alleles <- ref(subset_vcf)[[20]]
alt_alleles <- alt(subset_vcf)[[20]]
start(subset_vcf)[5]-start_position +1 
geno(subset_vcf)$GT[2,'1-1']
geno(subset_vcf)$GT


ref(subset_vcf)[] alt(subset_vcf)
# ??? alt(subset_vcf) ????????? DNAStringSet ?????????????????????????????????
alt_data <- sapply(alt(subset_vcf), as.character)

# ??? ref(subset_vcf) ????????? DNAStringSet ?????????????????????????????????
ref_data <- sapply(ref(subset_vcf), as.character)

# ??? geno(subset_vcf) ????????????????????????
geno_data <- as.character(geno(subset_vcf)$GT[,'1-1'])

# ????????????
combined_data <- cbind(geno_data, ref_data, alt_data)

# ????????????????????????
output_file <- "output.txt"

# ??????????????????
write.table(combined_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

getwd()
geno(subset_vcf)$GT
# ????????????
genename <- "in10396106R2"
dir.create(genename)
chromosome <- "16"
ref_genome <- genome[["16"]]

start_position <- 10396106 - 2500
end_position <- 10396106 + 2500


chromosome1 <- "contig16"
# ??????????????????VCF??????
subset_vcf <- vcf[seqnames(vcf) == chromosome1 &
                    start(vcf) >= start_position &
                    end(vcf) <= end_position]
################################################


for (sample_name in sample_names) {
  # ?????????????????????
  sample_sequence1 <- character(end_position - start_position + 1)
  sample_sequence2 <- character(end_position - start_position + 1)
  cumulative_difference1 <- 0 
  cumulative_difference2 <- 0 
  replace_position1 <- 0
  replace_position2 <- 0
  variant <- 0
  ref_alleles <- 0
  for (position in seq(start_position, end_position)){
      
    
       
    if (position %in% start(subset_vcf)) {
      variant <-  variant +1
      allele <- strsplit(geno(subset_vcf)$GT[variant, sample_name], "|")[[1]]
      ref_alleles <- ref(subset_vcf)[[variant]]
      alt_alleles <- alt(subset_vcf)[[variant]]
      
      # ???????????????
      if (allele[1] == "0" && allele[3] == "0")
        {
        if (length(ref_alleles) == 0) {
          replace_position1 <- replace_position1 + 1
          replace_position2 <- replace_position2 + 1
          sample_sequence1[replace_position1] <- "N" 
          sample_sequence2[replace_position2] <- "N" 
        } 
        else if (length(ref_alleles) > 0) { 
          
          for (a in 1:length(ref_alleles)) {
            replace_position1 <- replace_position1 + 1
            sample_sequence1[replace_position1] <- as.character(ref_alleles[a]) 
            replace_position2 <- replace_position2 + 1
            
            sample_sequence2[replace_position2] <- as.character(ref_alleles[a])
          }
        }
      } 
      else if (allele[1] == "0" && allele[3] == "1") {
        # ?????????????????? "0|1" ?????????
        if (length(ref_alleles)[[1]] == 0 && width(alt_alleles)[[1]] == 0)
          {
          replace_position1 <- replace_position1 + 1
          replace_position2 <- replace_position2 + 1
          sample_sequence1[replace_position1] <- "N" 
          sample_sequence2[replace_position2] <- "N" 
          
          
        }
        if(length(ref_alleles)[[1]] == 0 && width(alt_alleles)[[1]] > 0)
        {
          replace_position1 <- replace_position1 + 1
          
          sample_sequence1[replace_position1] <- "N" 
  
          
          # cumulative_difference2 <- cumulative_difference2 + difference2
          
          for (a in 1:width(alt_alleles)[[1]]) {
            replace_position2 <- replace_position2 + 1
            sample_sequence2[replace_position2] <- as.character(alt_alleles[[1]][a])
          
          }
       }
        
        if(length(ref_alleles)[[1]] > 0 && width(alt_alleles)[[1]] == 0)
        {  
          # difference2 <- length(ref_alleles)[[1]] - 1
          # # cumulative_difference2 <- cumulative_difference2 + difference2
          replace_position2 <- replace_position2 + 1
          sample_sequence2[replace_position2] <- "N" 
          
          for (a in 1:length(ref_alleles)[[1]]) {
            replace_position1 <- replace_position1 + 1
            sample_sequence1[replace_position1] <- as.character(ref_alleles[a])
            
          }
        }
        
        
        if(length(ref_alleles)[[1]] > 0 && width(alt_alleles)[[1]] > 0)
        {
          # difference1 <- length(ref_alleles)[[1]] - 1
          # # cumulative_difference1 <- cumulative_difference1 + difference1
          # difference2 <- width(alt_alleles)[[1]] - 1
          # # cumulative_difference2 <- cumulative_difference2 + difference2
          # 
          for (a in 1:length(ref_alleles)[[1]]) {
            replace_position1 <- replace_position1 + 1
            sample_sequence1[replace_position1] <- as.character(ref_alleles[a])
            
          }   
          for (a in 1:width(alt_alleles)[[1]]) {
            replace_position2 <- replace_position2 + 1
            sample_sequence2[replace_position2] <- as.character(alt_alleles[[1]][a])
            
          }
          
        }
        
      }
      else if (allele[1] == "1" && allele[3] == "1") {
        # ?????????????????? "1|1" ?????????
        if ( width(alt_alleles)[[1]] == 0) {
          replace_position1 <- replace_position1 + 1
          replace_position2 <- replace_position2 + 1
          sample_sequence1[replace_position1] <- "N" 
          sample_sequence2[replace_position2] <- "N" 
          
        }
        else if  (width(alt_alleles)[[1]] > 0) {  
          # difference1 <- width(alt_alleles)[[1]] - 1
          # # cumulative_difference1 <- cumulative_difference1 + difference1
          # difference2 <- width(alt_alleles)[[1]] - 1
          # # cumulative_difference2 <- cumulative_difference2 + difference2
          for (a in 1:width(alt_alleles)[[1]]) {
            replace_position1 <- replace_position1 + 1
            sample_sequence1[replace_position1] <- as.character(alt_alleles[[1]][a])
            
            replace_position2 <- replace_position2 + 1
            sample_sequence2[replace_position2] <- as.character(alt_alleles[[1]][a])
            
           
          } 
      }
        
      }
      
      else if (allele[1] == "1" && allele[3] == "0") {
        # ?????????????????? "1|0" ?????????
        if (length(ref_alleles)[[1]] == 0 && width(alt_alleles)[[1]] == 0 )
          {
          replace_position1 <- replace_position1 + 1 
          replace_position2 <- replace_position2 + 1 
          sample_sequence1[replace_position1] <- "N" 
          sample_sequence2[replace_position2] <- "N" 
          
        }
        if(length(ref_alleles)[[1]] == 0 && width(alt_alleles)[[1]] > 0)
        {
          replace_position2 <- replace_position2 + 1
          sample_sequence2[replace_position2] <- "N" 
          
          # difference1 <- width(alt_alleles)[[1]] - 1
          
          # cumulative_difference1 <- cumulative_difference1 + difference1
          
          for (a in 1:width(alt_alleles)[[1]]) {
            replace_position1 <- replace_position1 + 1
            sample_sequence1[replace_position1] <- as.character(alt_alleles[[1]][a])
          
          }
        }
        if(length(ref_alleles)[[1]] > 0 && width(alt_alleles)[[1]] == 0)
        {  
          # difference2 <- length(ref_alleles)[[1]] - 1
          # # cumulative_difference2 <- cumulative_difference2 + difference2
          replace_position1 <- replace_position1 + 1
          sample_sequence1[replace_position1] <- "N" 
          
          for (a in 1:length(ref_alleles)[[1]]) {
            replace_position2 <- replace_position2 + 1
            sample_sequence2[replace_position2] <- as.character(ref_alleles[a])
            
          }
        }
          
       
        if(length(ref_alleles)[[1]] > 0 && width(alt_alleles)[[1]] > 0)
        {
          # difference1 <- length(ref_alleles)[[1]] - 1
          # # cumulative_difference1 <- cumulative_difference1 + difference1
          # difference2 <- width(alt_alleles)[[1]] - 1
          # # cumulative_difference2 <- cumulative_difference2 + difference2
          # 
          for (a in 1:length(ref_alleles)[[1]]) {
            replace_position2 <- replace_position2 + 1
            sample_sequence2[replace_position2] <- as.character(ref_alleles[a])
            
          }   
          for (a in 1:width(alt_alleles)[[1]]) {
            replace_position1 <- replace_position1 + 1
            sample_sequence1[replace_position1] <- as.character(alt_alleles[[1]][a])
           
          }
          
      }
      
      
    
      
    }
    }  
    else {
      if ( length(ref_alleles) - 1  > 0  ){
        skipped <- 0
        ref_alleles <- ref_alleles[-length(ref_alleles)]
        skipped <- skipped + 1
        next
      }
      else {
          replace_position1 <- replace_position1  + 1
          replace_position2 <- replace_position2  + 1
          sample_sequence1[replace_position1] <- as.character(ref_genome[position])
          sample_sequence2[replace_position2]  <- as.character(ref_genome[position])
        }
    }
}
  # 
  sample_sequence1 <- gsub("N", "", sample_sequence1)
   # # sample_sequence1 <- sample_sequence1[!is.na(sample_sequence1)]
  sample_sequence2 <- gsub("N", "", sample_sequence2)
  # sample_sequence2 <- sample_sequence2[!is.na(sample_sequence2)]
 
  sample_sequence_file <- paste0(sample_name, "_sequence.fasta")
  # ?????????????????????????????????????????????
  output_path <- file.path(genename, sample_sequence_file)
  
  write(paste(">10396106_allele1_",sample_name,"_sequence1\n",paste(sample_sequence1,collapse =""),
              "\n>10396106_allele2_",sample_name,"_sequence2\n",paste(sample_sequence2,collapse ="")), file = output_path) 
}

setwd("E:/NEWRNA/RNASEQ/in10396106R2")
sample1 <- "10396106"
getwd()
# ??????????????????????????? _sequence.fasta ??????????????????
file_names <- list.files(pattern = "_sequence.fasta")
# ????????????????????????????????????????????????
merged_content <- character()
# ???????????????????????????????????????????????????????????? merged_content ???
for (file in file_names) {
# ???????????????????????????
file_content <- readLines(file)
prefixed_content <- paste0(file_content)
# ?????????????????? merged_content ???
merged_content <- c(merged_content, prefixed_content)
}
# ????????????????????????????????????????????????
writeLines(merged_content, paste0("merged_sequences_", sample1, ".fasta"))

