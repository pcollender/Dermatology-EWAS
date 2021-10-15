#Main difference here is that we do not want to remove low detect p loci or do any normalization
##We DO want to remove the 'tube' sample
library(pryr) # for monitoring memory use
library(matrixStats) # for calculating summary statistics
library(reshape)
library(scales)
library(shinyMethyl)
library(minfi)
library(ENmix) # probe type adjustment "r
library(limma) # for MDS plots
library(sva) # for addressing batch effects
library(RColorBrewer)

setwd('Dermatology/')
#' importing idat files, result is a RGChannelSet

sheet <- read.metharray.sheet(getwd(), pattern = "csv$")
sheet$recorded.sex = substr(sapply(strsplit(sheet$Sample_Name,'_'),'[',2),1,1)

#table(sheet$Slide)
#source('../Berkeley-CHAMACOS-methylation/readmetharrayexp_fixed.R')
#WB <- read.metharray.exp(targets = sheet, recursive = T, verbose = T)
#dir.create('Processed')
#save(WB, file = 'Processed/WB_obj.Rdata')
load('Processed/WB_obj.Rdata')

WB = WB[,-1]
sheet = sheet[-1,] #remove 'tube' sample

# Calculate sex based on data
GRset <- mapToGenome(WB)
Sex<-getSex(GRset, cutoff = -2)
GRsex = addSex(GRset)

sheet$predicted.sex = Sex$predictedSex
sheet$sexFlag = sheet$recorded.sex != sheet$predicted.sex

#' QC plots

MSet <- preprocessRaw(WB) 
qc <- getQC(MSet)

sheet$Sample_Name2 = sapply(lapply(strsplit(sheet$Sample_Name,'_'),'[',1:2), paste, collapse = '_')

#dropping remaining duplicate samples (1)
ord = order(sheet$sexFlag, sheet$Sample_Name2)
drop = duplicated(sheet$Sample_Name2[ord])[order(ord)]
drop[grepl('jurkat|jurkart',sheet$Sample_Name2,ignore.case = T)] = F

WB = WB[,!drop]
sheet = sheet[!drop,]

summaryqc <- shinySummarize(WB)

betas_raw = getBeta(WB)
save(betas_raw, file = 'Processed/betas_raw_dnamage.Rdata')

