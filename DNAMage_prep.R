#Prepping data for DNAMage analysis
library(tidyverse)
library(readxl)
load('Processed/beta.rcp.ComBat.RData')

colnames(betas.rcp_full.adj) = sheet_full.adj$Sample_Name2[match(colnames(betas.rcp_full.adj), sheet_full.adj$basename)]
betas.rcp_full.adj = cbind(data.frame(ProbeID = rownames(betas.rcp_full.adj)),betas.rcp_full.adj)

dir.create('Processed/DnaMAge')

#restrict to necessary CpGs
inc_cpgs = read.csv('Processed/DnaMAge/datMiniAnnotation3.csv')

#some cpgs missing from epic data, add in with NA values
missed_cpgs = setdiff(inc_cpgs$Name, betas.rcp_full.adj$ProbeID)

add.data = betas.rcp_full.adj[1:length(missed_cpgs), ]

add.data$ProbeID = missed_cpgs

add.data[,-1] = NA

betas.rcp_full.adj.use = betas.rcp_full.adj[betas.rcp_full.adj$ProbeID %in% inc_cpgs$Name,]
betas.rcp_full.adj.use = rbind(betas.rcp_full.adj.use, add.data)
betas.rcp_full.adj.use = betas.rcp_full.adj.use[match(inc_cpgs$Name, betas.rcp_full.adj.use$ProbeID),]
#Note - missing ~ 2575 probes... will this cause a problem?

write.csv(betas.rcp_full.adj.use, file = 'Processed/DnaMAge/betas.csv', row.names = F)

agedat = read_xlsx('Ages.xlsx') 

annot = data.frame(Age = agedat$AGE[match(colnames(betas.rcp_full.adj.use)[-1], agedat$`Material ID`)],
                   Female = as.numeric(sheet_full.adj$predicted.sex == 'F')[match(colnames(betas.rcp_full.adj.use)[-1], sheet_full.adj$Sample_Name2)],
                   Tissue = "Epidermis")

ages = strsplit(gsub(' m', '',annot$Age), ' y, ')

ages = sapply(ages, function(x) as.numeric(x[1]) + as.numeric(x[2])/12)

annot$Age = ages

write.csv(annot, file = 'Processed/DnaMAge/annot.csv', row.names = F)

#Using raw beta values, no preprocessing at all

load('Processed/betas_raw_dnamage.Rdata')
colnames(betas_raw) = sheet_full.adj$Sample_Name2[match(colnames(betas_raw), sheet_full.adj$basename)]

betas_raw = betas_raw[,!is.na(colnames(betas_raw))]  #forgot to drop a duplicate

betas_raw = cbind(data.frame(ProbeID = rownames(betas_raw)),betas_raw)

missed_cpgs = setdiff(inc_cpgs$Name, betas_raw$ProbeID)

add.data = betas_raw[1:length(missed_cpgs), ]

add.data$ProbeID = missed_cpgs

add.data[,-1] = NA

betas_raw = betas_raw[betas_raw$ProbeID %in% inc_cpgs$Name,]
betas_raw = rbind(betas_raw, add.data)
betas_raw = betas_raw[match(inc_cpgs$Name, betas_raw$ProbeID),]
#Note - missing ~ 2575 probes... will this cause a problem?

write.csv(betas_raw, file = 'Processed/DnaMAge/betas_raw.csv', row.names = F)
