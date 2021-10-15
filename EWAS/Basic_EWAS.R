#Preliminaries
library(limma)
library(qvalue)
library(readxl)
setwd('Dermatology')

# Load Phenotype

phenotype = read_xlsx('Ages.xlsx')

phenotype = subset(phenotype, !is.na(TYPE))
phenotype$condition = as.numeric(grepl('hidradenitis', phenotype$TYPE))

ids = phenotype$`Material ID`

phenotype$age = sapply(strsplit(gsub(' m','',phenotype$AGE),' y, '), function(x) as.numeric(x[1]) + as.numeric(x[2])/12)

#Load methylation data
load('Processed/beta.rcp.ComBat.RData')
sheet = sheet_full.adj
betas = betas.rcp_full.adj

colnames(betas) = sheet$Sample_Name2[match(colnames(betas), sheet$basename)]
#Add observed sex to phenotype (estimated via minfi)
phenotype$sex = sheet$predicted.sex[match(ids, sheet$Sample_Name2)]

#align betas and cell types with phenotype data
betas = betas[,colnames(betas) %in% ids]
betas = betas[,match(ids,colnames(betas))]

#Categorical analysis
Z = model.matrix(~ age + sex + condition, data = phenotype)

getinds = grep('condition',colnames(Z))

fit = lmFit(log(betas/(1-betas)), design = Z)

res = eBayes(fit, robust = T)

results<-data.frame(pval = res$p.value[,'condition'],
                    coef = res$coef[,'condition'],
                    se = res$coef[,'condition']/res$t[,'condition'])

results$CpG = rownames(betas)
#' Save Results from EWAS
dir.create('EWAS')
save(results,file=file.path('EWAS',"basic_EWAS_limma.RData"))

#get predictions of difference for modal / average values of covariates
contrmat = matrix(0, nrow = ncol(Z), ncol = 2)
colnames(contrmat) = c('normal','HS')  
rownames(contrmat) = colnames(Z)

contrmat[1:3,] = colMeans(Z)[1:3]
contrmat['sexM',] = round(contrmat['sexM',]) #most subjects are male (barely)
contrmat['condition',2] = 1

#reference comparison is for 28.814 year old male
contrfit2 = contrasts.fit(fit, contrmat)

res2 = eBayes(contrfit2, robust = T)

invlogit = function(x) 1/(1+exp(-x))

library(MASS)
library(parallel)
library(WGCNA)
library(pbmcapply)

#simulating coefficients from multivariate normal distribution
results2<-data.frame(absVal.ML = invlogit(res2$coef),
                     absVal.low = invlogit(res2$coef - 1.96 * res2$stdev.unscaled * res2$sigma),
                     absVal.high = invlogit(res2$coef + 1.96 * res2$stdev.unscaled * res2$sigma),
                     absChange.ML = invlogit(res2$coef[,-1]) - invlogit(res2$coef[,1]))

CIfun = function(i){
  variates = mvrnorm(3e3, res2$coef[i,], res2$cov.coefficients * res2$sigma[i]^2)
  contrasts = invlogit(variates[,-1]) - invlogit(variates[,1])
  low = colQuantileC(contrasts, 0.025)
  high = colQuantileC(contrasts, 0.975)
  names(low) = names(high) = colnames(contrasts)
  return(list(low = low,high = high))
}

CIdat = pbmclapply(1:nrow(res2$coefficients), CIfun, mc.preschedule = T) #parallel apply with progress bar
temp = data.frame(absChange.l95 = do.call(rbind, lapply(CIdat,'[[','low')),
                  absChange.u95 = do.call(rbind, lapply(CIdat,'[[','high')))
results2 = cbind(results2, temp)

save(results2,file=file.path('EWAS',"basic_EWAS_limma_absdiff.RData"))
