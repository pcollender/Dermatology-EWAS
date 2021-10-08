#' List of idat files
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
WB <- read.metharray.exp(targets = sheet, recursive = T, verbose = T)

dir.create('Processed')
save(WB, file = 'Processed/WB_obj.Rdata')

# Calculate sex based on data
GRset <- mapToGenome(WB)
Sex<-getSex(GRset, cutoff = -2)
GRsex = addSex(GRset)

write.csv(as.data.frame(Sex), file = file.path('Processed', 'Sex_Estimate_Minfi.csv'), row.names = T)  
png(file.path('Processed','Sex_estimate.png'))
plotSex(GRsex)
dev.off()

sheet$predicted.sex = Sex$predictedSex
sheet$sexFlag = sheet$recorded.sex != sheet$predicted.sex

#' QC plots

MSet <- preprocessRaw(WB) 
qc <- getQC(MSet)

png(file.path('Processed','qcplot.png'))
plotQC(qc)
dev.off()

include = (qc$mMed + qc$uMed) / 2 > 10.5
sheet$Sample_Name2 = sapply(lapply(strsplit(sheet$Sample_Name,'_'),'[',1:2), paste, collapse = '_')

WB = WB[,include]
sheet = sheet[include,]

#dropping remaining duplicate samples (1)
ord = order(sheet$sexFlag, sheet$Sample_Name2)
drop = duplicated(sheet$Sample_Name2[ord])[order(ord)]
drop[grepl('jurkat|jurkart',sheet$Sample_Name2,ignore.case = T)] = F

WB = WB[,!drop]
sheet = sheet[!drop,]

summaryqc <- shinySummarize(WB)

#' Step #2

setwd('Processed')
qcReport(WB, sampNames=sheet$Sample_Name2, sampGroups = sheet$Slide, pdf =  "qcReport.pdf")
setwd('../')

QC_minfi <- minfiQC(GRset, fixOutliers = TRUE)

#' Densities
table(pData(WB)$Sample_Plate)
table(pData(WB)$Slide)
table(pData(WB)$Array)
#' by Slide
png('Processed/Slide_densities.png')
densityPlot(MSet, sampGroups = pData(WB)$Slide)
dev.off()
#' by Row
png('Processed/Array_densities.png')
densityPlot(MSet, sampGroups = pData(WB)$Array)
dev.off()

#' Preprocess Funnorm
WB.FunNorm <- preprocessFunnorm(WB)
#Warning message:
#In .getSex(CN = CN, xIndex = xIndex, yIndex = yIndex, cutoff = cutoff) :
#  An inconsistency was encountered while determining sex. One possibility is that only one sex is present. We recommend further checks, for example with the plotSex function.

save(WB.FunNorm, file = 'Processed/WB_FunNorm.Rdata')

#' Distribution of betas before and After normalization

png('Processed/FunNorm_Density.png')
densityPlot(WB, main = "density plots before and after preprocessing", pal="#440154FF", ylim=c(0,4.5))
densityPlot(getBeta(WB.FunNorm), add = F, pal = "#FDE725FF")
# Add legend
legend("topleft", c("FunNorm","Raw"), 
       lty=c(1,1), title="Normalization", 
       bty='n', cex=1.3, col=c("#FDE725FF","#440154FF"))

dev.off()

#' Detection pvalues
detect.p <- detectionP(WB, type = "m+u")
save(detect.p, file = 'Processed/detection_p.Rdata')

png('Processed/mean_detection_p.png')
barplot(colMeans(detect.p), col=rainbow(dim(detect.p)[2]), las=2,# ylim=c(0,7e-4),
        cex.names=0.7, ylab="Mean detection P",cex.axis=0.9)
dev.off()

#' Detection by sample
table(t(as.matrix(colSums(detect.p > 0.01))))

png('Processed/n_nondetect.png')
barplot(colSums(detect.p > 0.01), col=rainbow(dim(detect.p)[2]), las=2, 
        cex.names=0.7, ylab="Number of non-detects per sample",cex.axis=0.9)
dev.off()
#'
detect.p[detect.p > 0.01] <- NA
detect.p <- na.omit(detect.p) #dropping loci with any non-detects
intersect <- intersect(rownames(getAnnotation(WB)), rownames(detect.p))
length(intersect)

nrow(WB.FunNorm)
## [1] 865859
WB.FunNorm <- WB.FunNorm[rownames(getAnnotation(WB.FunNorm)) %in% intersect,]
nrow(WB.FunNorm)
## [1] 858344: 7515 loci dropped if dropping those with any non-detects

#########################################################################
##
##            Probe-type
##
#' Preprocess Quantile
betas.FunNorm<-getBeta(WB.FunNorm)
range(betas.FunNorm)
#[1]  0.0000000 0.9952753
median(betas.FunNorm)
#[1] 0.6516782

## Get values from FunNorm Data
mindetval = min(betas.FunNorm[betas.FunNorm != 0])
betas.FunNorm_nz = replace(betas.FunNorm, which(betas.FunNorm==0), mindetval/2)
raw.M <- logit2(betas.FunNorm_nz)
range(raw.M)
#[1]      -13.522302   7.718724
## Other objects needed
dist=25
quantile.grid=seq(0.001,0.999,by=0.001)
qcscore=NULL
nbthre=3
detPthre=0.000001

# find therby pairs of type I probes and type II probes;; harmonizing type II probes to type I
#PAC: what is a therby pair?
annotation<-getAnnotation(GRset)
annotation=annotation[intersect(rownames(betas.FunNorm),rownames(annotation)),]
dim(annotation);dim(betas.FunNorm)

probe.II.Name=annotation$Name[annotation$Type=="II"]
annotation=annotation[order(annotation$chr,annotation$pos),]
anno1=annotation[1:(nrow(annotation)-1),]
anno2=annotation[2:nrow(annotation),]
flag=(abs(anno1$pos-anno2$pos)<dist & anno1$chr==anno2$chr & 
        anno1$Relation_to_Island==anno2$Relation_to_Island & anno1$Type !=
        anno2$Type)
anno1=anno1[flag,]
anno2=anno2[flag,]
probe.I=anno1$Name
probe.II=anno2$Name
probe.I[anno2$Type=="I"]=anno2$Name[anno2$Type=="I"]
probe.II[anno1$Type=="II"]=anno1$Name[anno1$Type=="II"]

raw.M.t=raw.M[c(probe.I,probe.II),]

#remove low quality data
if(is.null(qcscore)){}else if((sum(!(rownames(raw.M.t) %in% 
                                     rownames(qcscore$detP))) +
                               sum(!(colnames(raw.M.t) %in% colnames(qcscore$detP))))>0){
  stop("Wrong qcscore matrix, please check...\n")}else{
    temp <- qcscore$nbead<nbthre | qcscore$detP>detPthre
    temp=temp[rownames(raw.M.t),]
    temp=temp[,colnames(raw.M.t)]
    raw.M.t[temp]=NA
  }

#linear regression
M.II<-raw.M.t[probe.II,]
M.I<-raw.M.t[probe.I,]


##quantile.grid=seq(0.001,0.999,by=0.001)
qtl<-function(x) quantile(x, quantile.grid, na.rm=TRUE)
M.I=apply(M.I,2,qtl)
M.II=apply(M.II,2,qtl)

beta.est<-mat.or.vec(2,ncol(betas.FunNorm))

for (i in 1:ncol(betas.FunNorm)){
  index<-(M.II[,i]!=Inf & M.II[,i]!=-Inf & M.I[,i]!=Inf & M.I[,i]!=-Inf)
  X<-cbind(rep(1,sum(index)),M.II[index,i]); Y<-M.I[index,i]
  beta.est[,i]<-solve(t(X)%*%X)%*%t(X)%*%Y
}

M.II.all<-raw.M[probe.II.Name,]
M.II.new<-mat.or.vec(nrow(M.II.all),ncol(M.II.all))
for (i in 1:ncol(M.II.all)){
  M.II.new[,i]<-beta.est[1,i]+beta.est[2,i]*M.II.all[,i]
}
M.II.new[M.II.all==Inf]<-Inf; M.II.new[M.II.all==-Inf]<-(-Inf)

betas.FunNorm[probe.II.Name,]<-ilogit2(M.II.new)
beta.rcp<-betas.FunNorm
range(beta.rcp)
sum(is.na(beta.rcp))


## Remove non-detects
beta.rcp <- beta.rcp[rownames(beta.rcp) %in% intersect,]
sheet$basename = sapply(strsplit(sheet$Basename,'/'),'[',7)
sheet<-sheet[order(sheet$basename),]
beta.rcp<-beta.rcp[,sort(colnames(beta.rcp))]
colnames(beta.rcp)==sheet$basename
identical(colnames(beta.rcp),sheet$basename)

typeI <-   minfi::getProbeInfo(WB,type="I")$Name
typeII <-  minfi::getProbeInfo(WB,type="II")$Name
onetwo <- rep(1, nrow(beta.rcp))
onetwo[rownames(beta.rcp) %in% typeII] <- 2
# almost 84% of our probes are type II: PAC - true
knitr::kable(t(table(onetwo)))

#" Figure

png('Processed/rcp_adjustment_plot.png')
par(mfrow=c(1,2)) # Side-by-side density distributions 
sampGroups = rownames(getAnnotation(WB.FunNorm)) %in% typeI
sampGroups = c('Type I','Type II')[as.numeric(!sampGroups) + 1]
densityPlot(getBeta(WB.FunNorm), sampGroups = sampGroups,pal = c("#FDE725FF","#440154FF"),
            main='Beta density',ylim=c(0,6.5))

densityPlot(beta.rcp, sampGroups = sampGroups,pal = c("#FDE725FF","#440154FF"),
            main='Beta density probe-type adjusted',ylim=c(0,6.5))

dev.off()

save(beta.rcp, sheet, file = 'Processed/beta.rcp.RData')

#' Chip issue

qual_col_pals  <- brewer.pal.info[brewer.pal.info$category == 'qual',]
pal <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

par(mfrow=c(1,1)) 

png('Processed/MDS_plots.png')
plotMDS(beta.rcp, top=10000, gene.selection="common",
        pch=24,  #labels = sheet$Slide,
        bg=alpha(pal[factor(sheet$Array)],0.75),
        col="black", dim=c(2,3),cex=3)
dev.off()

#'Combat by Chip ;;; Adjusting for chip / batch-wide differences in fluorescence
#'Maybe some residual differences by chip / by lane
range(beta.rcp)
#replacing zero values with minimum nonzero/2
bmin = min(beta.rcp[beta.rcp!=0])
beta.rcp.nz = replace(beta.rcp,beta.rcp==0,bmin/2)
Mvals <- log2(beta.rcp.nz)-log2(1-beta.rcp.nz)


Mvals.ComBat <- ComBat(Mvals, batch = sheet$Sample_Plate)

betas.rcp_m.adj <- 2^Mvals.ComBat/(1+2^Mvals.ComBat)

colnames(betas.rcp_m.adj)==sheet$basename
identical(colnames(betas.rcp_m.adj),sheet$basename)

#Drop sample in isolated group 'Tube', otherwise ComBat only adjusts for mean differences and no rescaling
Mvals.ComBat <- ComBat(Mvals[,-1], batch = sheet$Sample_Plate[-1])

betas.rcp_full.adj <- 2^Mvals.ComBat/(1+2^Mvals.ComBat)
sheet_full.adj = sheet[-1,]

colnames(betas.rcp_full.adj)==sheet_full.adj$basename
identical(colnames(betas.rcp_full.adj),sheet_full.adj$basename)

save(betas.rcp_m.adj, 
     betas.rcp_full.adj,
     sheet, sheet_full.adj, 
     file = 'Processed/beta.rcp.ComBat.RData')
