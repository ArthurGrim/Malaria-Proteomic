library(MSnbase)
library(BiocParallel)

library(limma)
library(qvalue)


setwd("~/Bureau/Master/Semestre 9/BioStatisticsII/Final assignment")


#PSM report in a dataframe:
tpsms <- list()
tpsms[[1]] <- read.csv("quant-TMT6-1-clean.txt", sep = "\t",stringsAsFactors = F)
tpsms[[2]] <- read.csv("quant-TMT6-2-clean.txt", sep = "\t",stringsAsFactors = F)
psms <- rbind(tpsms[[1]],tpsms[[2]])

#Keep only high confidence and remove decoy hits.
psms <- psms[psms$Confidence.... == 100,]
psms$isDecoy <- grepl("_REVERSED",psms$Sequence)

#Prepare data frame for usage in MSnbase.
psms$rank <- 1
psms$desc <- psms$Protein.s.
psms$spectrumID <- psms$Spectrum.Title 
psms$spectrumFile <- psms$Spectrum.File
psms$idFile <- "Testfile"
psms$sequence <- psms$Modified.Sequence

#Reading spectra from .filtered files.
myExp <- Specs <- list()
mgf_file <- c("mgf/TMT6-1.mgf","mgf/TMT6-2.mgf")

#Reading mgf-files can take a while.
for (i in 1:2) {
  myExp[[i]] <- readMgfData(mgf_file[i], verbose=F)
  
}

#To be able to start from this point again, save MSExp objects (of even entire workspace).
save(myExp, file="ReadMgfs.RData")
load("ReadMgfs.RData")

# Match identification data to spectra, this also takes time
for(i in 1:length(myExp)) {
  
  Specs[[i]] <- addIdentificationData(myExp[[i]],psms,decoy="isDecoy",rank="rank",acc="Protein.s.",
                                      icol="Spectrum.Title",fcol="TITLE",desc="desc",
                                      pepseq="Modified.Sequence",verbose=T)
}


save(Specs, file="Specs.RData")
load("Specs.RData")

#Distribution of peptides per protein group
par(mfrow=c(1,2))
for (i in 1:length(Specs)) {
  
  Peps <- unique(fData(Specs[[i]])[,c("desc","sequence")])
  hist(table(Peps[,1]),200)
  
}
par(mfrow=c(1,1))


cl <- MulticoreParam(1) #For Linux and Mac
cl <- SnowParam(workers = 4, type = "SOCK") # For windows. It can work for Linux and Mac too by type = "FORK".

#Retrieving quantifications (reporter ion intensities) from spectra (this also takes time).

quant.method <- TMT6

qnt <- list()
for (i in 1:length(myExp)){
  
  qnt2 <- quantify(Specs[[2]], 
                       method = "sum", 
                       reporters = TMT6, 
                       strict = F, verbose = T, BPPARAM = cl)
  
}

                                    save(qnt2, file="Qnts2.RData")

                                    
load(file="Qnts.RData")
load(file="Qnts2.RData")
qnt[[2]] = qnt2

#Further necessary data operations.
prots <- list()
for (i in 1:length(qnt)) {
  
  tqnt <- qnt[[i]]
  #Reporters are impure. This needs to be correct by solving a linear equation system.
  imp<-makeImpuritiesMatrix(length(TMT6),edit=F)
  tqnt <- purityCorrect(tqnt,imp)
  #Transformation to log-scale
  exprs(tqnt) <- log2(exprs(tqnt))
  #Normalization by median (columns)
  exprs(tqnt) <- t(t(exprs(tqnt)) - apply(exprs(tqnt), 2, median, na.rm=T))   
  #Remove quantifications without identified peptide.
  tqnt <- removeNoId(tqnt)
  #Update feature names in preparation of merging.
  tqnt <- updateFeatureNames(tqnt,label = paste("Sample",i))

  #Data imputation
  #x <- impute(tqnt, "MLE")
  #processingData(x)
  #Potein summarization
  prots[[i]] <- combineFeatures(tqnt, group=strsplit(fData(tqnt)[["Protein.s."]], ", "), method="NTR",verbose=T)
  #Calculation of log-ratios versus mean to make different iTRAQ runs comparable.
  exprs(prots[[i]]) <- exprs(prots[[i]]) - rowMeans(exprs(prots[[i]]), na.rm=T)
  
}

fData(tqnt)

save(prots, file = "prots.RData")
load("prots.RData")





######################################################################################################################################

####################################################################################################################################


exprs(prots[[1]])


quant1Pshort = data.frame(exprs(prots[[1]]))
quant2Pshort = data.frame(exprs(prots[[2]]))


quant1Pshort[,1:6] <- sapply(quant1Pshort[,1:6],as.numeric)
quant2Pshort[,1:6] <- sapply(quant2Pshort[,1:6],as.numeric)

colnames(quant1Pshort) = c("NI-1","NI-2","d3-1","d3-2","ECM-1","ECM-2")
colnames(quant2Pshort) = c("ECM-1","ECM-2","NI-1","NI-2","d3-1","d3-2")

quant2Pshort = quant2Pshort[,c(3,4,5,6,1,2)] #order conditions as in quant1Pshort

#scaling 
quant2Pshort[,2:7] = scale(x = quant2Pshort[,2:7],center = quant2Pshort[quant2Pshort == "P02754",2:7])

quant1Pshort = quant1Pshort[complete.cases(quant1Pshort), ] #remove NA values
quant2Pshort = quant2Pshort[complete.cases(quant2Pshort), ] #remove NA values

#mergedExperiments

quantTot = merge.data.frame(quant1Pshort,quant2Pshort, by = 0)
row.names(quantTot) = quantTot[,1]
quantTot$Row.names = NULL
quantTot = quantTot[,c(1,2,7,8,3,4,9,10,5,6,11,12)]
colnames(quantTot) = c("NI-1","NI-2","NI-3","NI-4","d3-1","d3-2","d3-3","d3-4","ECM-1","ECM-2","ECM-3","ECM-4")
quantTot <- quantTot - rowMeans(quantTot, na.rm=T)

#remove extreme values

quant1Pshort = quant1Pshort[apply(quant1Pshort, 1, min) > -2.5,]
quant1Pshort = quant1Pshort[apply(quant1Pshort, 1, max) < 2.5,]


##################### show ratio in TMT6 channels #########################

  
par(mfrow=c(1,1))

plot(1, type="n", xlab="Conditions", ylab="Ratios", xlim=c(1, 6), ylim=c(-4,4 ),xaxt = "n")
axis(1, at=1:6, labels=c("NI-1","NI-2","d3-1","d3-2","ECM-1","ECM-2"))


for (i in 1:length(quant1Pshort[,1])){
  lines(t(quant1Pshort[i,1:6]))
}


#########################  PCA separate exp  ##############################

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

pca.out <- princomp(quantTot)
plot(pca.out$loadings, col = rep(1:3, each = 4), pch = 19, main = "Merged TMT6")
legend("topright", legend = c("NI","d3","ECM"), col = 1:3, pch = 19)


pca.out <- princomp(quant1Pshort)
plot(pca.out$loadings, col = rep(1:3, each = 2), pch = 19, main = "TMT6#1")
legend("topright", legend = c("NI","d3","ECM"), col = 1:3, pch = 19)

pca.out <- princomp(quant2Pshort)
plot(pca.out$loadings, col = rep(1:3, each = 2), pch = 19, main = "TMT6#2")
legend("topright", legend = c("NI","d3","ECM"), col = 1:3, pch = 19)




################### Statistical testing using ANOVA  #################

ExpDesign <- c(1,1,2,2,3,3)

pvals <- matrix(NA, nrow=nrow(quantToTest), ncol=3)

quantToTest = quant1Pshort

afp = quantToTest

for (i in 1:nrow(quantToTest)) {
  
  testIn <- data.frame(y=as.vector(t(quantToTest[i,1:6])), exp=factor(ExpDesign))
  lm.out <- lm(y~exp, testIn)
  pvals[i,1:length(summary(lm.out)$coefficients[,4])] <- summary(lm.out)$coefficients[,4]
  print(pvals[,2][i]  < 0.01 )
  
}

afp$pval.d3.NI = pvals[,2]
afp$pval.ECM.NI = pvals[,3]

ExpDesign <- c(1,1,2,2)
pvals <- matrix(NA, nrow=nrow(quantToTest), ncol=2)


for (i in 1:nrow(quantToTest)) {
  
  testIn <- data.frame(y=as.vector(t(quantToTest[i,3:6])), exp=factor(ExpDesign))
  lm.out <- lm(y~exp, testIn)
  pvals[i,1:length(summary(lm.out)$coefficients[,4])] <- summary(lm.out)$coefficients[,4]
  print(pvals[,2][i] < 0.01 )
  
}

afp$pval.ECM.d3 = pvals[,2]

afp<-afp[!(is.na(afp$pval.ECM.NI)),] #remove NA values


############### ANOVA ECM vs ni+d3 ##################

ExpDesign <- c(1,1,1,1,2,2)

pvals <- matrix(NA, nrow=nrow(quantToTest), ncol=3)

quantToTest = quant1Pshort

afp = quantToTest

for (i in 1:nrow(quantToTest)) {
  
  testIn <- data.frame(y=as.vector(t(quantToTest[i,1:6])), exp=factor(ExpDesign))
  lm.out <- lm(y~exp, testIn)
  pvals[i,1:length(summary(lm.out)$coefficients[,4])] <- summary(lm.out)$coefficients[,4]
  print(pvals[,2][i]  < 0.01 )
  
}

afp$pval.ECM.NI = pvals[,2]



#Statistical testing using Limma vs. the first condition.

quantToTest = quant1Pshort
afp = quantToTest
ExpDesign <- c(1,1,1,1,2,2)
fit <- lmFit(quantToTest, design=model.matrix(~ ExpDesign))
fit <- eBayes(fit)

#Extract all p-values.
pvalues.limma <- fit$p.value

qvals.limma <- list()

qvals.limma <-  qvalue(pvalues.limma[,2])$qvalues

print(sum(pvalues.limma < 0.01, na.rm = T))
print(sum(qvals.limma < 0.01, na.rm = T))
afp$pval.ECM.NI = qvals.limma

  

#visualization of p-vals:

par(mfrow=c(1,3))
hist(afp$pval.d3.NI, main = "NI vs d3", xlab = "pVal")
hist(afp$pval.ECM.d3, main = "ECM vs d3", xlab = "pVal")
hist(afp$pval.ECM.NI, main = "ECM vs NI", xlab = "pVal")


#To remove rows based on pvalues

afp<-afp[!(afp$pval.ECM.NI>=0.01),]


######## HEATMAP ########
afpHeat = data.frame( (afp[,1]+afp[,2] + afp[,3]+afp[,4])/4 , (afp[,5]+afp[,6])/2 )
rownames(afpHeat) = rownames(afp)
colnames(afpHeat) = c("NI, D3","ECM")

res = uniprot_mapping(rownames(afpHeat))
row.names(res) = res$Entry
afpHeat = merge(afpHeat,res,by=0)

rownames(afpHeat) = afpHeat$Protein.names
rownames(afpHeat) = gsub("\\(.*","",rownames(afpHeat))
rownames(afpHeat) = substr(rownames(afpHeat),1,35)

library(pheatmap)
library(RColorBrewer)
library(viridis)

breaksList = seq(-2, 2, by = 0.04)
pheatmap(afpHeat[,2:3], color = inferno(10))



#retrieve prot name:

uniprot_mapping <- function(ids) {
  uri <- 'http://www.uniprot.org/uniprot/?query='
  idStr <- paste(ids, collapse="+or+")
  format <- '&format=tab'
  fullUri <- paste0(uri,idStr,format)
  dat <- read.delim(fullUri)
  dat
}

res = uniprot_mapping(rownames(afpHeat))

#Clustering


library(factoextra)
set.seed(3)

fviz_nbclust(afp, FUNcluster = kmeans) #determine optimal n of clusters

numberClusters = 3

km.res <- kmeans(afp, numberClusters, nstart = 100)

afp$clust = km.res$cluster

par(mfrow=c(2,numberClusters))

#Visualizing clustering results

for (i in 1:numberClusters){
  afpClust = afp[afp$clust==i,]
  
  plot(1, type="n", xlab="Conditions", ylab="Ratios", xlim=c(1, 6), ylim=c(-5, 5),xaxt = "n")
  axis(1, at=1:6, labels=c("NI-1","NI-2","d3-1","d3-2","ECM-1","ECM-2"))
  
  
  for (i in 1:nrow(afpClust)){
    lines(t(afpClust[i,1:6]))
  }
}


for (i in 1:numberClusters){
  afpClust = afp[afp$clust==i,]
  boxplot(afpClust$pval.ECM.NI,ylim=c(0,0.01),xlab = "p-Value ECM vs NI") 
  abline(h = 0.01, col = "red")
}


### Create Latex Table ###

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


grp = unique(unlist(strsplit(row.names(afp), ", ")))

grp = strsplit(row.names(afp), ", ")

ResultHH = row.names(afp[afp$clust==3 , ]) #extract over expressed (highly)
ResultH = row.names(afp[afp$clust==2 , ]) #extract over expressed
ResultL = row.names(afp[afp$clust==1 , ])


ResultHH = unique(unlist(strsplit(ResultHH, ", ")))
ResultH = unique(unlist(strsplit(ResultH, ", ")))
ResultL = unique(unlist(strsplit(ResultL, ", ")))


colnames(TableResultH) = c("Protein","ECM/NI ratio","pValue ECM-NI")
colnames(TableResultL) = c("Protein","ECM/NI ratio","pValue ECM-NI")

sq <- seq(max(length(ResultHH), length(ResultH)))
data.frame(ResultHH[sq], ResultH[sq])


xtable(data.frame(ResultHH,ResultH,ResultL))




############   MISC ##################

cor(quant2Pshort$`ECM-2`, quant2Pshort$`ECM-1`)
cor(quant1Pshort$`ECM-2`, quant1Pshort$`ECM-1`)

cor(quantTot$`ECM-2`,quantTot$`ECM-3`)

  cor(rowMeans(quant1Pshort[,1:2]), rowMeans(quant1Pshort[,3:4]))


