library(BiocSingular)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(GO.db)
library(limma)
library(org.Hs.eg.db)
library(scDblFinder)
library(DropletUtils)
library(DropletTestFiles)
library(uwot)
library(rtracklayer)
library(PCAtools)
library(celldex)
library(SingleR)
library(batchelor)
library(bluster)

# No of cells filtered after QC
original: 2136
gene expression + antibody: 227
doublet: 47
final: 1862

# read data into sce object
fw.d0.sce <- read10xCounts("data/")

# split gene expression antibody libraries
fw.d0.sce <- splitAltExps(fw.d0.sce, rowData(fw.d0.sce)$Type)

# QC

unfiltered.fw.d0.sce<-fw.d0.sce # save a copy

mito <- grep("^MT-", rowData(fw.d0.sce)$Symbol)
  
ge.qc <- perCellQCMetrics(fw.d0.sce, subsets=list(Mito=mito))

reasons <- perCellQCFilters(ge.qc, sub.fields=c("subsets_Mito_percent"))

discard.ge <- reasons$discard

> table(discard.ge)
discard.ge
FALSE  TRUE 
1949   187 

colData(fw.d0.sce) <- cbind(colData(fw.d0.sce), ge.qc)

fw.d0.sce$discard.ge <- discard.ge

# look at discarded cells

  # total UMI  
  plotColData(fw.d0.sce, y="sum", colour_by="discard.ge")+ 
    scale_y_log10() + 
    ggtitle("Total count")    
# detected genes
  plotColData(fw.d0.sce, y="detected", colour_by="discard.ge")+ 
    scale_y_log10() + 
    ggtitle("Detected features")
 # mito percent
  plotColData(fw.d0.sce, y="subsets_Mito_percent", colour_by="discard.ge")+ 
  scale_y_log10() + 
    ggtitle("Mito percent")   
  
  plotColData(fw.d0.sce, x="sum", y="subsets_Mito_percent", colour_by="discard.ge")+
    scale_x_log10()
  
  plotColData(fw.d0.sce, x="detected", y="subsets_Mito_percent", colour_by="discard.ge")+
    scale_x_log10()

# QC for antibodies
rowData(altExp(fw.d0.sce))
  controls <- grep("^T", rownames(altExp(fw.d0.sce)))
  qc.stats <- cleanTagCounts(altExp(fw.d0.sce), controls=controls)
  print(summary(qc.stats$zero.ambient))
  colData(altExp(fw.d0.sce)) <- cbind(colData(altExp(fw.d0.sce)), qc.stats)
  discard.adt <- qc.stats$zero.ambient # use zero ambient as the only criteria
  altExp(fw.d0.sce)$discard.adt <- discard.adt
  
  > table(discard.adt)
  discard.adt
  FALSE  TRUE 
  2073    63 
  
  plotColData(altExp(fw.d0.sce), y="sum.controls", colour_by="discard.adt")+ 
    scale_y_log10() + 
    ggtitle("Total count")
  
  discard.all <- discard.adt | discard.ge
  > table(discard.all)
  discard.all
  FALSE  TRUE 
  1909   227 

fw.d0.sce<- fw.d0.sce[,!discard.all]
dim(fw.d0.sce)

# Normalization for Gene expression
set.seed(100)
clust <- quickCluster(fw.d0.sce)
fw.d0.sce <- computeSumFactors(fw.d0.sce, cluster=clust)
fw.d0.sce <- logNormCounts(fw.d0.sce)
> assayNames(fw.d0.sce)
[1] "counts"    "logcounts"RC

######################################## need to be re-done after doublet removal
# Finding variable features

dec.d0<-modelGeneVar(fw.d0.sce)

fit <- metadata(dec.d0)

plot(fit$mean, fit$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

rownames(fw.d0.sce)<-rowData(fw.d0.sce)$Symbol

hvgs.d0<-getTopHVGs(fw.d0.sce,prop=0.2)

# PCA

set.seed(100) 
fw.d0.sce<- fixedPCA(fw.d0.sce, rank=100, subset.row=hvgs.d0,name="PCA.1")

# nearest neighbor clustering
nn.clust.1 <- clusterCells(fw.d0.sce, use.dimred="PCA.1",BLUSPARAM=NNGraphParam(k=20,cluster.fun="louvain",type="jaccard"))
  
colLabels(fw.d0.sce) <- nn.clust.1

# UMAP
set.seed(100)
fw.d0.sce <- runUMAP(fw.d0.sce, dimred="PCA.1",n_neighbors=20,name="UMAP.1")

colLabels(fw.d0.sce) <- nn.clust.1
plotReducedDim(fw.d0.sce, dimred="UMAP.1",colour_by="label")

### Doublet removal
library(scDblFinder)
library(BiocSingular)

set.seed(100)
fw.d0.sce<-scDblFinder(fw.d0.sce, clusters=colLabels(fw.d0.sce),dbr=0.02,dbr.sd=0.01,
                          artificialDoublets = 4000)
head(colData(fw.d0.sce))

plotReducedDim(fw.d0.sce, dimred="UMAP.1",colour_by="scDblFinder.score")

# remove doublet
fw.d0.sce<-fw.d0.sce[,fw.d0.sce$scDblFinder.class=="singlet"]
dim(fw.d0.sce)
1909 -> 1862

### ADT export to map to FACS phenotypes
### normalize and log10 transformation
logtrans<-function(x){log10(x+1)}

# Day0
sz.f<-colData(altExp(fw.d0.sce))$sizeFactor
adt.norm <- normalizeCounts(altExp(fw.d0.sce), size.factors =sz.f, 
                            assay.type="counts", log=FALSE,center.size.factors=TRUE)
adt.norm.m<-t(adt.norm)
adt.norm.m<-as.matrix(adt.norm.m)
adt.lg10.d0<-apply(adt.norm.m, 2, logtrans)

save(adt.lg10.d0, file="~/sc-citeseq_data/Output/adt.lg10.d0.RData")
save(fw.d0.sce,ge.qc,qc.stats,file="~/sc-citeseq_data/Output/sce.d0.RData")
