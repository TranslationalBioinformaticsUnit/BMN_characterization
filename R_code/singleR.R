library(Seurat)
library(scales)
library(biomaRt)
library(dplyr)
library(SingleR)
library(ggplot2)
library(patchwork)


human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <-  useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#### mouse MSC 
load(file="/Data/MSC/MSC.integration_sub_20200624.RData") 
load(file = "/Data/20201009_hBMN.MSC.sub.RData") 

MSC.markers <- read.csv("/Markers/MSCs/Final_12_clusters_MSC.csv", row.names = 1)


hBMN.MSC.sub@active.assay ###"RNA"
hBMN.MSC.sub <- FindVariableFeatures(hBMN.MSC.sub, nfeatures = 3000)
genes <- VariableFeatures(hBMN.MSC.sub)
hBMN.MSC.sub@active.assay <- "integrated"
hBMN.MSC.sub<- ScaleData(hBMN.MSC.sub,features = genes)

dim(hBMN.MSC.sub@assays$integrated@scale.data)
hMSC.genes_conveted <- getLDS(attributes = "hgnc_symbol",filters = "hgnc_symbol",values =rownames(hBMN.MSC.sub@assays$integrated@scale.data), mart = human, attributesL = "mgi_symbol", martL = mouse)
hMSC.genes_conveted.sub <- hMSC.genes_conveted[hMSC.genes_conveted$MGI.symbol%in%rownames(MSC.integration_sub@assays$integrated@data),]

hMSC <- GetAssayData(hBMN.MSC.sub,assay = "integrated",slot = "scale.data")
hMSC.sub <- hMSC[rownames(hMSC)%in%unique(hMSC.genes_conveted.sub$HGNC.symbol),]
hMSC.sub <- data.frame(hMSC.sub)
hMSC.sub.names <- hMSC.genes_conveted.sub$MGI.symbol[match(rownames(hMSC.sub),hMSC.genes_conveted.sub$HGNC.symbol)]
which(duplicated(hMSC.sub.names))
hMSC.sub <- hMSC.sub[-which(duplicated(hMSC.sub.names)),]
hMSC.sub.names<- hMSC.sub.names[-which(duplicated(hMSC.sub.names))]
rownames(hMSC.sub) <- hMSC.sub.names

mMSC <- GetAssayData(MSC.integration_sub, assay = "integrated", slot = "scale.data")
pred.grun.MSC <- SingleR(test=hMSC.sub, ref=mMSC, labels=MSC.integration_sub$bootstrapped, de.method="wilcox")
table(pred.grun.MSC$labels)

plotScoreHeatmap(pred.grun.MSC)
plotDeltaDistribution(pred.grun.MSC, ncol = 3)
plab.MSC <- data.frame(pred.grun.MSC$pruned.labels)
rownames(plab.MSC) <- rownames(pred.grun.MSC)
#plab.EC$pred.grun.pruned.labels[is.na(plab.EC$pred.grun.EC.pruned.labels)] <- "Not Assigned"
rownames(plab.MSC) <- gsub("\\.", "-",rownames(plab.MSC))
hBMN.MSC.sub<- AddMetaData(hBMN.MSC.sub,metadata = plab.MSC)
DimPlot(hBMN.MSC.sub, group.by = "pred.grun.MSC.pruned.labels")


#### mouse EC
###load(/Data/20210228_hEC.singleR.obj.for.plots.RData")
load(file ="/Data/EC.integration_614.RData") 
hBMN.EC.rmC4 <-readRDS("/Data/hBMN.EC.rmC4_20201026.Rds")
EC.markers <- read.csv("/Markers/ECs/all_ECs_14_clusters.csv", row.names = 1)

hBMN.EC.rmC4@active.assay ###"RNA"
hBMN.EC.rmC4 <- FindVariableFeatures(hBMN.EC.rmC4, nfeatures = 3000)
genes <- VariableFeatures(hBMN.EC.rmC4)
hBMN.EC.rmC4@active.assay <- "integrated"
hBMN.EC.rmC4 <- ScaleData(hBMN.EC.rmC4,features = genes)


dim(hBMN.EC.rmC4@assays$integrated@scale.data)
hEC.genes_conveted <- getLDS(attributes = "hgnc_symbol",filters = "hgnc_symbol",values =rownames(hBMN.EC.rmC4@assays$integrated@scale.data), mart = human, attributesL = "mgi_symbol", martL = mouse)
hEC.genes_conveted.sub <- hEC.genes_conveted[hEC.genes_conveted$MGI.symbol%in%rownames(EC.integration@assays$integrated@data),]

hEC <- GetAssayData(hBMN.EC.rmC4,assay = "integrated",slot = "scale.data")
hEC.sub <- hEC[rownames(hEC)%in%unique(hEC.genes_conveted.sub$HGNC.symbol),]
hEC.sub <- data.frame(hEC.sub)
hEC.sub.names <- hEC.genes_conveted.sub$MGI.symbol[match(rownames(hEC.sub),hEC.genes_conveted.sub$HGNC.symbol)]
which(duplicated(hEC.sub.names))
hEC.sub <- hEC.sub[-which(duplicated(hEC.sub.names)),]
hEC.sub.names<- hEC.sub.names[-which(duplicated(hEC.sub.names))]
rownames(hEC.sub) <- hEC.sub.names

mEC <- GetAssayData(EC.integration, assay = "integrated", slot = "scale.data")
pred.grun.EC <- SingleR(test=hEC.sub, ref=mEC, labels=EC.integration$bootstrapped, de.method="wilcox")
table(pred.grun.EC$labels)

plotScoreHeatmap(pred.grun.EC)
plotDeltaDistribution(pred.grun.EC, ncol = 3)
plab.EC <- data.frame(pred.grun.EC$pruned.labels)
rownames(plab.EC) <- rownames(pred.grun.EC)
#plab.EC$pred.grun.pruned.labels[is.na(plab.EC$pred.grun.EC.pruned.labels)] <- "Not Assigned"
rownames(plab.EC) <- gsub("\\.", "-",rownames(plab.EC))
hBMN.EC.rmC4 <- AddMetaData(hBMN.EC.rmC4,metadata = plab.EC)
DimPlot(hBMN.EC.rmC4, group.by = "pred.grun.EC.pruned.labels")

load(file = "/Data/hBMN_4_sample_filtered_IG_integrated.RData")
hBMN.EC <- subset(hBMN.integrated, cells = colnames(hBMN.EC.rmC4))

hBMN.EC <- AddMetaData(hBMN.EC,metadata = plab.EC)
DimPlot(hBMN.EC, group.by = "pred.grun.EC.pruned.labels")


out.dir <- "/Figures/"
prefix <- "20210603_Fig6_SingleR_EC"
file_name <- paste0(out.dir,prefix,"1.pdf")
pdf(file = file_name, width = 8, height = 8)
plotScoreHeatmap(pred.grun.EC)
dev.off()

file_name <- paste0(out.dir,prefix,"2.pdf")
pdf(file = file_name, width = 8, height = 8)
plotDeltaDistribution(pred.grun.EC, ncol = 3)
dev.off()

file_name <- paste0(out.dir,prefix,"3.pdf")
pdf(file = file_name, width = 8, height = 8)
DimPlot(hBMN.EC, group.by = "pred.grun.EC.pruned.labels")
dev.off()


prefix <- "20210603_Fig6_SingleR_MSC"
file_name <- paste0(out.dir,prefix,"1.pdf")
pdf(file = file_name, width = 8, height = 8)
plotScoreHeatmap(pred.grun.MSC)
dev.off()

file_name <- paste0(out.dir,prefix,"2.pdf")
pdf(file = file_name, width = 8, height = 8)
plotDeltaDistribution(pred.grun.MSC, ncol = 3)
dev.off()

file_name <- paste0(out.dir,prefix,"3.pdf")
pdf(file = file_name, width = 8, height = 8)
DimPlot(hBMN.MSC.sub, group.by = "pred.grun.MSC.pruned.labels")
dev.off()

### check the markers shared between mouse and human MVG
MSC.markers.sub <- filter(MSC.markers, cluster!="MSC_C7")
mMSC.sub <- mMSC[rownames(mMSC)%in%MSC.markers.sub.100$gene,]

Human_mouse_shared_genes <- function(hEC.sub, EC.markers.100.ranked) {
  hEC_genes_in_mosue <- data.frame(gene = character(),cluster=character(),rank=numeric(), avg_logFC = numeric())
  for (i in unique(EC.markers.100.ranked$cluster)){
    EC.markers.cluster<- filter(EC.markers.100.ranked,cluster==i)
    print(i)
    print(length(intersect(rownames(hEC.sub),EC.markers.cluster$gene)))
    print(intersect(rownames(hEC.sub),EC.markers.cluster$gene))
    shared_genes <- intersect(rownames(hEC.sub),EC.markers.cluster$gene)
    shared_gene_matrix <- EC.markers.cluster[which(EC.markers.cluster$gene%in%shared_genes),c("gene","cluster","rank","avg_logFC")]
    shared_gene_matrix <- shared_gene_matrix[order(shared_gene_matrix$rank),]
    hEC_genes_in_mosue <- rbind(hEC_genes_in_mosue, shared_gene_matrix)
  }
  return(hEC_genes_in_mosue)
}

MSC.markers.ranked <- MSC.markers.sub%>% group_by(cluster) %>% mutate(rank=rank(desc(avg_logFC),cluster))
hMSC_genes_in_mouse <- Human_mouse_shared_genes(hEC.sub = hMSC.sub,EC.markers.100.ranked = MSC.markers.ranked)

View(hMSC_genes_in_mouse)
table(hMSC_genes_in_mouse$cluster)

length(intersect(rownames(hMSC.sub), unique(MSC.markers.ranked$gene)))



