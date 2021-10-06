library(dplyr)
library(Seurat)
library(ggplot2)
set.seed(666)

#### loading the data and filtering####3
filelist <- list.files("/Data/GSE128423_RAW/std")
filelist

d10x.data <- sapply(filelist, function(i){
  cat(i,": ", list.files(paste0("Data/GSE128423_RAW/std/",i)),"\n")
  d10x <- Read10X(data.dir=(paste0("Data/GSE128423_RAW/std/",i)))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

bm_stroma <- do.call("cbind", d10x.data)


nibone <- read.table("/Data/GSE108891_subset_Iannis/GSE108891_control/niche-col23/GSM2915575_niche-col23.txt.gz")
nipv<- read.table("/Data/GSE108891_subset_Iannis/GSE108891_control/niche-lepr/GSM2915576_niche-lepr.txt.gz")
nivas<- read.table("/Data/GSE108891_subset_Iannis/GSE108891_control/niche-vecad/GSM2915577_niche-vecad.txt.gz")


m_NICHE <-  Read10X(data.dir = "/Data/BMNiche_mouse_samples/NICHE/counts")
mouse_NICHE <- CreateSeuratObject(counts = m_NICHE, project = "mouse sample NICHE", min.cells = 3, min.features = 200)

nCount0.1_0.9 <- quantile(mouse_NICHE$nCount_RNA, c(0.1, 0.9))
nCount0.1_0.9
nFeature0.1_0.9 <- quantile(mouse_NICHE$nFeature_RNA, c(0.1, 0.9))
nFeature0.1_0.9
mouse_NICHE[["percent.mito"]] <-PercentageFeatureSet(mouse_NICHE, pattern = "^mt-")
mouse_NICHEs<- subset(mouse_NICHE, subset = nFeature_RNA > nFeature0.1_0.9[1]  & nFeature_RNA < nFeature0.1_0.9[2] &nCount_RNA > nCount0.1_0.9[1] & nCount_RNA< nCount0.1_0.9[2]& percent.mito < 10)
mouse_NICHEs



bm_str <- CreateSeuratObject(bm_stroma, min.cells = 3, min.features = 200,project = "niche and stroma sample intergration")


BM_list <- list("bm_stroma", "nibone", "nipv", "nivas")
for (i in 1:length(BM_list)) {
  BM_list[[i]] <-eval(substitute(CreateSeuratObject(meta,min.cells = 3, min.features = 200,project = "niche and stroma sample intergration"), list(meta=as.name(BM_list[[i]]))))
}
for (i in 1:length(BM_list)) {
  nCount0.1_0.9 <- quantile(BM_list[[i]]$nCount_RNA, c(0.1, 0.9))
  nCount0.1_0.9
  nFeature0.1_0.9 <- quantile(BM_list[[i]]$nFeature_RNA, c(0.1, 0.9))
  nFeature0.1_0.9
  BM_list[[i]][["percent.mito"]] <-PercentageFeatureSet(BM_list[[i]], pattern = "^mt-")
  cat("dim of the ", i, "-th file is ", dim(BM_list[[i]]) , "\n")
  BM_list[[i]]<- subset(BM_list[[i]], subset = nFeature_RNA > nFeature0.1_0.9[1]  & nFeature_RNA < nFeature0.1_0.9[2] &nCount_RNA > nCount0.1_0.9[1] & nCount_RNA< nCount0.1_0.9[2]& percent.mito < 5)
  cat("dim of the ", i, "-th subset file is ",dim(BM_list[[i]]), "\n")
}

for (i in 1:length(BM_list)) {
  BM_list[[i]] <- SCTransform(BM_list[[i]], verbose = F)
}


Iannis_merged_sub <- merge(nibone, y=c(nipv, nivas))

#### pairwise integratoin #####
### using Iannis dataset as reference: Scadden+Iannis; In-house+Iannis

bm.merge.features <- SelectIntegrationFeatures(object.list = list(bm_str2,Iannis_merged_sub), nfeatures = 3000)
bm.2Da.merged <- PrepSCTIntegration(object.list = list(bm_str2,Iannis_merged_sub), anchor.features =  bm.merge.features) 
bm.2Da.anchors <- FindIntegrationAnchors(object.list = bm.2Da.merged, anchor.features =  bm.merge.features, normalization.method = "SCT")

bm.2Da.merged <- IntegrateData(anchorset = bm.2Da.anchors,normalization.method = "SCT")
bm.2Da.merged <- RunPCA(bm.2Da.merged) %>% RunUMAP(dims =1:50, verbose = F)
bm.2Da.merged <- FindNeighbors(bm.2Da.merged, dims = 1:50) %>% FindClusters(resolution = seq(0.2,1.0,0.2))



Ian_pamp.features <- SelectIntegrationFeatures(object.list = list(mouse_NICHEs,Iannis_merged_sub), nfeatures = 3000)
Ian_pamp.merged <- PrepSCTIntegration(object.list = list(mouse_NICHEs,Iannis_merged_sub), anchor.features =  Ian_pamp.features) 
Ian_pamp.anchors <- FindIntegrationAnchors(object.list = Ian_pamp.merged, anchor.features =  Ian_pamp.features, normalization.method = "SCT")
Ian_pamp.merged <- IntegrateData(anchorset = Ian_pamp.anchors,normalization.method = "SCT")
Ian_pamp.merged <- RunPCA(Ian_pamp.merged) %>% RunUMAP(dims =1:50, verbose = F)
Ian_pamp.merged <- FindNeighbors(Ian_pamp.merged,dims = 1:50) %>% FindClusters(resolution = seq(0.2,1.0,0.2))

#### clusters were annotated manually and EC, MSC clusters were selected ####

### get the ID of EC cells and re-normalize and integrate EC cells only####

Sca_Ian_EC_ID <- names(bm.2Da.merged$orig.ident)[bm.2Da.merged$integrated_snn_res.0.4%in%c(0,4)]
length(Sca_Ian_EC_ID)
Ian_pamp_EC_ID <- names(Ian_pamp.merged$orig.ident)[Ian_pamp.merged$integrated_snn_res.0.4%in%c(0,1,3)]
length(Ian_pamp_EC_ID)

total_ID <- unique(c(Sca_Ian_EC_ID, Ian_pamp_EC_ID))
length(total_ID)

dataset_3_EC_ID <- names(sample_NICHE_merged$orig.ident)[sample_NICHE_merged$integrated_snn_res.0.4%in%c(0,3,10)]
length(intersect(dataset_3_EC_ID, total_ID))
common <- intersect(dataset_3_EC_ID, total_ID)


length(total_ID)
SCadden.EC.ID <- names(sample_NICHE_merged@active.ident)[sample_NICHE_merged$orig.ident=="bm_stroma Scadden"&names(sample_NICHE_merged@active.ident)%in%total_ID]
length(SCadden.EC.ID)
Iannis.EC.ID <- names(sample_NICHE_merged@active.ident)[sample_NICHE_merged$orig.ident%in%c("Col2.3","LEPR","VE-Cad")&names(sample_NICHE_merged@active.ident)%in%total_ID]
cat("Iannis.EC.No", "\n",length(Iannis.EC.ID), "\n")
Pamp.EC.ID<- names(sample_NICHE_merged@active.ident)[sample_NICHE_merged$orig.ident=="mouse sample NICHE"&names(sample_NICHE_merged@active.ident)%in%total_ID]
cat("mouse sample pamplona EC No.","\n",length(Pamp.EC.ID), "\n")
sum(length(Pamp.EC.ID)+length(Iannis.EC.ID)+length(SCadden.EC.ID)-length(total_ID))


###  normalization and integration of EC cells

Scadden.EC <- subset(bm_str, cells = Scadden.EC.ID)
Pamp.EC <- subset(mouse_NICHE, cells = Pamp.EC.ID)
Iannis.EC <- subset(Iannis_merged, cells = Iannis.EC.ID)

EC.integration <- list(Scadden.EC,Iannis.EC,Pamp.EC)
EC.features <- SelectIntegrationFeatures(object.list = EC.integration, nfeatures = 3000)
EC.integration <- PrepSCTIntegration(object.list = EC.integration, anchor.features = EC.features, verbose = F)
EC.integration <- PrepSCTIntegration(object.list = EC.integration, anchor.features = EC.features, verbose = F)
EC.anchors <- FindIntegrationAnchors(object.list = EC.integration, normalization.method = "SCT", anchor.features = EC.features, verbose = F)
EC.integration <- IntegrateData(anchorset = EC.anchors,normalization.method = "SCT", verbose = F)
EC.integration <- RunPCA(EC.integration,verbose = F, npcs = 100) 
EC.integration <- JackStraw(EC.integration,num.replicate = 100, dims = 50)
EC.integration <- ScoreJackStraw(EC.integration, dims = 1:50)


#### level of clustering
EC.integration <- IKAP(EC.integration, out.dir = "../../Results/IKAP_0127_EC.itegration", scale.data = F,find.var.features = F,  r.kmax.est =1)
## best result have 2 cluster, subclsutering of clusterA

EC.clusters <- SplitObject(EC.integration, split.by = "PCK102")

ClusterA <- IKAP(EC.clusters[[1]], out.dir = "../../Results/IKAP_0127_EC.itegration_clusterA2_rerunPCAfor_this_cluster", scale.data = F,find.var.features = F,  r.kmax.est =1)

### repeat subclustering for subclsuters of cluster A until the cluster reach louvain high resolution; see methods in https://www.biorxiv.org/content/10.1101/2021.07.17.452614v1 
### repeat for MSCs from mouse
