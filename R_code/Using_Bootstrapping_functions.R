source("Bootstrapping_functions.R")
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)

EC.integration@active.assay <- "integrated" ## The target seurat object
## removing unwanted assays to save memory
EC.integration@assays$RNA <- NULL
EC.integration@assays$SCT <- NULL

# Step 1: get the rm result from pair wise comparison with the paramters below
n.bin =5 ### using 20% as testing, 80% as trainning
n.gene =20 ### using 20 genes to classify 2 clusters with rf
n.runs =10 ### repeat the process for 10 times

### Function 1 Pair selection, select the cluster pairs ####
##  obj: a seurat object
##  idents: Identities of each cell from seurat object for testing
mEC_rf_step1 <- Pair_select(EC.integration, idents=EC.integration$bootstrapped)  

### Function 2 add_recall_per_cell: add plotting matrix "recall_per_cell" to the suerat object
EC.integration <- add_recall_per_cell(rf_step1_result = mEC_rf_step1,obj = EC.integration,
                                      seurat_colname = "recall_per_cell")

### function 3 Dominant_cluster: determine the dominant cluster for each run #####
mEC_Dominant <- Dominant_cluster(mEC_rf_step1)

###function 4: how many times a clusters remains dominant cluster for cells inside in all runs #####
EC.integration<- stable_times(obj = EC.integration,obj_dominant =mEC_Dominant,meta.names = "stable_times") ## adding no. of stable times (#correct) to Seurat object meta.names: seurat col names to store no. of stable times

# PLOTING

##plot UMAP
DimPlot(EC.integration,group.by = "bootstrapped", label = T )+RotatedAxis()+ labs(title = paste0("subclustering of ", prefix))
#plot recall per cell in UMAP
FeaturePlot(EC.integration, features = "recall_per_cell", min.cutoff = 0, max.cutoff = 1)+scale_colour_viridis_c(begin = 0.2)+ RotatedAxis()
#plot the stable times
EC.integration@meta.data %>% dplyr::select(stable_times,bootstrapped) %>% table() %>% 
  data.frame() %>% mutate(stable_times= factor(sstable_times)) %>% 
  ggplot()+geom_bar(aes(fill=stable_times, y = Freq, x=bootstrapped),position="fill", stat = "identity")+
  scale_fill_viridis_d(begin = 0.2)+theme_bw()+RotatedAxis()

### dominant plot:
Plot_dominant_umap(pdf_prefix = "mEC", EC.integration, idents = EC.integration$bootstrapped, obj_dominant  = mEC_Dominant, output.dir = "~/Downloads")


####  Assign unstable cluster 
####  using cluster A2 as an example: subcluster A2.2 and cluster A2.3 is not stable, will be assigned to the rest of the clusters 
#obj: seurat object 
#unstable_clusters: the clusters that are not stable and you wish to asign to the rest 
# idents: cell annoation of the all clusters including the unsatble ones
# no bins=5: using 20% as testing, 80% as trainning
#no.genes  = 50: using 50 genes to train a random forest classifier to assgin the unstable clusters one by one.
Cluster_A2 <- assign_unstable_clusters(obj = Cluster_A2,unstable_clusters = c("A2.2", "A2.3"),idents = Cluster_A2@active.ident,no.bins = 5,no.genes  = 50)




