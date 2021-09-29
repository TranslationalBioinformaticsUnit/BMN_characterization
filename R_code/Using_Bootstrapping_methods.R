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

### Function 2 rf_preplot: to get metrix iike recall of every cell for plotting #####
##  rf_step1_result: result from Pair select function 
##  obj: same Seurat object for Pair_selection
mEC_rf_preplot <- rf_preplot(mEC_rf_step1, EC.integration)

#add recall to the seurat object
EC.integration$recall_per_cell <- mEC_rf_report$pred_real_sum$recall 

### function 3 Dominant_cluster: determine the dominant cluster for each run #####

mEC_Dominant <- Dominant_cluster(mEC_rf_step1)

###function 5: how many times a clusters remains dominant cluster for cells inside in all runs #####
EC.integration<- stable_times(obj = EC.integration,obj_dominant =mEC_Dominant) ## adding no. of stable times to Seurat object


# PLOTING

##plot UMAP
DimPlot(EC.integration,group.by = "bootstrapped", label = T )+RotatedAxis()+ labs(title = paste0("subclustering of ", prefix))
#plot recall per cell in UMAP
FeaturePlot(EC.integration, features = "recall_per_cell", min.cutoff = 0, max.cutoff = 1)+scale_colour_viridis_c(begin = 0.2)+ RotatedAxis()
#plot the stable times
EC.integration@meta.data %>% dplyr::select(stable_times_of_13_comparisons,bootstrapped) %>% table() %>% 
  data.frame() %>% mutate(stable_times_of_13_comparisons= factor(stable_times_of_13_comparisons)) %>% 
  ggplot()+geom_bar(aes(fill=stable_times_of_13_comparisons, y = Freq, x=bootstrapped),position="fill", stat = "identity")+
  scale_fill_viridis_d(begin = 0.2)+theme_bw()+RotatedAxis()

### dominant plot:
Plot_dominant_umap(pdf_prefix = "mEC", EC.integration, idents = EC.integration$bootstrapped, obj_dominant  = mEC_Dominant, output.dir = "~/Downloads")

####  Assign unstable cluster 
####  using cluster A2 as an example: subcluster A2.2 and cluster A2.3 is not stable, will be assigned to the rest of the clusters 

names <- levels(Cluster_A2@active.ident)
#names[1] is A2.3, names[5] is A2.2
obj_sub_train <- subset(Cluster_A2, idents=names[c(2,3,4,6,7)])
obj_sub_test <- subset(Cluster_A2, idents=names[1])
n.bin=5
n.gene = 50 ### use 50 genes to perform random forest


### function 6: assign unstable clusters ######
##  obj_sub_train: clusters use as a referencing (training dataset)
##  obj_sub_test: unstable cluster to be assigned to clusters in obj_sub_train
##  n1: cluster names in obj_sub_train
##  n2: cluster names in obj_sub_test

### Assign subcluster A2.3 to the rest 
A2_3_list <- pair_rf_assign(obj_sub_test = obj_sub_test,obj_sub_train=obj_sub_train, n1=paste0(names[c(2,3,4,6,7)],collapse=""), n2=names[1])
A2_3_list$pred.table$most_classified  <- apply(A2_3_list$pred.table, 1, function(x)names(which.max(table(t(x))))) ### chose the new cluster names that has been assigned most to each cell 

### Assign subcluster A2.3 to the rest 
obj_sub2 <- subset(Cluster_A2, idents=names[5])
A2_2_list <- pair_rf_assign(obj_sub2,obj_sub_train, n1=paste0(names[c(2,3,4,6,7)],collapse=""), n2=names[5])
A2_2_list$pred.table$most_classified  <- apply(A2_2_list$pred.table, 1, function(x)names(which.max(table(t(x)))))

#### adding metadata to cluster A2 
PC65 <- c(A2_3_list$pred.table$most_classified, A2_2_list$pred.table$most_classified)
names(PC65) <- c(rownames(A2_3_list$pred.table), rownames(A2_2_list$pred.table)) ### new identity for cluster A2.2 and clusterA2.3
PC67 <- Cluster_A2$PC67[-which(colnames(Cluster_A2)%in%names(PC65))] ### remove old identity for cells from cluster A2.2 and cluster A2.3
PC65 <- c(PC65, PC67) ### adding new identity of cells from cluster A2.2 and cluster A2.3
Cluster_A2 <- AddMetaData(Cluster_A2, PC65, col.name = "PC65") ### adding new identity to the seurat object

