library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(randomForest)
require(caTools)
library(zoo)
library(scales)

set.seed(666)

### function 1 Pair selection, select the cluster pairs #########
##obj: a seurat object
##idents: Idententies of each cell from seurat object for testing
Pair_select <- function(obj, idents, fix.training.size = FALSE,n.bin =5,n.gene =20,n.runs =10) {
  Idents(obj) <- idents
  
  names <- levels(obj@active.ident)
  prediction.list <- list()
  cm.all <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(cm.all) <- c("real", "pred", "Freq", "run", "pair") 
  gene.all <- data.frame(matrix(ncol = 0, nrow =n.gene*n.bin))
  if (fix.training.size) {
    if (min(table(obj@active.ident)) < 100) {
      training_size <- floor(min(table(obj@active.ident))*0.8*1.5)
    } else {training_size = min(table(obj@active.ident))*0.8}
  } #### if use fixed training size, determine the size 
  for (i in 1:(length(names)-1)){
    for (j in (i+1):length(names)){
      cat(i,j, "\n")
     cat(names[i],"\t", names[j],"\n")
      obj_sub <- subset(obj, idents  =c(names[i], names[j]))   #### subset the seurat object into two pairs 
      
      if (fix.training.size) {
        all.list <- pair_rf_fix_no(obj_sub, n1=names[i], n2=names[j], training_size) ##### if fixed trainning size has been chosen, choose pair_rf_fix_no
      } else { 
        all.list <- pair_rf(obj_sub, n1=names[i], n2=names[j])#### otherwise, choose pair_rf 
      }
      
      
      all.list$pred.table$real<-obj_sub@active.ident
      predictions <- split(all.list$pred.table, f=all.list$pred.table$real)
      prediction.list <- c(prediction.list,predictions)
      cm.all <- rbind(cm.all, all.list$cm.list)
      #gene.all <- cbind(gene.all, all.list$gene.list)
      cat("@@@@@@@@@@@@@@@@@ finished ",names[i]," ", names[j],"@@@@@@@@@@@@@@\n")
    }
  }
  
  prediction.list.clean <- list()
  
  for (n in names ) {
    print(n)
    index <- names(prediction.list)==n
    one_cluster <- do.call(cbind.data.frame, prediction.list[index])
    colnames(one_cluster) <- gsub(paste0(n,"\\."),"",colnames(one_cluster))
    no_col_rm <- sum(index)-1
    # cols_rm <- n.runs*seq(1,no_col_rm)+1
    if(no_col_rm != 0){
      cols_rm <- (n.runs+1)*seq(from=1,to=no_col_rm)
      one_cluster <- one_cluster[,-cols_rm]
    }
    prediction.list.clean<- c(prediction.list.clean, list(one_cluster))
  }
  names(prediction.list.clean) <- names
  
  pred.final <- list(prediction.list.clean,cm.all,gene.all)
  names(pred.final) <- c("prediction.list.clean", "cm.all","gene.all")
  return(pred.final)
}


pair_rf <- function(obj_sub, n1, n2) { #### build rm model on one pair of clusters
  pred.table <- data.frame(matrix(ncol = n.runs, nrow =dim(obj_sub)[2])) ### predition table, cell as rows and the prediction as columns 
  rownames(pred.table) <- colnames(obj_sub)
  cm.list <- data.frame(matrix(ncol = 5, nrow = 0))    #### store the 5 rm parameters: "real" identity, "pred" identity, "Freq", which "run", "pair" of which two clusters
  colnames(cm.list) <- c("real", "pred", "Freq", "run", "pair") 
  
  for (j in 1:n.runs) {  
   # cat("********* run", j, "************\n")
    bins  <- sample(1:n.bin, dim(obj_sub)[2],replace=T) #### sample 5 bins for testing
    pred <- data.frame(pred = character())
    cm <- data.frame(real=character(),
                     pred=character(), 
                     Freq=numeric()) 
    gene_table <- character()
    for (i in 1:n.bin) {
      #cat("---------------bin", i, "-------------\n")
      obj_sub_train <- subset(obj_sub,cells = colnames(obj_sub)[bins!=i])  #### get the features from training set (80%)
      markers <- FindAllMarkers(object= obj_sub_train, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) #### find all markers from the 
      genes <- markers %>% group_by(cluster) %>% top_n(n=(n.gene/2), wt = avg_log2FC)  ####using 20 genes for clustering, 10 from each cluster
      cat("dim of training dataset:",dim(obj_sub_train), "\n")
      cat("dim of training markers:",dim(markers), "\n")
      rm(obj_sub_train)
      gc()
      index <- which(rownames(obj_sub)%in%genes$gene) ### row numbers of the 20 genes
      sample <- GetAssayData(obj_sub, slot = "scale.data", assay = "integrated")[index,]  ##### get scale dataset from integrated object (only 20 genes) -------can be change later--------------------
      sample <- t(sample)
      sample <- as.data.frame(sample)
      if (ncol(sample)<20){
        cat("Less than 20 marker genes find, use",ncol(sample)," genes instead\n")
        #n.gene=ncol(sample)
      }
      sample <- cbind(sample,Idents(obj_sub))   ###########append cluster identity to the training matrix
    
      dimnames(sample)[[2]][ncol(sample)] <- "Cluster"    #### Cluster is the column for rm to classify
      cat("dim of the sample:", dim(sample),"\n")
      cat("dim of the genes:", dim(genes),"\n")
      colnames(sample)[1:ncol(sample)-1]<- paste0("gene_",colnames(sample))[1:ncol(sample)-1]  ##### incase genes started with numbers 
      if (length(grep("-", colnames(sample)))!=0) {   ###### change - to "_"
        index.gene <- grep("-", colnames(sample))
        cat("         !!!gene names has '-' :", colnames(sample)[index.gene], "!!!         \n")
        colnames(sample)[index.gene] <- gsub("-", "_",colnames(sample)[index.gene])
      }
      test <- sample[bins==i,]
      train <-sample[bins!=i,]
      cat("dim of testing dataset:",dim(test), "\n")
      
      
      rf <- randomForest(Cluster ~.,ntrees = 1000,data=train)  #### use 1000 trees 
      pred.temp<- predict(rf, newdata = test[-ncol(test)]) 
      cm <- rbind(cm,as.data.frame(table(test[order(match(test, names(pred.temp)))][,(ncol(test))], pred.temp)))
      gene_table <- c(gene_table, genes$gene)
      pred <- rbind(pred,as.data.frame(pred.temp))
      #cat("---------------finished bin", i, "-------------\n")
    }
    colnames(pred) <- "pred"
    pred$cells <- rownames(pred)
    pred <- pred[order(match(pred$cells, rownames(pred.table))),]
    pred.table[,j]<- pred$pred
    colnames(pred.table)[j] <- paste0(n1, n2, "_run", j)
    
    colnames(cm) <- c("real", "pred", "Freq") 
    cm <- aggregate(Freq ~ real + pred, data = cm, sum)
    cm$run <- j
    cm$pair <- paste0(n1, n2)
    
    cm.list <- rbind(cm.list, cm)
    
   
    cat("\n\n********* finished run", j, "************\n", "\n")
  }
  
  pred.all <- list(pred.table,cm.list)
  
  names(pred.all) <- c("pred.table", "cm.list")
  return(pred.all)
}

pair_rf_fix_no <- function(obj_sub, n1, n2, traning_size) {
  pred.table <- data.frame(matrix(ncol = n.runs, nrow =dim(obj_sub)[2]))
  rownames(pred.table) <- colnames(obj_sub)
  cm.list <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(cm.list) <- c("real", "pred", "Freq", "run", "pair") 
  
  
  for (j in 1:n.runs) {  
    cat("********* run", j, "************\n")
    bins  <- sample(1:n.bin, dim(obj_sub)[2],replace=T)
    pred <- data.frame(pred = character())
    cm <- data.frame(real=character(),
                     pred=character(), 
                     Freq=numeric()) 
    gene_table <- character()
    for (i in 1:n.bin) {
      cat("---------------bin", i, "-------------\n")
      obj_sub_train <- subset(obj_sub,cells = colnames(obj_sub)[bins!=i])
      obj_sub_train_list <- SplitObject(obj_sub_train, split.by = "ident") 
      train_list <- lapply(obj_sub_train_list, function(x){sample(1:dim(x)[2], size=training_size, replace=T)})
      
      count1 <- obj_sub_train_list[[1]]@assays$integrated@scale.data[, train_list[[1]]]
      colnames(count1) <- paste0(colnames(count1),"_", seq(1:dim(count1)[2]))
      meta1 <- as.data.frame(obj_sub_train_list[[1]]@active.ident[train_list[[1]]])
      colnames(meta1) <- "ident"
      rownames(meta1) <- colnames(count1)
      
      count2 <- obj_sub_train_list[[2]]@assays$integrated@scale.data[, train_list[[2]]]
      colnames(count2) <- paste0(colnames(count2),"_", seq(1:dim(count2)[2]))
      meta2 <- as.data.frame(obj_sub_train_list[[2]]@active.ident[train_list[[2]]])
      colnames(meta2) <- "ident"
      rownames(meta2) <- colnames(count2)
      
      obj_sub_train_2 <- CreateSeuratObject(counts = cbind(count1, count2), meta.data = rbind(meta1, meta2))
      Idents(obj_sub_train_2) <- obj_sub_train_2$ident
      
      markers <- FindAllMarkers(object= obj_sub_train_2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, assay = "RNA", slot = "data")
      genes <- markers %>% group_by(cluster) %>% top_n(n=(n.gene/2), wt = avg_log2FC)
      cat("dim of training dataset:",dim(obj_sub_train_2), "\n")
      
      index <- which(rownames(obj_sub)%in%genes$gene)
      sample <- GetAssayData(obj_sub, slot = "scale.data", assay = "integrated")[index,]
      sample <- t(sample)
      sample <- as.data.frame(sample)
      sample <- cbind(sample,Idents(obj_sub))
      dimnames(sample)[[2]][21] <- "Cluster"
      
      colnames(sample)[1:n.gene]<- paste0("gene_",colnames(sample))[1:n.gene]
      if (length(grep("-", colnames(sample)))!=0) {
        index.gene <- grep("-", colnames(sample))
        cat("         !!!gene names has '-' :", colnames(sample)[index.gene], "!!!         \n")
        colnames(sample)[index.gene] <- gsub("-", "_",colnames(sample)[index.gene])
      }
      
      test <- sample[bins==i,]
    
      
      train <- GetAssayData(obj_sub_train_2, slot = "data", assay = "RNA")[index,]
      train <- t(train)
      train <- as.data.frame(train)
      train <- cbind(train,Idents(obj_sub_train_2))
      colnames(train) <- colnames(test)
      
      rm(obj_sub_train, count1, count2,obj_sub_train_2, obj_sub_train_list)
      gc()
      
      cat("dim of testing dataset:",dim(test), "\n")
      
      
      rf <- randomForest(Cluster ~.,ntrees = 1000,data=train)
      pred.temp<- predict(rf, newdata = test[-(n.gene+1)])
      cm <- rbind(cm,as.data.frame(table(test[order(match(test, names(pred.temp)))][,(n.gene+1)], pred.temp)))
      gene_table <- c(gene_table, genes$gene)
      pred <- rbind(pred,as.data.frame(pred.temp))
      cat("---------------finished bin", i, "-------------\n")
    }
    colnames(pred) <- "pred"
    pred$cells <- rownames(pred)
    pred <- pred[order(match(pred$cells, rownames(pred.table))),]
    pred.table[,j]<- pred$pred
    colnames(pred.table)[j] <- paste0(n1, n2, "_run", j)
    
    colnames(cm) <- c("real", "pred", "Freq") 
    cm <- aggregate(Freq ~ real + pred, data = cm, sum)
    cm$run <- j
    cm$pair <- paste0(n1, n2)
    
    cm.list <- rbind(cm.list, cm)
    
    
    
    cat("\n\n********* finished run", j, "************\n", "\n")
  }
  
  pred.all <- list(pred.table,cm.list)
  names(pred.all) <- c("pred.table", "cm.list")
  
  return(pred.all)
}

### Function 2 rf_preplot: to get metrix lile reall of every cell for plotting ########
### rf_step1_result: result from Pair select function 
### obj: same Seurat object for Pair_selection
add_recall_per_cell <- function(rf_step1_result,obj,seurat_colname) {
  names <- unlist(sapply(rf_step1_result$prediction.list.clean, rownames, simplify = T))
  pred_real_sum <- data.frame(matrix(nrow = dim(obj)[2], ncol = 3))
  rownames(pred_real_sum) <- names
  pred_real_sum[, 1] <- names
  recall_all_avg_pcell <- numeric()
  recall_avg_ppari_pcell_all <- list()
  sd_all <- numeric()
  for (i in rf_step1_result$prediction.list.clean) {
    counts <- apply(i[,1:(dim(i)[2]-1)],1,function(x)length(grep(unique(i$real),x))) # how many times a cell has been correctly assigned
    recall_avg_pcell <- counts/length(colnames(i)) # how many times a cell has been correctly classified/total runs
    recall_all_avg_pcell <- c(recall_all_avg_pcell, recall_avg_pcell)
    cat(length(recall_avg_pcell), "\n")
    
    if(dim(i)[2] <20){
      counts_per_pair <- data.frame(apply(i[,1:(dim(i)[2]-1)],1,function(y){rollapply(y, width=10, FUN=function(x)length(grep(unique(i$real),x)), by=10)}))
    } else {
      counts_per_pair <- data.frame(t(apply(i[,1:(dim(i)[2]-1)],1,function(y){rollapply(y, width=10, FUN=function(x)length(grep(unique(i$real),x)), by=10)})))
    }
    colnames(counts_per_pair) <- unique(sapply(strsplit(colnames(i[,1:(dim(i)[2]-1)]),"_" ), function(x)x[[1]], simplify = T))
    recall_avg_ppari_pcell <- counts_per_pair/10
    sd_per_cell_all_run <- apply(recall_avg_ppari_pcell,1, sd)
    sd_all <- c(sd_all, sd_per_cell_all_run)
    
    recall_avg_ppari_pcell$sd <- sd_per_cell_all_run
    recall_avg_ppari_pcell_all <- c(recall_avg_ppari_pcell_all, recall_avg_ppari_pcell)
    
  }
  pred_real_sum[,2] <- recall_all_avg_pcell
  pred_real_sum[,3]  <- sd_all
  colnames(pred_real_sum) <- c("cell", "recall", "sd")
  pred_real_sum <- pred_real_sum[order(match(pred_real_sum$cell, colnames(obj))),]
  
  stat_rf <- list(pred_real_sum, recall_avg_ppari_pcell_all)
  names(stat_rf) <- c("pred_real_sum", "recall_avg_ppari_pcell_all")
  obj[[seurat_colname]] <-  stat_rf$pred_real_sum$recall
  return(obj)
}

### function 3 Dominant_cluster
######### determine the dominat cluster for each run  ### need no.runs for width
Dominant_cluster <- function(rf_step1_result, n.runs=10) {
  dominant_cluster_name_all <- list()
  dominant_cluster_count_all <- list()
  for (i in rf_step1_result$prediction.list.clean) { 
    if (length(rf_step1_result$prediction.list.clean) <3){
      dominant_cluster_name <- data.frame(apply(i[,1:(dim(i)[2]-1)],1,function(y){rollapply(y, width=n.runs, FUN=function(x)names(which.max(table(t(x)))), by=n.runs)}))
      colnames(dominant_cluster_name) <- unique(sapply(strsplit(colnames(i[,1:(dim(i)[2]-1)]),"_run" ), function(x)x[[1]], simplify = T))
      dominant_cluster_count <- data.frame(apply(i[,1:(dim(i)[2]-1)],1,function(y){rollapply(y, width=n.runs, FUN=function(x)max(table(t(x))), by=n.runs)}))
      colnames(dominant_cluster_count) <- unique(sapply(strsplit(colnames(i[,1:(dim(i)[2]-1)]),"_run" ), function(x)x[[1]], simplify = T))
    } else {
      dominant_cluster_name <- data.frame(t(apply(i[,1:(dim(i)[2]-1)],1,function(y){rollapply(y, width=n.runs, FUN=function(x)names(which.max(table(t(x)))), by=n.runs)})))
      colnames(dominant_cluster_name) <- unique(sapply(strsplit(colnames(i[,1:(dim(i)[2]-1)]),"_run" ), function(x)x[[1]], simplify = T))
      dominant_cluster_count <- data.frame(t(apply(i[,1:(dim(i)[2]-1)],1,function(y){rollapply(y, width=n.runs, FUN=function(x)max(table(t(x))), by=n.runs)})))
      colnames(dominant_cluster_count) <- unique(sapply(strsplit(colnames(i[,1:(dim(i)[2]-1)]),"_run" ), function(x)x[[1]], simplify = T))
    }
    dominant_cluster_name_all <- c(dominant_cluster_name_all, list(dominant_cluster_name))
    dominant_cluster_count_all <- c(dominant_cluster_count_all, list(dominant_cluster_count_all))
  }
  dominant_cluster_all <- list(dominant_cluster_name_all,dominant_cluster_count_all)
  names(dominant_cluster_all) <- c("dominant_cluster_name_all","dominant_cluster_count_all")
  names(dominant_cluster_all$dominant_cluster_name_all) <- names(rf_step1_result$prediction.list.clean)
  names(dominant_cluster_all$dominant_cluster_count_all) <- names(rf_step1_result$prediction.list.clean)
  return(dominant_cluster_all)
}


###  Function 4 plot the dominant clusters in different UMAP plots ##### 
##   pdf_prefix: name of the seurat obj
##   output.dir : output directory
##   annotation bar and text position can be adjusted manually if not showing properly
Plot_dominant_umap <- function(pdf_prefix, obj, idents, obj_dominant,output.dir) {
  out.dir <- paste0(output.dir,"/dominant_cluster_plots_", pdf_prefix)
  dir.create(out.dir)
  
  Idents(obj) <- idents
  pal <- as.character(hue_pal()(length(levels(Idents(obj)))))
  colors <- data.frame(cluster = sort(levels(Idents(obj))),cols= pal)
  for (j in 1:length(obj_dominant$dominant_cluster_name_all)) {
    
    plot.list <- list()
    for (c in 1:ncol(obj_dominant$dominant_cluster_name_all[[j]])) {
      cluster <-sort(unique(as.character(obj_dominant$dominant_cluster_name_all[[j]][,c])))
      if(length(cluster)==1){
        cluster1_id <- rownames(obj_dominant$dominant_cluster_name_all[[j]])[obj_dominant$dominant_cluster_name_all[[j]][,c] ==cluster[1]]
        cluster1_col <- as.character(colors$cols[grep(cluster[1], colors$cluster)])
        plot.list[[c]] <-  DimPlot(obj, label = T, label.size = 4, cells.highlight =list(cluster1_id), pt.size = 0.8, repel = T)+
          scale_color_manual(labels =c("others",as.character(cluster[1])), values = c("grey", cluster1_col)) +
          labs(title = colnames(obj_dominant$dominant_cluster_name_all[[j]])[c])

      }
      else{
        cluster1_id <- rownames(obj_dominant$dominant_cluster_name_all[[j]])[obj_dominant$dominant_cluster_name_all[[j]][,c] ==cluster[1]]
        cluster2_id <- rownames(obj_dominant$dominant_cluster_name_all[[j]])[obj_dominant$dominant_cluster_name_all[[j]][,c] ==cluster[2]]
        cluster1_col <- as.character(colors$cols[grep(cluster[1], colors$cluster)])
        cluster2_col <- as.character(colors$cols[grep(cluster[2], colors$cluster)])
        plot.list[[c]] <-  DimPlot(obj, label = T, label.size = 4, cells.highlight =list(cluster1_id,cluster2_id), pt.size = 0.8, repel =T)+
          scale_color_manual(labels =c("others",as.character(cluster[2]), as.character(cluster[1])), values = c("grey",cluster2_col, cluster1_col)) +
          labs(title = colnames(obj_dominant$dominant_cluster_name_all[[j]])[c])+
          annotation_raster(c(cluster1_col, cluster2_col)[c(rep(1, length(cluster1_id)), rep(2,length(cluster2_id)))],xmin = 6.5, xmax =Inf,ymin=3,ymax = 5.8)+ ### annotation bar position: can be adjusted manually 
          annotate("text", label = round(length(cluster1_id)/(length(cluster1_id)+length(cluster2_id)),2), x = 7.14, y = 6.05, size = 3, colour = cluster1_col)+ ### annotation text position: can be adjusted manually 
          annotate("text", label = round(length(cluster2_id)/(length(cluster1_id)+length(cluster2_id)),2), x = 7.14, y = 2.75, size = 3, colour = cluster2_col)### annotation text position: can be adjusted manually 
      }
    }
    
    file_name <- paste0(out.dir,"/",pdf_prefix, "_",names(obj_dominant$dominant_cluster_name_all[j]),".pdf")
    if (length(plot.list) < 2){
      nrows = 1
    }else {
      nrows=2
    }
    pdf(file=file_name, width = 5*length(plot.list)/nrows, height = 4.85*nrows)
    do.call("grid.arrange", c(plot.list, nrow=nrows, top=names(obj_dominant$dominant_cluster_name_all[j])))
    dev.off()
  }
}


### function 5: how many times a clusters remains dominant cluster for cells inside in all runs ####### 
stable_times <- function(obj_dominant, obj,meta.names) {
  obj_stable_times <-numeric() 
  for (i in 1:length(obj_dominant$dominant_cluster_name_all)){
    obj_dominant$dominant_cluster_name_all[[i]]$stable_times <- apply(obj_dominant$dominant_cluster_name_all[[i]], 1, function(x)sum(x==names(obj_dominant$dominant_cluster_name_all)[i]))
    obj_stable_times <- c(obj_stable_times, t(obj_dominant$dominant_cluster_name_all[[i]]$stable_times))
  }
  names_stable_times <- unlist(sapply(obj_dominant$dominant_cluster_name_all, rownames, simplify = T))
  names(obj_stable_times) <- names_stable_times
  #meta.names <- paste0("stable_times_of_",length(obj_dominant$dominant_cluster_name_all)-1,"_comparisons")
  obj <- AddMetaData(obj, obj_stable_times, col.name = meta.names)
  return(obj)
}

### function 6: assign unstable clusters ######
##  obj_sub_train: clusters use as a referencing (training dataset)
##  obj_sub_test: unstable cluster to be assigned to clusters in obj_sub_train
##  n1: cluster names in obj_sub_train
##  n2: cluster names in obj_sub_test
pair_rf_assign <- function(obj_sub_train,obj_sub_test,n1, n2,n.gene=no.genes,n.bin=no.bins) {
  pred.table <- data.frame(matrix(ncol = n.runs, nrow =dim(obj_sub_test)[2]))
  rownames(pred.table) <- colnames(obj_sub_test)
  
  
  #############fixed training set############ 
  markers <- FindAllMarkers(object= obj_sub_train, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  genes <- markers %>% group_by(cluster) %>% top_n(n=n.gene/length(levels(obj_sub_train@active.ident)), wt = avg_log2FC)
  cat("dim of training dataset:",dim(obj_sub_train), "\n")
  index <- which(rownames(obj_sub_train)%in%genes$gene)
  sample <- GetAssayData(obj_sub_train, slot = "scale.data", assay = "integrated")[index,]
  sample <- t(sample)
  sample <- as.data.frame(sample)
  colnames(sample)<- paste0("gene_",colnames(sample))
  Cluster <- Idents(obj_sub_train)
  sample <- cbind(sample,Cluster)
  
  if (n.gene > length(index)) {
    n.gene <- length(index)
    cat("at least on cluster has less than",ceiling(n.gene/length(levels(obj_sub_train@active.ident))), " genes found, use ", n.gene," instead")
    
  }
  
  
  cm.list <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(cm.list) <- c("real", "pred", "Freq", "run", "pair") 
 
  
  
  if (length(grep("-", colnames(sample)))!=0) {
    index.gene <- grep("-", colnames(sample))
    cat("         !!!gene names has '-' :", colnames(sample)[index.gene], "!!!         \n")
    colnames(sample)[index.gene] <- gsub("-", "_",colnames(sample)[index.gene])
  }
  
  train <-sample 
  rf <- randomForest(Cluster ~.,ntrees = 1000,data=train) ##### fixed classifier
  ###########################end#############################################
  for (j in 1:n.runs) {  
    cat("********* run", j, "************\n")
    bins  <- sample(1:n.bin, dim(obj_sub_test)[2],replace=T)
    pred <- data.frame(pred = character())
    cm <- data.frame(real=character(),
                     pred=character(), 
                     Freq=numeric()) 
   
    
    for (i in 1:n.bin) {
      cat("---------------bin", i, "-------------\n")
      
      test_all <- GetAssayData(obj_sub_test, slot = "scale.data", assay = "integrated")[index,] %>% t() %>% data.frame() 
      test_all$Cluster <- Idents(obj_sub_test)
      colnames(test_all) <- colnames(train)
      test <- test_all[bins==i,]
      cat("dim of testing dataset:",dim(test), "\n")
      
      pred.temp<- predict(rf, newdata = test[-(n.gene+1)])
      cm <- rbind(cm,as.data.frame(table(test[order(match(test, names(pred.temp)))][,(n.gene+1)], pred.temp)))
 
      pred <- rbind(pred,as.data.frame(pred.temp))
      cat("---------------finished bin", i, "-------------\n")
    }
    colnames(pred) <- "pred"
    pred$cells <- rownames(pred)
    pred <- pred[order(match(pred$cells, rownames(pred.table))),]
    pred.table[,j]<- pred$pred
    colnames(pred.table)[] <- paste0(n1, n2, "_run", j)
    
    colnames(cm) <- c("real", "pred", "Freq") 
    cm <- aggregate(Freq ~ real + pred, data = cm, sum)
    cm$run <- j
    cm$pair <- paste0(n1, n2)
    
    cm.list <- rbind(cm.list, cm)
    
   
    cat("\n\n********* finished run", j, "************\n", "\n")
  }

  pred.all <- list(pred.table,cm.list)
  
  names(pred.all) <- c("pred.table", "cm.list")
  return(pred.all)
}

assign_unstable_clusters <- function(obj, idents,unstable_clusters, no.bins=5, no.genes =50,new_ident_name="assigned_new_ident") {
  Idents(obj) <- idents
  names <- levels(obj@active.ident)
  i_train <- !levels(obj@active.ident) %in% unstable_clusters
  obj_sub_train <- subset(obj, idents=names[i_train])
  
  test_ident_all <- character(0)
  for (x in unstable_clusters) {
    cat("assigning cluster",x," \n")
    i_test <- levels(obj@active.ident) %in% x
    obj_sub_test <- subset(obj, idents=names[i_test])
    
    ### function 6: assign unstable clusters ######
    ##  obj_sub_train: clusters use as a referencing (training dataset)
    ##  obj_sub_test: unstable cluster to be assigned to clusters in obj_sub_train
    ##  n1: cluster names in obj_sub_train
    ##  n2: cluster names in obj_sub_test
    
    ### Assign subcluster A2.3 to the rest 
    assign_list <- pair_rf_assign(obj_sub_train=obj_sub_train, obj_sub_test = obj_sub_test,n1=paste0(names[i_train],collapse=""), n2=names[i_test],n.gene=no.genes,n.bin=no.bins)
    
    assign_list$pred.table$most_classified  <- apply(assign_list$pred.table, 1, function(x)names(which.max(table(t(x))))) ### chose the new cluster names that has been assigned most to each cell 
    
    test_ident <- assign_list$pred.table$most_classified
    names(test_ident) <- rownames(assign_list$pred.table)
    test_ident_all <- c(test_ident_all,test_ident)
  }
  
  train_ident <- as.character(obj@active.ident)[obj@active.ident%in%names[i_train]]
  names(train_ident) <- names(obj@active.ident[obj@active.ident%in%names[i_train]])
  New_ident <- c(train_ident,test_ident_all) 
  obj <- AddMetaData(obj, New_ident, col.name = new_ident_name)
  return(obj)
}
