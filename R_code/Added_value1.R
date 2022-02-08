library(Seurat)
library(scales)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(openxlsx)


### define markers identified by on dataset
EC.list$`bm_stroma Scadden`@active.assay <- "RNA"
Idents(EC.list$`bm_stroma Scadden`)<-EC.list$`bm_stroma Scadden`$bootstrapped
#EC.list$`bm_stroma Scadden` <- subset(EC.list$`bm_stroma Scadden`, features = )
Scadden.EC.marker <- FindAllMarkers(EC.list$`bm_stroma Scadden`, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,assay = "SCT", slot = "scale.data")

EC.integration@active.assay <- "integrated"
EC.integration@active.ident <- as.factor(EC.integration$bootstrapped)
EC.markers <- FindAllMarkers(EC.integration,only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

EC.integration@active.assay <- "integrated"
EC.integration@active.ident <- as.factor(EC.integration$bootstrapped)
EC.markers <- FindAllMarkers(EC.integration,only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
overlap_matrix_EC_scadden <- matrix(nrow = 14, ncol = 4)
overlap_matrix_EC_scadden <- data.frame(overlap_matrix_EC_scadden)
colnames(overlap_matrix_EC_scadden) <-c("overlapped", "all EC", "single dataset", "cluster") 
overlap_matrix_EC_scadden$cluster<- unique(EC.markers$cluster)



overlap_matrix_EC_scadden <- matrix(nrow = 14, ncol = 4)
overlap_matrix_EC_scadden <- data.frame(overlap_matrix_EC_scadden)
colnames(overlap_matrix_EC_scadden) <-c("overlapped", "all EC", "single dataset", "cluster") 
overlap_matrix_EC_scadden$cluster<- unique(EC.markers$cluster)

for (i in 1:14){
  cl <- overlap_matrix_EC_scadden$cluster[i]
  overlap_matrix_EC_scadden[i, 1] <- length(intersect(Scadden.EC.marker$gene[Scadden.EC.marker$cluster==cl],EC.markers$gene[EC.markers$cluster==cl]))
  overlap_matrix_EC_scadden[i, 2] <- sum(EC.markers$cluster==cl)
  overlap_matrix_EC_scadden[i, 3] <- sum(Scadden.EC.marker$cluster==cl)
}

overlap_matrix_EC_scadden$shared_pt <- overlap_matrix_EC_scadden$overlapped/overlap_matrix_EC_scadden$`all EC`
overlap_matrix_EC_scadden$dataset <- "Scadden"

overlap_matrix_EC <- overlap_matrix_EC_scadden

EC.list$`VE-Cad`@active.assay <- "RNA"
Idents(EC.list$`VE-Cad`)<-EC.list$`VE-Cad`$bootstrapped
Aifantis.EC.marker <- FindAllMarkers(EC.list$`VE-Cad`, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,assay = "SCT", slot = "scale.data")

for (i in 1:14){
  cl <- as.character(overlap_matrix_EC_scadden$cluster[i])
  overlap_matrix_EC_scadden[i, 1] <- length(intersect(Aifantis.EC.marker$gene[Aifantis.EC.marker$cluster==cl],EC.markers$gene[EC.markers$cluster==cl]))
  overlap_matrix_EC_scadden[i, 2] <- sum(EC.markers$cluster==cl)
  overlap_matrix_EC_scadden[i, 3] <- sum(Aifantis.EC.marker$cluster==cl)
}

overlap_matrix_EC_scadden$shared_pt <- overlap_matrix_EC_scadden$overlapped/overlap_matrix_EC_scadden$`all EC`
overlap_matrix_EC_scadden$dataset <- "Aifantis"

overlap_matrix_EC <- rbind(overlap_matrix_EC, overlap_matrix_EC_scadden)

EC.list$`mouse sample NICHE`@active.assay <- "RNA"
Idents(EC.list$`mouse sample NICHE`)<-EC.list$`mouse sample NICHE`$bootstrapped
Pamplona.EC.marker <- FindAllMarkers(EC.list$`mouse sample NICHE`, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,assay = "SCT", slot = "scale.data")

for (i in 1:14){
  cl <- overlap_matrix_EC_scadden$cluster[i]
  overlap_matrix_EC_scadden[i, 1] <- length(intersect(Pamplona.EC.marker$gene[Pamplona.EC.marker$cluster==cl],EC.markers$gene[EC.markers$cluster==cl]))
  overlap_matrix_EC_scadden[i, 2] <- sum(EC.markers$cluster==cl)
  overlap_matrix_EC_scadden[i, 3] <- sum(Pamplona.EC.marker$cluster==cl)
}

overlap_matrix_EC_scadden$shared_pt <- overlap_matrix_EC_scadden$overlapped/overlap_matrix_EC_scadden$`all EC`
overlap_matrix_EC_scadden$dataset <- "Pamplona"

overlap_matrix_EC <- rbind(overlap_matrix_EC, overlap_matrix_EC_scadden)


### Added value 1 for MSC

MSC.list <- SplitObject(MSC.integration_sub,split.by = "orig.ident")

MSC.markers <- FindAllMarkers(MSC.integration_sub,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

MSC.list$`bm_stroma Scadden`@active.assay <- "SCT"
Idents(MSC.list$`bm_stroma Scadden`)<-MSC.list$`bm_stroma Scadden`$bootstrapped
Scadden.MSC.marker <- FindAllMarkers(MSC.list$`bm_stroma Scadden`, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, verbose = F,assay = "SCT", slot = "scale.data")
MSC.list$LEPR@active.assay <- "SCT"
Idents(MSC.list$LEPR)<-MSC.list$LEPR$bootstrapped
Aifantis.MSC.marker <- FindAllMarkers(MSC.list$LEPR, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, verbose = F,assay = "SCT", slot = "scale.data")
MSC.list$`mouse sample NICHE`@active.assay <- "SCT"
Idents(MSC.list$`mouse sample NICHE`)<-MSC.list$`mouse sample NICHE`$bootstrapped
Pamplona.MSC.marker <- FindAllMarkers(MSC.list$`mouse sample NICHE`, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, verbose = F,assay = "SCT", slot = "scale.data")

overlap_matrix_EC_scadden <- matrix(nrow = 11, ncol = 4)
overlap_matrix_EC_scadden <- data.frame(overlap_matrix_EC_scadden)
colnames(overlap_matrix_EC_scadden) <-c("overlapped", "all EC", "single dataset", "cluster") 
overlap_matrix_EC_scadden$cluster<- unique(MSC.markers.sub$cluster)
for (i in 1:11){
  cl <- as.character(overlap_matrix_EC_scadden$cluster[i])
  overlap_matrix_EC_scadden[i, 1] <- length(intersect(Scadden.MSC.marker$gene[Scadden.MSC.marker$cluster==cl],MSC.markers$gene[MSC.markers$cluster==cl]))
  overlap_matrix_EC_scadden[i, 2] <- sum(MSC.markers$cluster==cl)
  overlap_matrix_EC_scadden[i, 3] <- sum(Scadden.MSC.marker$cluster==cl)
}
overlap_matrix_EC_scadden$shared_pt <- overlap_matrix_EC_scadden$overlapped/overlap_matrix_EC_scadden$`all EC`
overlap_matrix_EC_scadden$dataset <- "Scadden"
overlap_matrix_MSC <- overlap_matrix_EC_scadden

for (i in 1:11){
  cl <- as.character(overlap_matrix_EC_scadden$cluster[i])
  overlap_matrix_EC_scadden[i, 1] <- length(intersect(Aifantis.MSC.marker$gene[Aifantis.MSC.marker$cluster==cl],MSC.markers$gene[MSC.markers$cluster==cl]))
  overlap_matrix_EC_scadden[i, 2] <- sum(MSC.markers$cluster==cl)
  overlap_matrix_EC_scadden[i, 3] <- sum(Aifantis.MSC.marker$cluster==cl)
}
overlap_matrix_EC_scadden$shared_pt <- overlap_matrix_EC_scadden$overlapped/overlap_matrix_EC_scadden$`all EC`
overlap_matrix_EC_scadden$dataset <- "Aifantis"
overlap_matrix_MSC <- rbind(overlap_matrix_MSC, overlap_matrix_EC_scadden)

for (i in 1:11){
  cl <- as.character(overlap_matrix_EC_scadden$cluster[i])
  overlap_matrix_EC_scadden[i, 1] <- length(intersect(Pamplona.MSC.marker$gene[Pamplona.MSC.marker$cluster==cl],MSC.markers$gene[MSC.markers$cluster==cl]))
  overlap_matrix_EC_scadden[i, 2] <- sum(MSC.markers$cluster==cl)
  overlap_matrix_EC_scadden[i, 3] <- sum(Pamplona.MSC.marker$cluster==cl)
}
overlap_matrix_EC_scadden$shared_pt <- overlap_matrix_EC_scadden$overlapped/overlap_matrix_EC_scadden$`all EC`
overlap_matrix_EC_scadden$dataset <- "Pamplona"
overlap_matrix_MSC <- rbind(overlap_matrix_MSC, overlap_matrix_EC_scadden)



####Venn Diagram#####

library(VennDiagram)
out.dir <- paste0("/added_value/point1/Venn_diagram_v2/EC/")
EC.markers.in <- EC.markers
for (i in unique(EC.markers.in$cluster)){
  cat (i)
  venn.diagram(
    x = list(EC.markers.in$gene[EC.markers.in$cluster==i],Pamplona.EC.marker$gene[Pamplona.EC.marker$cluster==i],Scadden.EC.marker$gene[Scadden.EC.marker$cluster==i],Aifantis.EC.marker$gene[Aifantis.EC.marker$cluster==i]),
    category.names = c("integrated" , "In-house", "Baryawana", "Tikhonova"),
    filename = paste0(out.dir,"/",prefix, "cluster_",i, ".png"),
    output=TRUE,
    imagetype="png" ,
    main = paste0("added_value1: Cluster ", i),
    main.cex =1,
    main.fontfamily = "Arial",
    height = 2000, 
    width = 3000, resolution = 500,compression = "lzw",
    lwd = 1,
    col=c("wheat", 'light blue', 'light pink','#21908dff'),
    fill = c(alpha("wheat",0.3), alpha('light blue',0.3), alpha('light pink',0.3),alpha('#21908dff',0.3)),
    cex = 0.5,
    fontfamily = "sans",cat.cex = 1 ,
    cat.default.pos = "outer", cat.fontfamily = "sans",cat.col = c("wheat", 'light blue', 'light pink','#21908dff')
  )
}

MSC.markers.in <- MSC.markers
out.dir <- paste0("~/added_value/point1/Venn_diagram_v2/MSC/")

for (i in unique(MSC.markers.in$cluster)){
  cat (i)
  venn.diagram(
    x = list(MSC.markers.in$gene[MSC.markers.in$cluster==i],Pamplona.MSC.marker$gene[Pamplona.MSC.marker$cluster==i],Scadden.MSC.marker$gene[Scadden.MSC.marker$cluster==i],Aifantis.MSC.marker$gene[Aifantis.MSC.marker$cluster==i]),
    category.names = c("integrated" , "In-house", "Baryawana", "Tikhonova"),
    filename = paste0(out.dir,"/",prefix, "cluster_",i, ".png"),
    output=TRUE,
    imagetype="png" ,
    main = paste0("added_value1: Cluster ", i),
    main.cex =1,
    main.fontfamily = "Arial",
    height = 2000, 
    width = 3000, resolution = 500,compression = "lzw",
    lwd = 1,
    col=c("wheat", 'light blue', 'light pink','#21908dff'),
    fill = c(alpha("wheat",0.3), alpha('light blue',0.3), alpha('light pink',0.3),alpha('#21908dff',0.3)),
    cex = 0.5,
    fontfamily = "sans",cat.cex = 1 ,
    cat.default.pos = "outer", cat.fontfamily = "sans",cat.col = c("wheat", 'light blue', 'light pink','#21908dff')
  )
}


#### False positive, False negative

overlap_matrix_EC$False_postive <- 1-(overlap_matrix_EC$overlapped/overlap_matrix_EC$`single dataset`)
overlap_matrix_EC$False_negative <- 1- overlap_matrix_EC$shared_pt

save(overlap_matrix_EC, overlap_matrix_MSC, file = "~/added_value1_overlapped_matixs.RData")
overlap_matrix_MSC$False_postive <- 1-(overlap_matrix_MSC$overlapped/overlap_matrix_MSC$`single dataset`)
overlap_matrix_MSC$False_negative <- 1- overlap_matrix_MSC$shared_pt
overlap_matrix_MSC$False_postive[is.na(overlap_matrix_MSC$False_postive)] <- 1

p1 <- ggplot(overlap_matrix_MSC,aes_string("cluster", "dataset")) +
  theme_bw() +
  geom_tile(aes(fill = False_negative), color='white') +
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab') +
  #geom_text(aes(label = shared_pt), color='#93db69')+
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),panel.grid.major = element_blank()
        #panel.grid.major=element_line(color='#eeeeee')
  )+scale_x_discrete(position = "top") 
out.dir <- paste0("~/added_value/point1/")
prefix <- "20210404"
file_name <- paste0(out.dir,"/",prefix, "Added_value_1_MSC_False_negative.jpg")
jpeg(file=file_name, width = 7, height = 7*0.37, res= 300,units = "in", bg="transparent")
grid.arrange(p1,nrow=1)
dev.off()
file_name <- paste0(out.dir, prefix,"Added_value_1_MSC_False_negative.pdf")
pdf(file=file_name, width = 7, height = 7*0.37)
grid.arrange(p1,nrow=1)
dev.off()

p1 <- ggplot(overlap_matrix_EC,aes_string("cluster", "dataset")) +
  theme_bw() +
  geom_tile(aes(fill = False_negative), color='white') +
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab',na.value ='darkblue' ) +
  #geom_text(aes(label = shared_pt), color='#93db69')+
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),panel.grid.major = element_blank()
        #panel.grid.major=element_line(color='#eeeeee')
  )+scale_x_discrete(position = "top") 
out.dir <- paste0("~//added_value/point1/")
prefix <- "20210404"
file_name <- paste0(out.dir,"/",prefix, "Added_value_1_EC_False_negative.jpg")
jpeg(file=file_name, width = 7, height = 7*0.37, res= 300,units = "in", bg="transparent")
grid.arrange(p1,nrow=1)
dev.off()
file_name <- paste0(out.dir, prefix,"Added_value_1_EC_False_negative.pdf")
pdf(file=file_name, width = 7, height = 7*0.37)
grid.arrange(p1,nrow=1)
dev.off()

p1 <- ggplot(overlap_matrix_EC,aes_string("cluster", "dataset")) +
  theme_bw() +
  geom_tile(aes(fill = False_postive), color='white') +
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab',na.value ='darkblue' ) +
  #geom_text(aes(label = shared_pt), color='#93db69')+
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),panel.grid.major = element_blank()
        #panel.grid.major=element_line(color='#eeeeee')
  )+scale_x_discrete(position = "top") 
out.dir <- paste0("~/Figures/added_value/point1/")
prefix <- "20210404"
file_name <- paste0(out.dir,"/",prefix, "Added_value_1_EC_False_postivive.jpg")
jpeg(file=file_name, width = 7, height = 7*0.37, res= 300,units = "in", bg="transparent")
grid.arrange(p1,nrow=1)
dev.off()
file_name <- paste0(out.dir, prefix,"Added_value_1_EC_False_postivive.pdf")
pdf(file=file_name, width = 7, height = 7*0.37)
grid.arrange(p1,nrow=1)
dev.off()
