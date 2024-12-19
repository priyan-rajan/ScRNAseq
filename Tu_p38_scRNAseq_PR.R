#libraries
library(scater)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(knitr)
library(kableExtra)
library(fgsea)
library(msigdbr)
library(dplyr)
library(Seurat)
library(tidyr)
library(ggsignif)
library(scCustomize)

#Sample key
#4T1-Tu indicates 4T1 tumors treated with vehicle control
#4T1-Tu-plus-p38i indicates 4T1 tumors treated with a p38 inhibitor
#p38a-ko-H9 indicates 4T1 tumors with p38 knockout.

# colors
palette(c("dodgerblue","tomato","forestgreen","orange","orchid",
                      "darkslategray2","firebrick","lightgreen","gold","hotpink3",
                      "lightblue","deeppink","purple","darkorange","plum1"))
                      
#Load data
setwd("/projects/rpci/avbakin/prajan3")

path <- "./Bakin/GSRdata/"
files <- list.files(path)
files <- files[c(1,2,3)]
samples <- gsub("^.{12}|.{18}$","",files)

sce <- vector("list", length(samples))
names(sce) <- samples

for (i in 1:length(samples)){
  matrix_dir = paste0(path, files[i],"/outs/filtered_feature_bc_matrix/")
  cellbarcodes <- read.table(paste0(matrix_dir, "barcodes.tsv.gz"))
  genenames <- read.table(paste0(matrix_dir, "features.tsv.gz"))
  countsmat <- Matrix::readMM(paste0(matrix_dir, "matrix.mtx.gz"))
  rownames(countsmat) <- genenames[,2]
  colnames(countsmat) <- cellbarcodes[,1]
  samplename <- data.frame(samplename = rep(samples[i], ncol(countsmat)))
  sce[[i]] <- SingleCellExperiment(assays = list(counts = as.matrix(countsmat)),
                                   colData = samplename)
  colnames(sce[[i]]) <- paste(samples[i], cellbarcodes[,1], sep="_")
}
rm(list = ls()[!grepl("my|samples|sce", ls())])

#png function
pngfunc <- function(filename = "mypng.png", width=7, height = 7, units = "in",
                    res = 300, ...){
  png(filename = filename, width = width, height = height, units = units,
      res = res, ...)
}

## Per Cell QC
# Cells that meet the following criteria will be retained for further analysis
# - > 1500 counts or reads AND
# - > 300 detcted genes AND
# - < 20% mitochondrial content


sce <- lapply(sce, function(x) {
  x <- addPerCellQC(x, subsets = list(Mito = grep("^mt-", rownames(x)),
                                      Ribo = grep("^Rp[l|s]", rownames(x)),
                                      Hemo = grep("^Hb[a|b]", rownames(x))))
  x
})

sce <- lapply(sce, function(x) {
  x$qc_count <- x$sum > 1500
  x$qc_gene <- x$detected > 300
  x$qc_mito <- x$subsets_Mito_percent < 20
  
  x$qc_pass <- x$qc_count & x$qc_gene & x$qc_mito
  x
})

#QC plots
pngfunc("Plots/QC/QC.png", height = 7, width = 9)
layout(t(matrix(c(0:19),4,5)), widths = c(3,3,3,3), heights = c(1,3,3,3,1))
par(mar = c(0,0,0,0))
plot.new(); text(0.5, 0.5, "No. of reads", cex = 2)
plot.new(); text(0.5, 0.5, "No. of genes", cex = 2)
plot.new(); text(0.5, 0.5, "% Mitochondrial genes", cex=2)
for(i in samples){
  par(mar = c(0,0,0,0)) 
  plot.new(); text(0.5, 0.5, i , cex=1.5)
  par(mar = c(2,2,0,0.5))
  
  hist(log10(sce[[i]]$sum), breaks=200, main="",
       xlim=c(2.5,5.5), xlab="", ylab="")
  abline(v=log10(1500), col=2)
  
  hist(log10(sce[[i]]$detected), breaks=200, main="",
       xlim=c(1,4), xlab="", ylab="")
  abline(v=log10(300), col=2)
  
  hist(log10(sce[[i]]$subsets_Mito_percent), breaks=200, main="",
       xlim=c(0,100), xlab="", ylab="")
  abline(v=20, col=2)
}
par(mar = c(0,0,0,0))
plot.new()
plot.new(); text(0.5,0.5,"log10 total reads")
plot.new(); text(0.5,0.5,"log10 total genes")
plot.new(); text(0.5,0.5,"% Mitochondria")

invisible(dev.off())

include_graphics("Plots/QC/QC.png", dpi=300)

#QC Table
QC <- c()
for(i in samples){
  QC <- rbind(QC, as.data.frame(colData(sce[[i]])))
}

df <- data.frame(Sample = samples)
df[,c("bqc","genes", "reads", "mito", "total", "aqc", "ribo")] <- NA

for(i in df$Sample){
  x<- QC[QC$samplename==i, ]
  df$bqc[df$Sample==i] <- nrow(x)
  df$genes[df$Sample==i] <- as.integer(sum(!x$qc_gene) / nrow(x) * 100)
  df$reads[df$Sample==i] <- as.integer(sum(!x$qc_count) / nrow(x) *100)
  df$mito[df$Sample==i] <- as.integer(sum(!x$qc_mito) / nrow(x) * 100)
  
  df$total[df$Sample==i] <- as.integer(sum(!(x$qc_gene & x$qc_count & x$qc_mito)) /
                                         nrow(x) * 100)
  df$aqc[df$Sample==i] <- sum(x$qc_pass)
}

df[,1:7] %>%
  kable(format = "html", align = "c", escape = F) %>%
  kable_styling("striped")

#Combine all data into one sce object
assign("sce", do.call("cbind", sce))
saveRDS(sce, "20220302.allcells.rds")
readRDS(sce, "20220302.allcells.rds")

## Per Gene QC

cpm_mat <- as.matrix(counts(sce))
lib_size <- colSums(cpm_mat)
cpm <- t(t(cpm_mat/lib_size)*1000000)
cpm <- cpm + 1
log2cpm <- log2(cpm)

assay(sce, "cpm") <- cpm
logcounts(sce) <- log2cpm

keep <- rowVars(logcounts(sce)) > 0.5

#Per Gene QC plots
pngfunc("Plots/QC/pergeneQC.png", width = 5, height=4)
layout(t(matrix(c(0:5),2,3)), widths = c(1,4), heights = c(1,4,2))
par(mar=c(0,0,0,0)) 
plot.new();text(0.5,0.5, "Mean counts vs Variance of counts- All samples", cex=1.5)
plot.new(); text(0.5,0.5,"Variance of Gene counts (log2CPM)", srt = 90, cex=1)
par(mar = c(2,2,3,1))
plot(rowMeans(logcounts(sce)), rowVars(logcounts(sce)), xlab = "",
     ylab = "", 
     col = ifelse(rowVars(logcounts(sce))>0.5, "orange", "blue"), pch=19, cex=0.5)

plot.new();
plot.new(); text(0.5, 0.5, "Mean Gene counts (log2CPM)")

invisible(dev.off())

include_graphics("Plots/QC/pergeneQC.png", dpi=300)

#Per gene QC table
df_gene <- data.frame(Sample=samples)
df[,c("BQC","Filtered", "AQC")] <- NA
df$BQC <- nrow(sce)
df$Filtered <- as.integer(sum(!keep) / nrow(sce) * 100)
df$AQC <- sum(keep)

#Filter sce object based on QC criteria
#Per cell QC filter
sce <- sce[ ,sce$qc_pass]

#Based on per gene QC 
sce <- sce[keep, ]

## *t*- distributed stochastic neighborhood embedding 
# Assessing distribution of data points present in the samples.
#Determining tSNE dimensions
set.seed(100)
sce <- runPCA(sce)
sce <- runTSNE(sce, dimred = "PCA", perplexity = 40)
sce$TSNEx <- reducedDims(sce)$TSNE[,1]
sce$TSNEy <- reducedDims(sce)$TSNE[,2]

pngfunc(filename = "Plots/QC/tSNE.png", width= 9, height = 7)
layout(matrix(c(1,4,2,5,3,6),2,3), widths = c(3,3,3))
par(mar = c(1,1,2,1))
plot.new()
text(0.5,0.5, "All samples combined", cex=2)
plot(sce$TSNEx, sce$TSNEy, pch = 19, cex = 0.2, xlim=c(-40,40),col = factor(sce$samplename), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = "Cell types in all samples combined")
plot.new()
legend("center", legend = unique(sce$samplename), col = palette(), bty="n", ncol = 1, pch = 19, 
       title = "Samples", cex = 2)
for(i in unique(sce$samplename)){
  a <- ifelse(sce$samplename == i, 1, 0)
  b <- order(a)
  c <- ifelse(sce$samplename == i, which(unique(sce$samplename)==i), "gray")
  par(mar = c(1,1,2,1))
  plot(sce$TSNEx[b], sce$TSNEy[b], pch = 19, cex = 0.2, xlim=c(-40,40),
       col = c[b], 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", main = i)
}

invisible(dev.off())

include_graphics("Plots/QC/tSNE.png", dpi = 300)


## Cell type annotation

#cell type annotation using reference datasets
library(celldex)
ref1 <- ImmGenData()
library(SingleR)
pred1 <- SingleR(test=sce, ref=ref1, labels=ref1$label.main)

#assigning cells with their cell types
sce$immgenlabel <- pred1$labels

#heatmap to see check for unique mapping to cell types.
pngfunc("Plots/QC/scoreheatmap.png", width = 11, height = 7)
plotScoreHeatmap(pred1)
invisible(dev.off())
include_graphics("Plots/QC/scoreheatmap.png", dpi = 300)

#TSNE plots

#cell types identified overall
pngfunc("Plots/QC/cell_types.png", width = 7, height =7)
layout(matrix(c(1,2,0,3),2,2), widths = c(4,2), heights = c(2,4))

plot.new(); text(0.5, 0.5, "Cell types identified overall", cex = 2)
par(mar=c(1,1,1,0))
plot(sce$TSNEx, sce$TSNEy, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     col = factor(sce$immgenlabel), pch = 19, cex = 0.2)
plot.new(); legend("center", legend = levels(factor(sce$immgenlabel)), ncol = 2, cex = 1, title = "Cell types",
                   pch = 19, col = palette(), bty = "n")
invisible(dev.off())
include_graphics("Plots/QC/cell_types.png", dpi = 300)

#cell types identified in each sample
pngfunc("Plots/QC/cell_types_each.png", width = 6, height = 9)
layout(matrix(c(1,3,5,2,4,6), 3,2), widths = c(3,2), heights = c(3,3,3))
par(mar = c(0.5,1,2,1))
for(i in unique(sce$samplename)) {
  s <- sce[ ,sce$samplename==i]
  plot(s$TSNEx, s$TSNEy, xlab="",ylab="",xaxt="n",yaxt="n",
       main = i, col = factor(s$immgenlabel), pch = 19, cex=0.2)
  plot.new(); legend("center", legend = levels(factor(s$immgenlabel)), ncol = 2, cex = 1,
                     title = "Cell Types",
                     pch=19, col= palette())
}
invisible(dev.off())
include_graphics("Plots/QC/cell_types_each.png", dpi = 300)

#Create tables with cell type percentages
mylist <- vector("list", length(samples))
names(mylist) <- samples
for (i in 1:length(samples)){
  df <- table(sce[ ,sce$samplename==samples[i]]$immgenlabel)
  df <- as.data.frame(df)
  df$Percent <- df$Freq/sum(df$Freq) * 100
  df$Sample <- rep(samples[i], nrow(df))
  names(df)[1] <- "Cell_Type"
  mylist[[i]] <- df
}

#Combine 3 lists onto one data frame
DF <- Reduce(function(...) merge(..., all=TRUE), mylist)

#Barplot of cell type percentages
pngfunc("Plots/QC/Cell_types_Percent.png", width = 11, height  = 6)
ggplot(DF, aes(fill=Sample, x=Cell_Type, y=Percent)) +
  geom_bar(position="dodge", stat="identity")+
  ggtitle("Proportions of cell types") +
  theme_bw() +
  scale_y_continuous(limits=c(0,100))

invisible(dev.off())
include_graphics("Plots/QC/Cell_types_Percent.png")

### Subset sce to retain cell types present more than 1%

#Analysis of cell type percentages revealed that the following cell types are more than 1% abundant. These will be further analysed while the rest will be filtered out

types <- c("DC", "T cells", "Monocytes", "Macrophages", "Neutrophils",
           "NKT", "NK cells", "ILC")
sce <- sce[ ,sce$immgenlabel %in% types]

#Create tSNE plots of subsetted sce object
pngfunc("Plots/QC/cell_types_after_subsetting.png", width = 6, height = 9)
layout(matrix(c(1,3,5,2,4,6), 3,2), widths = c(3,2), heights = c(3,3,3))

par(mar = c(0.5,1,2,1))
for(i in unique(sce$samplename)) {
  s <- sce[ ,sce$samplename==i]
  plot(s$TSNEx, s$TSNEy, xlab="",ylab="",xaxt="n",yaxt="n",
       main = i, col = factor(s$immgenlabel), pch = 19, cex=0.2)
  plot.new(); legend("center", legend = levels(factor(s$immgenlabel)), ncol = 2, cex = 1,
                     title = "Cell Types",
                     pch=19, col= palette())
}
invisible(dev.off())
include_graphics("Plots/QC/cell_types_after_subsetting.png", dpi = 300)

#Create tables with cell type percentages
mylist <- vector("list", length(samples))
names(mylist) <- samples
for (i in 1:length(samples)){
  df <- table(sce[ ,sce$samplename==samples[i]]$immgenlabel)
  df <- as.data.frame(df)
  df$Percent <- df$Freq/sum(df$Freq) * 100
  df$Sample <- rep(samples[i], nrow(df))
  names(df)[1] <- "Cell_Type"
  mylist[[i]] <- df
}

#Combine 3 lists onto one data frame
DF <- Reduce(function(...) merge(..., all=TRUE), mylist)

#Barplot of cell type percentages
pngfunc("Plots/QC/Cell_types_subsetted_Percent.png", width = 11, height  = 6)
ggplot(DF, aes(fill=Sample, x=Cell_Type, y=Percent)) +
  geom_bar(position="dodge", stat="identity")+
  ggtitle("Proportions of cell types") +
  theme_bw() +
  scale_y_continuous(limits=c(0,100))

invisible(dev.off())
include_graphics("Plots/QC/Cell_types_subsetted_Percent.png")

#After this, I performed further analysis on each cell type separately
#Here I have shown analysis performed for macrophages. Similar analysis was performed for other cell types.
sce <- readRDS("20220310_subsetted.rds")
sce.mac <- sce[ ,sce$immgenlabel=="Macrophages"]
#convert to seurat
mac.seurat <- as.Seurat(sce.mac, counts = "counts", data = "logcounts")
mac.seurat
Idents(mac.seurat) <- "samplename"

#perform dge
diff1 <- FindMarkers(mac.seurat, ident.1 = "4T1-Tu-plus-p38i", ident.2 = "4T1-Tu", assay = "originalexp", slot = "data", test.use = "MAST", min.pct = 0.05, logfc.threshold = 0)
write.table(diff1, "Macrophage analysis/i_vs_v.tsv", sep = "\t")
dge1 <- diff1[diff1$p_val<0.05, ]

diff2 <- FindMarkers(mac.seurat, ident.1 = "p38a-ko-H9", ident.2 = "4T1-Tu", assay = "originalexp", slot = "data", test.use = "MAST", min.pct = 0.05, logfc.threshold = 0)
write.table(diff, "Macrophage analysis/ko_vs_v.tsv", sep = "\t")
dge2 <- diff2[diff2$p_val<0.05, ]
max(dge2$p_val)

#from Mouse Gene Ontology database (MSIGDB-Broad Institute)- Biological processes is loaded 
m_gene_sets <-msigdbr(species = "mouse", category = "C5", subcategory = "GO:BP")
msigdb_list_go_bp <- split(x=m_gene_sets$gene_symbol, f=m_gene_sets$gs_name)

#ranks is a named numeric vector
ranks <- dge1[,2]
names(ranks) <- rownames(dge1)
ranks

#arrange in decreasing order of logFC values
ranks <- sort(ranks, decreasing = TRUE)

#running fGSEA
fgseares1 <- fgsea(msigdb_list_go_bp, stats=ranks)
fgseares1 <- fgseares1 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval<0.05)

df_res1 = data.frame(lapply(fgseares1, as.character), stringsAsFactors=FALSE)

write.table(df_res1, "Macrophage analysis/gsea_i_vs_v.tsv", sep = "\t")

# barplot
pngfunc("Plots/QC/Macrophage/i_vs_v_v4.png", width = 6.5, height = 3)
ggplot(fgseares1 %>% as_tibble() %>% arrange(desc(NES))
       %>% filter(pval<0.05) %>% head(n=7), 
       aes(reorder(pathway,NES), NES)) +
  geom_col(aes(fill=NES)) +
  coord_flip()+ labs(x="Pathway", y="Normalized Enrichment Score",
                     title="Pathways enriched in p38i") + 
  theme_minimal()
invisible(dev.off())

#heatmap markers

genes <- c("Tgfb1", "Vegfa", "Osm", "Ptgs2", "Il1b", "Lpl", "Cd36", "Cxcr4", "Cxcl2", "Dusp1", "Ano6", "Sphk1", "Slpi", "Tnfrsf9", "Cd244a","Mapkapk3", "Il7r", "Il1r2", "Rnase2a")


mean_mat <- AverageExpression(mac.seurat, group.by = "samplename",
                              layer = "data", features = genes)
mean_mat
mat <- mean_mat[[1]]
mat <- t(scale(t(mat)))
rownames(mat) <- genes
colnames(mat) <- c("Veh", "p38i", "p38ko")
mat
cluster_anno <- unique(mac.seurat$samplename)

quantile(mat, c(0.1, 0.95))

f1 = colorRamp2(c(-1.1,0,1.2), c("blue","white","red"), space = "RGB")
pngfunc("Plots/QC/Macrophage/Heatmap_Mac_v3.png", width = 2.5, height  = 3.5)
Heatmap(mat, 
        name="Expression", 
        cluster_rows = FALSE, 
        col = f1, 
        column_split = cluster_anno,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize=10),
        row_names_gp = gpar(fontsize=9),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4
)
invisible(dev.off())

#Violin plots

pngfunc("Plots/QC/Macrophage/with p38ko/Sphk1.png", width = 3.5, height = 3)

VlnPlot(mac.seurat, features = "Sphk1", pt.size = 0, group.by = "samplename", layer = "data") + labs(title = "Sphk1", x = "", y = "") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(face = "italic", hjust = 0.5)) + 
  geom_signif(comparisons = list(c("4T1-Tu-plus-p38i", "4T1-Tu"), c("p38a-ko-H9", "4T1-Tu")), test = "t.test" , test.args = list(alternative = "two.sided", var.equal = FALSE, paired = FALSE), map_signif_level = TRUE, textsize = 7, y_position = c(18, 22), tip_length = 0) + ylim(NA, 25) + scale_fill_manual(values = c("gray71", "deeppink", "forestgreen"))

