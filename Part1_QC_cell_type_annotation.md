QC and Cell type annotation
================
Priyanka Rajan
2/4/2025

## Introduction

Triple negative breast cancer has very poor prognosis (Leon-Ferre,
2022). While immunetherapy with immune checkpoint inhibitors has
improved survival outcomes, the benefit is limited to a small subset of
patients (10-15%) (Abdou Y., et al, 2022). Thus, there is a need to
develop effective therapeutic approaches for the treatment of TNBC.
Tumor microenvironment of TNBC consists of immune suppressive myeloid
populations such as TAMs and MDSCs which hinder anti-tumor CD8+ T cell
responses . Moreover, tumor induces the generation and recruitment of
TAMs and MDSCs. Identification of mechanisms that mediate crosstalk
between tumor and immune cells will provide actionable therapeutic
targets. Our group has identified that the activity of p38 alpha (p38a)
MAP kinase in the tumor cells promotes the generation and recruitment of
these populations (Rajan P., et al, 2025, under review at Can Imm Res).
In this study, we evaluated the effect of p38a blockade using
pharmacological inhibition and genetic depletion on the tumor immune
landscape in the 4T1 mouse mammary carcinoma model. Live CD45+ immune
populations were flow sorted from 4T1 tumors treated with vehicle
(4T1-Tu), p38 inhibitor LY2228820 (4T1-Tu-plus-p38i) or 4T1 tumors with
p38a knockout (p38a-ko-H9) and single cell RNA-seq was performed using
the 10X Chromium 3’ RACE technology. CellRanger was used to align reads
to mouse genome build mm10 and generate read-count matrices. In this
section (Part 1), quality control is performed to filter out low quality
cells and cells are visualized in 2D using dimensionality reduction
algorithms.

Setup:

``` r
knitr::opts_chunk$set(echo = TRUE, cache = FALSE)

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
library(RSpectra)
library(dplyr)
library(celldex)
library(SingleR)

# colors
palette(c("dodgerblue","tomato","forestgreen","orange","orchid",
          "darkslategray2","firebrick","lightgreen","gold","hotpink3",
          "lightblue","deeppink","purple","darkorange","plum1"))
```

Reading input data:

``` r
#Load data
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
```

``` r
#png function
pngfunc <- function(filename = "mypng.png", width=7, height = 7, units = "in",
                    res = 300, ...){
  png(filename = filename, width = width, height = height, units = units,
      res = res, ...)
}
```

## Per Cell QC

Cells that meet the following criteria will be retained for further
analysis \* \> 1500 read counts \* \> 300 detected genes \* \< 20%
mitochondrial content

``` r
#Using the addPerCellQC function from scater package to compute quality contorl metrics.
sce <- lapply(sce, function(x) {
  x <- addPerCellQC(x, subsets = list(Mito = grep("^mt-", rownames(x)),
                                      Ribo = grep("^Rp[l|s]", rownames(x)),
                                      Hemo = grep("^Hb[a|b]", rownames(x))))
  x
})

#Determining which cells pass the criteria
sce <- lapply(sce, function(x) {
  x$qc_count <- x$sum > 1500
  x$qc_gene <- x$detected > 300
  x$qc_mito <- x$subsets_Mito_percent < 20
  
  x$qc_pass <- x$qc_count & x$qc_gene & x$qc_mito
  x
  })

#Generating QC plots
pngfunc("Plots/QC2/percellQC.png", height = 7, width = 9)
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

include_graphics("Plots/QC2/percellQC.png", dpi=300)
```

![](Plots/QC2/percellQC.png)<!-- -->

``` r
#Generating QC Table
QC <- c()
for(i in samples){
  QC <- rbind(QC, as.data.frame(colData(sce[[i]])))
}

df <- data.frame(Sample = samples)
df[,c("BQC","genes", "reads", "mito", "total", "AQC", "ribo")] <- NA

for(i in df$Sample){
  x<- QC[QC$samplename==i, ]
  df$BQC[df$Sample==i] <- nrow(x)
  df$genes[df$Sample==i] <- as.integer(sum(!x$qc_gene) / nrow(x) * 100)
  df$reads[df$Sample==i] <- as.integer(sum(!x$qc_count) / nrow(x) *100)
  df$mito[df$Sample==i] <- as.integer(sum(!x$qc_mito) / nrow(x) * 100)
  
  df$total[df$Sample==i] <- as.integer(sum(!(x$qc_gene & x$qc_count & x$qc_mito)) /
    nrow(x) * 100)
  df$AQC[df$Sample==i] <- sum(x$qc_pass)
 }

df[,1:7] %>%
  kable(format = "html", align = "c", escape = F) %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:center;">
Sample
</th>
<th style="text-align:center;">
BQC
</th>
<th style="text-align:center;">
genes
</th>
<th style="text-align:center;">
reads
</th>
<th style="text-align:center;">
mito
</th>
<th style="text-align:center;">
total
</th>
<th style="text-align:center;">
AQC
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center;">
4T1-Tu
</td>
<td style="text-align:center;">
4588
</td>
<td style="text-align:center;">
1
</td>
<td style="text-align:center;">
15
</td>
<td style="text-align:center;">
3
</td>
<td style="text-align:center;">
16
</td>
<td style="text-align:center;">
3819
</td>
</tr>
<tr>
<td style="text-align:center;">
4T1-Tu-plus-p38i
</td>
<td style="text-align:center;">
6343
</td>
<td style="text-align:center;">
0
</td>
<td style="text-align:center;">
24
</td>
<td style="text-align:center;">
1
</td>
<td style="text-align:center;">
24
</td>
<td style="text-align:center;">
4791
</td>
</tr>
<tr>
<td style="text-align:center;">
p38a-ko-H9
</td>
<td style="text-align:center;">
5212
</td>
<td style="text-align:center;">
4
</td>
<td style="text-align:center;">
23
</td>
<td style="text-align:center;">
8
</td>
<td style="text-align:center;">
25
</td>
<td style="text-align:center;">
3857
</td>
</tr>
</tbody>
</table>

``` r
#Combine all data into one sce object
assign("sce", do.call("cbind", sce))
```

## Per Gene QC

The counts are normalized to library size and converted to counts per
million units in log2 scale. Ideally we would like to analyze genes that
vary upon p38 inhibition by p38i or p38KO as those are the most
interesting changes. Hence, genes that vary little between samples are
removed.

``` r
#Normalizing counts and converting to log2 counts per million
cpm_mat <- as.matrix(counts(sce))
lib_size <- colSums(cpm_mat)
cpm <- t(t(cpm_mat/lib_size)*1000000)
cpm <- cpm + 1
log2cpm <- log2(cpm)
  
assay(sce, "cpm") <- cpm
logcounts(sce) <- log2cpm
  
#Calculating which genes show variance greater than 0.5
keep <- rowVars(logcounts(sce)) > 0.5

#Generating per Gene QC plots
pngfunc("Plots/QC2/pergeneQC.png", width = 5, height=4)
layout(t(matrix(c(0:5),2,3)), widths = c(1,4), heights = c(1,4,2))
par(mar=c(0,0,0,0)) 
plot.new();text(0.5,0.5, "Mean vs Variance plot", cex=1.5)
plot.new(); text(0.5,0.5,"Variance of Gene counts (log2CPM)", srt = 90, cex=1.2)
par(mar = c(2,2,3,1))
plot(rowMeans(logcounts(sce)), rowVars(logcounts(sce)), xlab = "",
       ylab = "", 
     col = ifelse(rowVars(logcounts(sce))>0.5, "orange", "blue"), pch=19, cex=0.5)

plot.new();
plot.new(); text(0.5, 0.5, "Mean Gene counts (log2CPM)", cex = 1.2)

invisible(dev.off())

include_graphics("Plots/QC2/pergeneQC.png", dpi=300)
```

![](Plots/QC2/pergeneQC.png)<!-- -->

``` r
#Generating per gene QC table
df_gene <- data.frame(Sample=samples)
df_gene[,c("BQC","P_Filtered", "AQC")] <- NA
df_gene$BQC <- nrow(sce)
df_gene$P_Filtered <- as.integer(sum(!keep) / nrow(sce) * 100)
df_gene$AQC <- sum(keep)

df_gene[,1:4] %>%
  kable(format = "html", align = "c", escape = F) %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:center;">
Sample
</th>
<th style="text-align:center;">
BQC
</th>
<th style="text-align:center;">
P_Filtered
</th>
<th style="text-align:center;">
AQC
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:center;">
4T1-Tu
</td>
<td style="text-align:center;">
32285
</td>
<td style="text-align:center;">
65
</td>
<td style="text-align:center;">
11142
</td>
</tr>
<tr>
<td style="text-align:center;">
4T1-Tu-plus-p38i
</td>
<td style="text-align:center;">
32285
</td>
<td style="text-align:center;">
65
</td>
<td style="text-align:center;">
11142
</td>
</tr>
<tr>
<td style="text-align:center;">
p38a-ko-H9
</td>
<td style="text-align:center;">
32285
</td>
<td style="text-align:center;">
65
</td>
<td style="text-align:center;">
11142
</td>
</tr>
</tbody>
</table>

``` r
#Filtering sce object based on QC criteria
sce <- sce[keep, ]
```

## Visualizing in 2D space

Principal component analysis (PCA) is a commonly used method to compact
data and assess distribution of data points. The runPCA() function
determines PCs using the first 500 genes that show the highest
variation. *t*- distributed stochastic neighborhood embedding is used
for more complex single cell datasets. Here, I ran t-SNE using the
pre-existing PCA results as input to the tSNE algorithm. This is useful
as it improves speed by using a low-rank approximation of the expression
matrix; and reduces random noise, by focusing on the major factors of
variation.

``` r
#Determining tSNE dimensions
set.seed(100)
sce <- runPCA(sce)
sce <- runTSNE(sce, dimred = "PCA", perplexity = 40)

sce$TSNEx <- reducedDims(sce)$TSNE[,1]
sce$TSNEy <- reducedDims(sce)$TSNE[,2]

pngfunc(filename = "Plots/QC2/tSNE.png", width= 9, height = 7)
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

include_graphics("Plots/QC2/tSNE.png", dpi = 300)
```

![](Plots/QC2/tSNE.png)<!-- -->

## Annotation using ImmGen data

The immune populations present in the samples were identified using
SingleR. SingleR uses a reference dataset containing known labels. It
then labels each cell based on similarity to the reference. This is done
by comparing the expression of a set of marker genes and using spearman
correlation to assign a label to each cell. Immunological genome
database containing microarray based gene expression data for multiple
mouse immune cell populations was used as a reference. This was accessed
through the celldex package.

``` r
#cell type annotation using reference dataset from ImmGen db

ref1 <- celldex::ImmGenData()
pred1 <- SingleR(test=sce, ref=ref1, labels=ref1$label.main)

#assigning cells with their cell types
sce$immgenlabel <- pred1$labels

#heatmap to which label received the highest assignment score in each cell
pngfunc("Plots/QC2/scoreheatmap.png", width = 11, height = 7)
plotScoreHeatmap(pred1)
invisible(dev.off())
include_graphics("Plots/QC2/scoreheatmap.png", dpi = 300)
```

![](Plots/QC2/scoreheatmap.png)<!-- -->

``` r
#TSNE plots

#cell types identified overall
pngfunc("Plots/QC2/cell_types.png", width = 7, height =7)
layout(matrix(c(1,2,0,3),2,2), widths = c(4,2), heights = c(2,4))

plot.new(); text(0.5, 0.5, "Cell types identified overall", cex = 2)
par(mar=c(1,1,1,0))
plot(sce$TSNEx, sce$TSNEy, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     col = factor(sce$immgenlabel), pch = 19, cex = 0.2)
plot.new(); legend("center", legend = levels(factor(sce$immgenlabel)), ncol = 2, cex = 1, title = "Cell types",
       pch = 19, col = palette(), bty = "n")
invisible(dev.off())
include_graphics("Plots/QC2/cell_types.png", dpi = 300)
```

![](Plots/QC2/cell_types.png)<!-- -->

``` r
#cell types identified in each sample
pngfunc("Plots/QC2/cell_types_each.png", width = 6, height = 9)
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
include_graphics("Plots/QC2/cell_types_each.png", dpi = 300)
```

![](Plots/QC2/cell_types_each.png)<!-- -->

``` r
#Create tables with cell type percentages
mylist <- vector("list", length(samples))
names(mylist) <- samples
for (i in 1:length(samples)){
  df_type <- table(sce[ ,sce$samplename==samples[i]]$immgenlabel)
  df_type <- as.data.frame(df_type)
  df_type$Percent <- df_type$Freq/sum(df_type$Freq) * 100
  df_type$Sample <- rep(samples[i], nrow(df_type))
  names(df_type)[1] <- "Cell_Type"
  mylist[[i]] <- df_type
  }

#Combine 3 lists onto one data frame
DF_type <- Reduce(function(...) merge(..., all=TRUE), mylist)

#Barplot of cell type percentages

pngfunc("Plots/QC2/Cell_types_Percent.png", width = 11, height  = 6)
ggplot(DF_type, aes(fill=Sample, x=Cell_Type, y=Percent)) +
  geom_bar(position="dodge", stat="identity")+
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ggtitle("Proportions of cell types") +
  theme_bw() +
  scale_y_continuous(limits=c(0,100))

invisible(dev.off())
include_graphics("Plots/QC2/Cell_types_Percent.png")
```

![](Plots/QC2/Cell_types_Percent.png)<!-- -->

### Subset sce to retain cell types present more than 1%

Analysis of cell type percentages revealed that the following cell types
are more than 1% abundant. These will be further analysed while the rest
will be filtered out

``` r
types <- c("DC", "T cells", "Monocytes", "Macrophages", "Neutrophils",
            "NKT", "NK cells", "ILC")
sce <- sce[ ,sce$immgenlabel %in% types]

#Create tSNE plots of subsetted sce object
pngfunc("Plots/QC2/cell_types_after_subsetting.png", width = 6, height = 9)
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
include_graphics("Plots/QC2/cell_types_after_subsetting.png", dpi = 300)
```

![](Plots/QC2/cell_types_after_subsetting.png)<!-- -->

``` r
#Create tables with cell type percentages
mylist <- vector("list", length(samples))
names(mylist) <- samples
for (i in 1:length(samples)){
  df_subset <- table(sce[ ,sce$samplename==samples[i]]$immgenlabel)
  df_subset <- as.data.frame(df_subset)
  df_subset$Percent <- df_subset$Freq/sum(df_subset$Freq) * 100
  df_subset$Sample <- rep(samples[i], nrow(df_subset))
  names(df_subset)[1] <- "Cell_Type"
  mylist[[i]] <- df_subset
}

#Combine 3 lists onto one data frame
DF_subset <- Reduce(function(...) merge(..., all=TRUE), mylist)

#Barplot of cell type percentages
pngfunc("Plots/QC2/Cell_types_subsetted_Percent.png", width = 11, height  = 6)
ggplot(DF_subset, aes(fill=Sample, x=Cell_Type, y=Percent)) +
  geom_bar(position="dodge", stat="identity")+
  ggtitle("Proportions of cell types") +
  theme_bw() +
  scale_y_continuous(limits=c(0,100))

invisible(dev.off())
include_graphics("Plots/QC2/Cell_types_subsetted_Percent.png")
```

![](Plots/QC2/Cell_types_subsetted_Percent.png)<!-- -->

### Plotting cell types based on markers

Feature plots are a convenient way to represent distinct cell types as
well as their distribution relative to other cell types in tSNE space.
Here, two gene markers (features) representing each of the immune
populations is shown in red and blue. Combination of red and blue gives
purple, highlight the specific cell population.

``` r
#Feature plots
rgb <- vector("list", length = 8)
names(rgb) <- types
for(i in types){
  if(i == "T cells"){
    cd8 <- as.numeric(logcounts(sce["Cd8a",]))
    cd3d <- as.numeric(logcounts(sce["Cd3d",]))
    cd8 <- cd8/max(cd8)
    cd3d <- cd3d/max(cd3d)
    rgb[[i]] <- rgb(cd8,0,cd3d)
  } else if(i == "Monocytes"){
    ccr2 <- as.numeric(logcounts(sce["Ccr2",]))
    ly6c2 <- as.numeric(logcounts(sce["Ly6c2",]))
    ccr2 <- ccr2/max(ccr2)
    ly6c2 <- ly6c2/max(ly6c2)
    rgb[[i]] <- rgb(ccr2,0,ly6c2)
  } else if(i == "Neutrophils"){
    csf3r <- as.numeric(logcounts(sce["Csf3r",]))
    g0s2 <- as.numeric(logcounts(sce["G0s2",]))
    csf3r <- csf3r/max(csf3r)
    g0s2 <- g0s2/max(g0s2)
    rgb[[i]] <- rgb(csf3r,0,g0s2)
  } else if(i == "DC"){
    p2ry10 <- as.numeric(logcounts(sce["P2ry10",]))
    itgax <- as.numeric(logcounts(sce["Itgax",]))
    p2ry10 <- p2ry10/max(p2ry10)
    itgax <- itgax/max(itgax)
    rgb[[i]] <- rgb(p2ry10,0,itgax)
  } else if(i == "Macrophages"){
    arg1 <- as.numeric(logcounts(sce["Arg1",]))
    c1qa <- as.numeric(logcounts(sce["C1qa",]))
    arg1 <- arg1/max(arg1)
    c1qa <- c1qa/max(c1qa)
    rgb[[i]]<- rgb(arg1,0,c1qa)
  } else if(i == "NK cells"){
    nkg7 <- as.numeric(logcounts(sce["Nkg7",]))
    gzmb <- as.numeric(logcounts(sce["Gzmb",]))
    nkg7 <- nkg7/max(nkg7)
    gzmb <- gzmb/max(gzmb)
    rgb[[i]] <- rgb(nkg7,0,gzmb)
  } else if(i== "NKT"){
    icos <- as.numeric(logcounts(sce["Icos",]))
    lck <- as.numeric(logcounts(sce["Lck",]))
    icos <- icos/max(icos)
    lck <- lck/max(lck)
    rgb[[i]] <- rgb(icos,0,lck)
  } else if(i == "ILC"){
    tbx21 <- as.numeric(logcounts(sce["Tbx21",]))
    ncr1 <- as.numeric(logcounts(sce["Ncr1",]))
    tbx21 <- tbx21/max(tbx21)
    ncr1 <- ncr1/max(ncr1)
    rgb[[i]]<- rgb(tbx21,0,ncr1)
  } else {
    print("cell type not found")
  }
}

pngfunc("Plots/QC2/feature.png", width = 9, height = 5)
par(mfrow=c(2,4))
par(mar = c(1,1,2,1))
for(i in types){
  plot(sce$TSNEx, sce$TSNEy, col = rgb[[i]], pch = 19, cex=0.5, 
       xlab= "", ylab="", xaxt="n", yaxt="n", main= i)
  }

invisible(dev.off())
include_graphics("Plots/QC2/feature.png", dpi = 300)
```

![](Plots/QC2/feature.png)<!-- -->

## Session info

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 24.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /cvmfs/soft.ccr.buffalo.edu/versions/2023.01/easybuild/software/avx512/Compiler/gcc/11.2.0/flexiblas/3.0.4/lib/libflexiblas.so.3.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C         LC_TIME=C           
    ##  [4] LC_COLLATE=C         LC_MONETARY=C        LC_MESSAGES=C       
    ##  [7] LC_PAPER=C           LC_NAME=C            LC_ADDRESS=C        
    ## [10] LC_TELEPHONE=C       LC_MEASUREMENT=C     LC_IDENTIFICATION=C 
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] SingleR_2.0.0               celldex_1.8.0              
    ##  [3] dplyr_1.1.4                 RSpectra_0.16-2            
    ##  [5] kableExtra_1.4.0            knitr_1.38                 
    ##  [7] Matrix_1.6-5                circlize_0.4.14            
    ##  [9] ComplexHeatmap_2.14.0       scran_1.26.2               
    ## [11] scater_1.26.1               ggplot2_3.5.1              
    ## [13] scuttle_1.8.4               SingleCellExperiment_1.20.1
    ## [15] SummarizedExperiment_1.28.0 Biobase_2.58.0             
    ## [17] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
    ## [19] IRanges_2.32.0              S4Vectors_0.36.2           
    ## [21] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
    ## [23] matrixStats_0.62.0         
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] AnnotationHub_3.6.0           BiocFileCache_2.6.1          
    ##   [3] systemfonts_1.1.0             igraph_1.3.1                 
    ##   [5] BiocParallel_1.32.6           digest_0.6.29                
    ##   [7] foreach_1.5.2                 htmltools_0.5.2              
    ##   [9] viridis_0.6.2                 fansi_1.0.6                  
    ##  [11] magrittr_2.0.3                memoise_2.0.1                
    ##  [13] ScaledMatrix_1.6.0            cluster_2.1.3                
    ##  [15] doParallel_1.0.17             limma_3.54.2                 
    ##  [17] Biostrings_2.66.0             svglite_2.1.3                
    ##  [19] colorspace_2.1-1              blob_1.2.3                   
    ##  [21] rappdirs_0.3.3                ggrepel_0.9.1                
    ##  [23] xfun_0.30                     crayon_1.5.1                 
    ##  [25] RCurl_1.98-1.6                iterators_1.0.14             
    ##  [27] glue_1.8.0                    gtable_0.3.6                 
    ##  [29] zlibbioc_1.44.0               XVector_0.38.0               
    ##  [31] GetoptLong_1.0.5              DelayedArray_0.24.0          
    ##  [33] BiocSingular_1.14.0           shape_1.4.6                  
    ##  [35] scales_1.3.0                  pheatmap_1.0.12              
    ##  [37] DBI_1.1.2                     edgeR_3.40.2                 
    ##  [39] Rcpp_1.0.8.3                  viridisLite_0.4.2            
    ##  [41] xtable_1.8-4                  clue_0.3-60                  
    ##  [43] dqrng_0.4.1                   bit_4.0.4                    
    ##  [45] rsvd_1.0.5                    metapod_1.6.0                
    ##  [47] httr_1.4.2                    RColorBrewer_1.1-3           
    ##  [49] ellipsis_0.3.2                pkgconfig_2.0.3              
    ##  [51] farver_2.1.2                  dbplyr_2.1.1                 
    ##  [53] locfit_1.5-9.5                utf8_1.2.4                   
    ##  [55] tidyselect_1.2.1              labeling_0.4.3               
    ##  [57] rlang_1.1.4                   later_1.3.0                  
    ##  [59] AnnotationDbi_1.60.2          munsell_0.5.1                
    ##  [61] BiocVersion_3.16.0            tools_4.2.0                  
    ##  [63] cachem_1.0.6                  cli_3.6.3                    
    ##  [65] generics_0.1.3                RSQLite_2.2.12               
    ##  [67] ExperimentHub_2.6.0           evaluate_0.15                
    ##  [69] stringr_1.5.1                 fastmap_1.1.0                
    ##  [71] yaml_2.3.5                    bit64_4.0.5                  
    ##  [73] purrr_1.0.2                   KEGGREST_1.38.0              
    ##  [75] sparseMatrixStats_1.10.0      mime_0.12                    
    ##  [77] xml2_1.3.3                    compiler_4.2.0               
    ##  [79] rstudioapi_0.13               beeswarm_0.4.0               
    ##  [81] filelock_1.0.3                curl_4.3.2                   
    ##  [83] png_0.1-7                     interactiveDisplayBase_1.36.0
    ##  [85] tibble_3.2.1                  statmod_1.4.36               
    ##  [87] stringi_1.7.6                 highr_0.9                    
    ##  [89] lattice_0.20-45               bluster_1.8.0                
    ##  [91] vctrs_0.6.5                   pillar_1.9.0                 
    ##  [93] lifecycle_1.0.4               BiocManager_1.30.25          
    ##  [95] GlobalOptions_0.1.2           BiocNeighbors_1.16.0         
    ##  [97] bitops_1.0-7                  irlba_2.3.5                  
    ##  [99] httpuv_1.6.5                  R6_2.5.1                     
    ## [101] promises_1.2.0.1              gridExtra_2.3                
    ## [103] vipor_0.4.5                   codetools_0.2-18             
    ## [105] assertthat_0.2.1              rjson_0.2.21                 
    ## [107] withr_3.0.2                   GenomeInfoDbData_1.2.9       
    ## [109] parallel_4.2.0                beachmat_2.14.2              
    ## [111] rmarkdown_2.14                DelayedMatrixStats_1.20.0    
    ## [113] Rtsne_0.16                    shiny_1.7.1                  
    ## [115] ggbeeswarm_0.6.0
