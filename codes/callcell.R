#!/usr/bin/env Rscript
library(magrittr)
library(Signac)
library(GenomicRanges)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(DoubletFinder)


# RNA, ATAC and PET barcode selection uses the configured QC thresholds.
# Singlet filtering, RNA/ATAC WNN and plots follow the author's callcell.R/function.R.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("Usage: Rscript callcell.R SAMPLE RNA_MATRIX ATAC_MATRIX FRAGMENTS PET_STATS OUTPUT_DIR GENOME")
}
sampleid <- args[1]
rna.path <- normalizePath(args[2], mustWork = TRUE)
atac.path <- normalizePath(args[3], mustWork = TRUE)
fragment.path <- normalizePath(args[4], mustWork = TRUE)
pet.path <- normalizePath(args[5], mustWork = TRUE)
output.path <- args[6]
refgenome <- args[7]
if (!refgenome %in% c("mm10", "hg38")) stop("GENOME must be mm10 or hg38")
if (!file.exists(paste0(fragment.path, ".tbi"))) stop("Fragment .tbi index is required")
dir.create(output.path, recursive = TRUE, showWarnings = FALSE)
output.path <- normalizePath(output.path, mustWork = TRUE)
script.path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
source(file.path(dirname(normalizePath(script.path)), "function.R"))
setwd(output.path)
unlink("SUCCESS.txt")
dir.create("data", showWarnings = FALSE)
set.seed(as.integer(Sys.getenv("CHAIR_SEED", "1234")))
min.rna <- as.numeric(Sys.getenv("CHAIR_MIN_RNA", "400"))
min.atac <- as.numeric(Sys.getenv("CHAIR_MIN_ATAC", "4000"))
min.pet <- as.numeric(Sys.getenv("CHAIR_MIN_PET", "300"))
max.mt <- as.numeric(Sys.getenv("CHAIR_MAX_MT", "20"))
max.ribo <- as.numeric(Sys.getenv("CHAIR_MAX_RIBO", "50"))
doublet.rate <- as.numeric(Sys.getenv("CHAIR_DOUBLET_RATE", "0.054"))
if (any(!is.finite(c(min.rna, min.atac, min.pet, max.mt, max.ribo, doublet.rate)))) stop("Invalid QC configuration")
expected.md5 <- Sys.getenv("CHAIR_FRAGMENTS_MD5")
if (nzchar(expected.md5) && unname(tools::md5sum(fragment.path)) != expected.md5) stop("Fragment MD5 differs from the selected input")
expected.index.md5 <- Sys.getenv("CHAIR_FRAGMENTS_INDEX_MD5")
if (nzchar(expected.index.md5) && unname(tools::md5sum(paste0(fragment.path, ".tbi"))) != expected.index.md5) stop("Fragment index MD5 differs from the selected input")
input.paths <- c(RNA = rna.path, ATAC = atac.path, fragments = fragment.path, fragment_index = paste0(fragment.path, ".tbi"), PET = pet.path)
write.table(data.frame(input = names(input.paths), path = unname(input.paths), bytes = file.info(input.paths)$size), "data/input_paths.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read.matrix <- function(path, feature.type) {
  counts <- if (dir.exists(path)) Read10X(path) else Read10X_h5(path)
  if (is.list(counts)) counts <- counts[[feature.type]]
  if (is.null(counts)) stop("Requested feature type is missing from ", path)
  counts
}
count_rna <- read.matrix(rna.path, "Gene Expression")
count_atac <- read.matrix(atac.path, "Peaks")
count_rna <- count_rna[, Matrix::colSums(count_rna) > 0, drop = FALSE]
atac.ranges <- StringToGRanges(rownames(count_atac), sep = c(":", "-"))
standard.atac <- as.character(seqnames(atac.ranges)) %in% GenomeInfoDb::standardChromosomes(atac.ranges)
count_atac <- count_atac[standard.atac, , drop = FALSE]
count_atac <- count_atac[, Matrix::colSums(count_atac) > 0, drop = FALSE]
pet <- read.table(pet.path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
required.pet <- c("PET_total", "PET_intra", "PET_inter", "PET_intrale1000", "PET_intragt1000", "PET_intragt8000", "PET_intragt10000", "PET_intragt20000")
if (!all(required.pet %in% colnames(pet))) stop("PET statistics are missing required columns: ", paste(setdiff(required.pet, colnames(pet)), collapse = ", "))
if (anyDuplicated(rownames(pet))) stop("PET statistics contain duplicate barcodes")
pet.cells <- rownames(pet)[is.finite(pet$PET_intragt1000) & pet$PET_intragt1000 >= min.pet]
cells <- intersect(intersect(colnames(count_rna), colnames(count_atac)), pet.cells)
if (!length(cells)) stop("No RNA/ATAC/PET shared barcodes pass PET filtering")
writeLines(cells, "data/shared_barcodes.txt")
pet <- pet[cells, , drop = FALSE]
cat("RNA barcodes", ncol(count_rna), "ATAC barcodes", ncol(count_atac), "PET-selected shared barcodes", length(cells), "\n")

if (refgenome == "mm10") {
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  annotation.db <- EnsDb.Mmusculus.v79
  mtgene <- "^mt-"
  ribogene <- "^Rp[sl]"
} else {
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  annotation.db <- EnsDb.Hsapiens.v86
  mtgene <- "^MT-"
  ribogene <- "^RP[SL]"
}
annotations <- GetGRangesFromEnsDb(ensdb = annotation.db)
seqlevels(annotations) <- gsub("chrMT", "chrM", paste0("chr", seqlevels(annotations)))
genome(annotations) <- refgenome

rna.dims <- as.integer(Sys.getenv("CHAIR_RNA_DIMS", "30"))
atac.dims <- as.integer(Sys.getenv("CHAIR_ATAC_DIMS", "30"))
rna.resolution <- as.numeric(Sys.getenv("CHAIR_RNA_RESOLUTION", "0.8"))
atac.resolution <- as.numeric(Sys.getenv("CHAIR_ATAC_RESOLUTION", "0.8"))
wnn.resolution <- as.numeric(Sys.getenv("CHAIR_WNN_RESOLUTION", "0.5"))
rna.cluster <- paste0("RNA_snn_res.", rna.resolution)
atac.cluster <- paste0("ATAC_snn_res.", atac.resolution)
wnn.cluster <- paste0("rnaatac.snn_res.", wnn.resolution)
if (anyNA(c(rna.dims, atac.dims)) || rna.dims < 2 || atac.dims < 3 || max(rna.dims, atac.dims) > 50) stop("RNA/ATAC dimensions must be within the supplied reductions (2 to 50)")
if (any(!is.finite(c(rna.resolution, atac.resolution, wnn.resolution))) || min(rna.resolution, atac.resolution, wnn.resolution) <= 0) stop("Invalid clustering resolution")

## RNA
data <- CreateSeuratObject(counts = count_rna[, cells], project = sampleid, assay = "RNA", min.cells = 5)

data[['percent.mt']] <- PercentageFeatureSet(object = data, features = rownames(data[grepl(mtgene, rownames(data)), ]))
data[['percent.ribo']] <- PercentageFeatureSet(object = data, features = rownames(data[grepl(ribogene, rownames(data)), ]))
data <- data[!grepl(mtgene, rownames(data)) & !grepl(ribogene, rownames(data)), ]

## ATAC
data[["ATAC"]] <- CreateChromatinAssay(counts = count_atac[, cells], 
                                       sep = c(":", "-"), 
                                       genome = refgenome, 
                                       fragments = fragment.path, 
                                       min.cells = 5, 
                                       annotation = annotations)

DefaultAssay(data) <- "ATAC"
data <- TSSEnrichment(data)
data <- NucleosomeSignal(data, n = NULL)

blacklist.name <- paste0("blacklist_", refgenome)
utils::data(list = blacklist.name, package = "Signac", envir = environment())
data$blacklist_fraction <- FractionCountsInRegion(object = data, assay = "ATAC", regions = get(blacklist.name))

data <- AddMetaData(data, metadata = pet)

data@meta.data[is.na(data@meta.data)] <- 0

DefaultAssay(data) <- "RNA"
saveRDS(data, file = paste0("data/", sampleid, ".all.rds"))

dataflt <- subset(data, subset = nCount_RNA >= min.rna & nCount_ATAC >= min.atac & percent.mt < max.mt & percent.ribo < max.ribo)
cat(sprintf("There are %d high quality cell.\n", ncol(dataflt)))

# QC plots use all QC-passed cells before doublet removal.
qc.features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "PET_total", "PET_inter", "PET_intragt1000")
pdf(paste0(sampleid, ".raw.count.pdf"), width = 15, height = 9)
print(VlnPlot(dataflt, features = qc.features, pt.size = 0, ncol = 4))
print(VlnPlot(dataflt, features = qc.features, pt.size = 0, ncol = 4, log = TRUE))
p1 <- FeatureScatter(dataflt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 <- FeatureScatter(dataflt, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(dataflt, feature1 = "nCount_RNA", feature2 = "nCount_ATAC")
p4 <- FeatureScatter(dataflt, feature1 = "nCount_RNA", feature2 = "PET_intragt1000")
print(p1 + p2 + p3 + p4)
dev.off()

dataflt <- runRNA(dataflt)
dataflt <- runATAC(dataflt)

##### Remove doublet cells #####
dataflt <- finddoublet(dataflt, DoubletRate = doublet.rate)

saveRDS(dataflt, file = paste0("data/", sampleid, ".QC.rds"))
metadata <- as.data.frame(dataflt@meta.data) %>% dplyr::mutate(cell = rownames(.))
metadata <- metadata[, c("cell", setdiff(colnames(metadata), "cell")), drop = FALSE]
write.table(metadata, file = paste0("data/", sampleid, ".metadata.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

cluster_theme <- theme(axis.text = element_blank(), axis.ticks = element_blank())

pdf(paste0(sampleid, ".raw.RNAcluster.pdf"), width = 15)
print(DimPlot(dataflt, reduction = "rnatsne", group.by = c("RNADF.class", rna.cluster), label = TRUE, label.size = 5, repel = TRUE))
print(DimPlot(dataflt, reduction = "rnaumap", group.by = c("RNADF.class", rna.cluster), label = TRUE, label.size = 5, repel = TRUE))
dev.off()

pdf(paste0(sampleid, ".raw.ATACcluster.pdf"), width = 15)
print(DimPlot(dataflt, reduction = "atactsne", group.by = c("RNADF.class", atac.cluster), label = TRUE, label.size = 5, repel = TRUE))
print(DimPlot(dataflt, reduction = "atacumap", group.by = c("RNADF.class", atac.cluster), label = TRUE, label.size = 5, repel = TRUE))
dev.off()

pdf(paste0(sampleid, ".doublet.pdf"), width = 13)
p1 <- DimPlot(object = dataflt, reduction = "rnaumap", group.by = "RNADF.class", label = T, label.size = 5) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  ggtitle("RNA UMAP") + 
  NoLegend() + 
  cluster_theme
p2 <- DimPlot(object = dataflt, reduction = "rnatsne", group.by = "RNADF.class", label = T, label.size = 5) + 
  labs(x = "TSNE 1", y = "TSNE 2") + 
  ggtitle("RNA t-SNE") + 
  NoLegend() + 
  cluster_theme
print(p1 + p2)
dev.off()

dataflt <- subset(dataflt, subset = RNADF.class == "Singlet")
cat(sprintf("There are %d singlet cell.\n", ncol(dataflt)))

##### Cluster using singlet cells #####
meta_remove <- c("RNA.weight", "ATAC.weight", rna.cluster, atac.cluster, wnn.cluster, "seurat_clusters")
dataflt@meta.data <- dataflt@meta.data[, !colnames(dataflt@meta.data) %in% meta_remove]
dataflt@reductions <- list()

dataflt <- runRNA(dataflt)
dataflt <- runATAC(dataflt)
dataflt <- runWNN(dataflt, reduction.list = list("rnapca", "ataclsi"), dims.list = list(seq_len(rna.dims), seq.int(2, atac.dims)), reduction.name = "rnaatac")

metadata <- clusterdeeper(data = dataflt, idents = rna.cluster, outpath = "deepcluster", resultname = "rnacluster1")
dataflt <- AddMetaData(object = dataflt, metadata = metadata[, c("cluster")], col.name = c("rnacluster2"))

# s.genes <- gprofiler2::gorth(cc.genes.updated.2019$s.genes)$ortholog_name
# g2m.genes <- gprofiler2::gorth(cc.genes.updated.2019$g2m.genes)$ortholog_name
s.genes <- cc.genes.updated.2019$s.genes
if (refgenome == "mm10") s.genes <- stringr::str_to_title(s.genes)
g2m.genes <- cc.genes.updated.2019$g2m.genes
if (refgenome == "mm10") g2m.genes <- stringr::str_to_title(g2m.genes)

DefaultAssay(dataflt) <- "RNA"
DefaultAssay(dataflt) <- "RNA"
dataflt <- CellCycleScoring(dataflt, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

saveRDS(dataflt, file = paste0("data/", sampleid, ".final.rds"))

DefaultAssay(dataflt) <- "RNA"
idents <- c(rnacluster1 = rna.cluster, rnacluster2 = "rnacluster2")
findmarker(data = dataflt, idents = idents, outpath = "data", express = T)
plotcluster(data = dataflt, idents = idents, filename = paste0("deepcluster/", sampleid, ".rnacluster.pdf"), width = 24, plottheme = cluster_theme)

pdf(paste0(sampleid,".final.cluster1.pdf"), width = 10)
p1 <- DimPlot(object = dataflt, reduction = "rnaumap", group.by = rna.cluster, label = T, label.size = 5) + labs(x = "UMAP 1", y = "UMAP 2") + ggtitle("RNA UMAP") + NoLegend() + cluster_theme
p2 <- DimPlot(object = dataflt, reduction = "rnatsne", group.by = rna.cluster, label = T, label.size = 5) + labs(x = "t-SNE 1", y = "t-SNE 2") + ggtitle("RNA t-SNE") + NoLegend() + cluster_theme
p3 <- DimPlot(object = dataflt, reduction = "atacumap", group.by = atac.cluster, label = T, label.size = 5) + labs(x = "UMAP 1", y = "UMAP 2") + ggtitle("ATAC UMAP") + NoLegend() + cluster_theme
p4 <- DimPlot(object = dataflt, reduction = "atactsne", group.by = atac.cluster, label = T, label.size = 5) + labs(x = "t-SNE 1", y = "t-SNE 2") + ggtitle("ATAC t-SNE") + NoLegend() + cluster_theme
p5 <- DimPlot(object = dataflt, reduction = "rnaatacumap", group.by = wnn.cluster, label = T, label.size = 5) + labs(x = "UMAP 1", y = "UMAP 2") + ggtitle("RNA+ATAC UMAP") + NoLegend() + cluster_theme
print(p1 + p3 + p5 + p2 + p4)
# print(p1 + p2 + p3 + p4)
dev.off()

pdf(paste0(sampleid,".final.cluster2.pdf"), width = 15, height = 10)
p1 <- DimPlot(object = dataflt, reduction = "rnaumap", group.by = "rnacluster2", label = T, label.size = 5) + ggtitle("RNA UMAP") + labs(x = "UMAP 1", y = "UMAP 2") + NoLegend() + cluster_theme
p2 <- DimPlot(object = dataflt, reduction = "rnatsne", group.by = "rnacluster2", label = T, label.size = 5) + ggtitle("RNA t-SNE") + labs(x = "t-SNE 1", y = "t-SNE 2") + NoLegend() + cluster_theme
p3 <- DimPlot(object = dataflt, reduction = "atacumap", group.by = "rnacluster2", label = T, label.size = 5) + ggtitle("ATAC UMAP") + labs(x = "UMAP 1", y = "UMAP 2") + NoLegend() + cluster_theme
p4 <- DimPlot(object = dataflt, reduction = "atactsne", group.by = "rnacluster2", label = T, label.size = 5) + ggtitle("ATAC t-SNE") + labs(x = "t-SNE 1", y = "t-SNE 2") + NoLegend() + cluster_theme
p5 <- DimPlot(object = dataflt, reduction = "rnaatacumap", group.by = "rnacluster2", label = T, label.size = 5) + ggtitle("RNA+ATAC UMAP") + labs(x = "UMAP 1", y = "UMAP 2") + NoLegend() + cluster_theme
print(p1 + p3 + p5 + p2 + p4)
# print(p1 + p2 + p3 + p4)
dev.off()

stat_features <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "percent.mt", "percent.ribo")
pdf(paste0(sampleid, ".final.count.pre.origident.pdf"), width = 6)
p <- VlnPlot(dataflt, features = stat_features, group.by = "orig.ident", pt.size = 0, ncol = 2)
print(p)
p <- VlnPlot(dataflt, features = stat_features, group.by = "orig.ident", pt.size = 0, ncol = 2, log = T)
print(p)
dev.off()

pdf(paste0(sampleid, ".final.count.rnacluster.pdf"), width = 10)
p <- VlnPlot(dataflt, features = stat_features, group.by = rna.cluster, pt.size = 0, ncol = 3)
print(p)
p <- VlnPlot(dataflt, features = stat_features, group.by = rna.cluster, pt.size = 0, ncol = 3, log = T)
print(p)
dev.off()

# title is pearson correlation between the two features
pdf(paste0(sampleid, ".final.featurescatter.pdf"), width = 18)
p1 <- FeatureScatter(dataflt, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident") + NoLegend()
p2 <- FeatureScatter(dataflt, feature1 = "nCount_RNA", feature2 = "percent.ribo", group.by = "orig.ident") + NoLegend()
p3 <- FeatureScatter(dataflt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") + NoLegend()
print(p1 + p2 + p3)
p1 <- FeatureScatter(dataflt, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC", group.by = "orig.ident") + NoLegend()
p2 <- FeatureScatter(dataflt, feature1 = "nCount_ATAC", feature2 = "nCount_RNA", group.by = "orig.ident") + NoLegend()
p3 <- FeatureScatter(dataflt, feature1 = "nCount_ATAC", feature2 = "nFeature_RNA", group.by = "orig.ident") + NoLegend()
print(p1 + p2 + p3)
dev.off()

pdf(paste0(sampleid, ".final.densityscatter.pdf"))
DensityScatter(dataflt, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = T, quantiles = T) + ylim(c(0, 10))
dev.off()

marker1 <- read.table("data/degene.rnacluster1.txt", header = T, sep = "\t", check.names = F)
marker1 <- marker1[order(-marker1$avg_log2FC, -marker1$pct.1), ] %>% 
  group_by(cluster) %>% 
  slice_head(n = 5) %>% 
  as.data.frame()

marker2 <- read.table("data/degene.rnacluster2.txt", header = T, sep = "\t", check.names = F)
marker2 <- marker2[order(-marker2$avg_log2FC, -marker2$pct.1), ] %>% 
  group_by(cluster) %>% 
  slice_head(n = 5) %>% 
  as.data.frame()

mytheme <- theme(legend.position = "none",
                 axis.text = element_text(colour = "black"))

dataflt <- ScaleData(dataflt)
pdf(paste0(sampleid, ".degene.heatmap.rnacluster1.pdf"), height = 10, width = 14)
DoHeatmap(object = dataflt, features = marker1$gene, group.by = rna.cluster) + mytheme
dev.off()

pdf(paste0(sampleid, ".degene.heatmap.rnacluster2.pdf"), height = 20, width = 28)
DoHeatmap(object = dataflt, features = marker2$gene, group.by = "rnacluster2") + mytheme
dev.off()



dataflt <- readRDS(paste0("data/", sampleid, ".final.rds"))

stat_features <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "PET_total", "PET_intragt1000", "PET_inter")

mytheme <- theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        panel.grid = element_blank(), 
        axis.text = element_text(colour = "black"))

plotdata <- as.data.frame(dataflt@meta.data)

pdf(paste0(sampleid, ".final.count.origident.pdf"), height = 6, width = 4)
for (i in stat_features) {
  ymax <- max(plotdata[[i]])
  p <- ggplot(plotdata, aes(x = orig.ident, y = !!as.name(i))) +
    geom_violin(aes(fill = orig.ident)) +
    geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
    stat_summary(fun = "mean", geom = "point", shape = 8, size = 2.5, color = "black") +
    stat_summary(fun.data = function(x){
      median_value <- median(x)
      formatted_value <- if (all(x %% 1 == 0)) median_value else signif(median_value, 3)
      data.frame(y = ymax*0.9, label = paste("M=", formatted_value))
    }, geom = "text", size = 4) +
    stat_summary(fun.data = function(x) data.frame(y = ymax*0.8, label = paste("N=", length(x))), geom = "text", size = 4) +
    labs(title = "", x = "", y = i) +
    mytheme
  print(p)
}
dev.off()

stat_features <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC", "PET_total", "PET_intragt1000", "PET_inter")
for (i in stat_features) {
  tmp <- as.character(summary(dataflt[[i]]))
  tmpsummary <- data.frame(statistic = gsub(" *: *.*$", "", tmp), result = format(as.numeric(gsub(" *", "", gsub("^.*: *", "", tmp))), big.mark = ",", scientific = F))
  colnames(tmpsummary) <- c("statistic", i)
  if (i == "nCount_RNA") {
    finalsummary <- tmpsummary
  } else {
    finalsummary <- merge(finalsummary, tmpsummary, by = "statistic")
  }
}

finalsummary$statistic <- factor(finalsummary$statistic, levels = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
finalsummary <- finalsummary[order(finalsummary$statistic), ]
finalsummary$LibraryID <- paste0(sampleid, " (n=", ncol(dataflt), ")")
write.table(finalsummary[, c("LibraryID", "statistic", stat_features)], file = paste0("data/", sampleid, ".summary.tsv"), quote = F, sep = "\t", row.names = F)

writeLines(colnames(dataflt), "data/final_barcodes.txt")
write.table(data.frame(cells = ncol(dataflt), median_nucleosome_signal = median(dataflt$nucleosome_signal)), "data/nucleosome_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "data/sessionInfo.txt")
writeLines("complete", "SUCCESS.txt")
