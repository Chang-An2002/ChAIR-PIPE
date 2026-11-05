library(Seurat)
library(Signac)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)

checkdir <- function(dirpath) {
  if (!dir.exists(dirpath)) {
    dir.create(dirpath, recursive = T)
    print(paste0(dirpath, " is created!"))
  }
}

# Use the same selected, positive-variance scaled genes in Seurat v4 and v5.
# This is the matrix RunPCA route already validated in the author's revision run.
runSelectedPCA <- function(data, assay = "RNA", features = VariableFeatures(data[[assay]]), npcs = 50, reduction.name = "pca", reduction.key = "PC_") {
  scaled <- if (inherits(data[[assay]], "Assay5")) {
    LayerData(data, assay = assay, layer = "scale.data")
  } else {
    GetAssayData(data, assay = assay, slot = "scale.data")
  }
  features <- intersect(features, rownames(scaled))
  mat <- scaled[features, , drop = FALSE]
  variances <- matrixStats::rowVars(as.matrix(mat))
  features <- features[is.finite(variances) & variances > 0]
  npcs <- min(npcs, length(features) - 1L, ncol(data) - 1L)
  if (npcs < 2L) stop("Too few variable genes or cells for PCA")
  reduction <- Seurat::RunPCA(mat[features, , drop = FALSE], assay = assay, npcs = npcs, reduction.key = reduction.key, seed.use = 42, rev.pca = FALSE, weight.by.var = TRUE, approx = TRUE)
  stopifnot(identical(rownames(Loadings(reduction)), features))
  data[[reduction.name]] <- reduction
  command.name <- paste0("RunPCA.", assay)
  data[[command.name]] <- new("SeuratCommand", name = command.name, time.stamp = Sys.time(), assay.used = assay, call.string = "Seurat::RunPCA(selected_scaled_variable_gene_matrix)", params = list(assay = assay, features = features, npcs = npcs, rev.pca = FALSE, weight.by.var = TRUE, approx = TRUE, seed.use = 42, reduction.name = reduction.name, reduction.key = reduction.key))
  message("PCA used ", length(features), " selected genes")
  data
}

runRNA <- function(data, findcluster = T) {
  DefaultAssay(data) <- "RNA"
  data <- NormalizeData(object = data)
  data <- FindVariableFeatures(object = data, nfeatures = 3000)
  data <- ScaleData(object = data, features = rownames(data))
  data <- runSelectedPCA(data, assay = "RNA", features = VariableFeatures(data), reduction.name = "rnapca", reduction.key = "rnaPCA_")
  data <- RunUMAP(data, reduction = "rnapca", dims = seq_len(rna.dims), reduction.name = "rnaumap", reduction.key = "rnaUMAP_")
  data <- RunTSNE(data, reduction = "rnapca", dims = seq_len(rna.dims), reduction.name = "rnatsne", reduction.key = "rnaTSNE_", check_duplicates = F)
  if (findcluster == T) {
    data <- FindNeighbors(object = data, reduction = "rnapca", dims = seq_len(rna.dims))
    data <- FindClusters(data, resolution = rna.resolution)
  }
  
  return(data)
}

runATAC <- function(data, findcluster = T) {
  DefaultAssay(data) <- "ATAC"
  data <- RunTFIDF(data)
  data <- FindTopFeatures(data, assay = "ATAC", min.cutoff = 'q0')
  data <- RunSVD(data, assay = "ATAC", reduction.name = "ataclsi", reduction.key = "atacLSI_")
  data <- RunUMAP(object = data, reduction = 'ataclsi', dims = seq.int(2, atac.dims), reduction.name = "atacumap", reduction.key = "atacUMAP_")
  data <- RunTSNE(object = data, reduction = 'ataclsi', dims = seq.int(2, atac.dims), reduction.name = "atactsne", reduction.key = "atacSNE_", check_duplicates = F)
  if (findcluster == T) {
    data <- FindNeighbors(object = data, reduction = 'ataclsi', dims = seq.int(2, atac.dims))
    data <- FindClusters(object = data, resolution = atac.resolution, verbose = F, algorithm = 3)
  }
  
  return(data)
}

runPET <- function(data, findcluster = T) {
  DefaultAssay(data) <- "PET"
  data <- FindTopFeatures(data, assay = "PET", min.cutoff = "q1")
  data <- RunSVD(data, assay = "PET", reduction.name = "petlsi", reduction.key = "petLSI_")
  data <- RunUMAP(data, reduction = "petlsi", dims = 1:20, reduction.name = "petumap", reduction.key = "petUMAP_")
  data <- RunTSNE(data, reduction = "petlsi", dims = 1:20, reduction.name = "pettsne", reduction.key = "petTSNE_", check_duplicates = F)
  if (findcluster == T) {
    data <- FindNeighbors(data, reduction = 'petlsi', dims = 1:20, k.param = 30)
    data <- FindClusters(data, resolution = 0.5, algorithm = 3)
  }
  
  return(data)
}

runWNN <- function(data, reduction.list, dims.list, reduction.name) {
  data <- FindMultiModalNeighbors(data,
                                  reduction.list = reduction.list,
                                  dims.list = dims.list,
                                  knn.graph.name = paste0(reduction.name, ".knn"),
                                  snn.graph.name = paste0(reduction.name, ".snn"),
                                  weighted.nn.name = paste0(reduction.name, ".wnn"))
  data <- FindClusters(data, graph.name = paste0(reduction.name, ".snn"), algorithm = 3, resolution = wnn.resolution)
  data <- RunUMAP(data, nn.name = paste0(reduction.name, ".wnn"), reduction.name = paste0(reduction.name, "umap"), reduction.key = paste0(reduction.name, "UMAP_"))
  
  return(data)
}

baseals <- function(data, type, findcluster = T, integrate = T) {
  dims <- list(rna = 1:30, atac = 2:30, pet = 1:20)
  reductionnames <- c(rna = "rnapca", atac = "ataclsi", pet = "petlsi")
  
  for (type in types){
    if (type == "rna")
      data <- runRNA(data = data, findcluster = findcluster)
    else if (type == "atac")
      data <- runATAC(data = data, findcluster = findcluster)
    else if (type == "pet")
      data <- runPET(data = data, findcluster = findcluster)
    
    ## Without runRNA, runATAC, runPET functions
    # origreducname <- as.character(reductionnames[type])
    # DefaultAssay(data) <- toupper(type)
    # 
    # if (type == "rna") {
    #   data <- NormalizeData(object = data)
    #   data <- FindVariableFeatures(object = data, assay = toupper(type), nfeatures = 3000)
    #   data <- ScaleData(object = data, assay = toupper(type), features = rownames(data))
    #   data <- RunPCA(object = data, assay = toupper(type), features = VariableFeatures(data), reduction.name = origreducname, reduction.key = "rnaPCA_")
    # }
    # else if (type %in% c("atac", "pet")) {
    #   data <- RunTFIDF(data)
    #   data <- FindTopFeatures(object = data, assay = toupper(type), min.cutoff = 'q0')
    #   data <- RunSVD(object = data, assay = toupper(type), reduction.name = origreducname, reduction.key = paste0(type, "LSI_"))
    # }
    # 
    # data <- RunUMAP(object = data, reduction = origreducname, dims = dims, reduction.name = paste0(type, "umap"), reduction.key = paste0(type, "UMAP_"))
    # data <- RunTSNE(object = data, reduction = origreducname, dims = dims, reduction.name = paste0(type, "tsne"), reduction.key = paste0(type, "TSNE_"), check_duplicates = F)
    # 
    # if (findcluster == T) {
    #   k.param <- ifelse(type == "pet", 30, 20)
    #   resolution <- ifelse(type == "pet", 0.5, 0.8)
    #   verbose <- ifelse(type == "atac", TRUE, FALSE)
    #   algorithm <- ifelse(type == "rna", 1, 3)
    #   
    #   data <- FindNeighbors(object = data, reduction = origreducname, dims = dims, k.param = k.param)
    #   data <- FindClusters(object = data, resolution = resolution, verbose = verbose, algorithm = algorithm)
    # }
  }
  
  if (integrate == T) {
    alltypes <- c("rna", "atac", "pet")
    types <- alltypes[alltypes %in% types]
    
    combns <- unlist(lapply(2:length(types), function(n) {combn(types, n, simplify = F)}), recursive = F)
    for (i in combns) {
      data <- runWNN(data = data, reduction.name = paste0(i), reduction.list = as.vector(reductionnames[i]), dims.list = dims[i])
      
      ## Without runRNA, runATAC, runPET functions
      # reduction.name <- paste0(i)
      # reduction.list <- as.vector(reductionnames[i])
      # dims.list <- dims[i]
      # 
      # data <- runWNN(data = data, reduction.list = reduction.list, dims.list = dims.list, reduction.name = reduction.name)
      # data <- FindMultiModalNeighbors(object = data, 
      #                                 algorithm = 3, 
      #                                 resolution = 0.5, 
      #                                 reduction.list = reduction.list,
      #                                 dims.list = dims.list,
      #                                 knn.graph.name = paste0(reduction.name, ".knn"),
      #                                 snn.graph.name = paste0(reduction.name, ".snn"),
      #                                 weighted.nn.name = paste0(reduction.name, ".wnn"))
      # data <- FindClusters(data, graph.name = paste0(reduction.name, ".snn"))
      # data <- RunUMAP(data, nn.name = paste0(reduction.name, ".wnn"), reduction.name = paste0(reduction.name, "umap"), reduction.key = paste0(reduction.name, "UMAP_"))
    }
  }
}

finddoublet <- function(data, DoubletRate = 0.054) {
  # DoubletFinder reuses this ScaleData command when building artificial cells.
  data <- ScaleData(data, assay = "RNA", features = VariableFeatures(data[["RNA"]]))
  data@reductions$pca <- data@reductions$rnapca
  
  # find the best pk value
  pcs <- rna.dims
  sweep.function <- get(if ("paramSweep" %in% getNamespaceExports("DoubletFinder")) "paramSweep" else "paramSweep_v3", envir = asNamespace("DoubletFinder"))
  sweep <- sweep.function(data, PCs = 1:pcs, sct = F)
  sweep.stats <- summarizeSweep(sweep, GT = F)
  
  pdf("doubletfinder.pktest.pdf")
  bcmvn <- find.pK(sweep.stats)
  dev.off()
  
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  # detect doublet cell
  # DoubletRate <- ncol(data) * 8 * 1e-6 # if cell number is far more than 7000
  homotypic.prop <- modelHomotypic(data[[rna.cluster]][, 1])
  nExp_poi <- round(DoubletRate * ncol(data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  doublet.function <- get(if ("doubletFinder" %in% getNamespaceExports("DoubletFinder")) "doubletFinder" else "doubletFinder_v3", envir = asNamespace("DoubletFinder"))
  data <- doublet.function(data, PCs = 1:pcs, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  
  DF.name <- colnames(data@meta.data)[grepl("DF.classification", colnames(data@meta.data))]
  pANN.name <- colnames(data@meta.data)[grepl("pANN", colnames(data@meta.data))]
  data <- AddMetaData(data, data[[pANN.name]], col.name = "RNApANN")
  data <- AddMetaData(data, data[[DF.name]], col.name = "RNADF.class")
  
  data[[pANN.name]] <- NULL
  data[[DF.name]] <- NULL
  data@reductions$pca <- NULL
  
  return(data)
}

findmarker <- function(data, idents, outpath = NULL, express = F, return = F, save = T) {
  if (!is.null(outpath)) {
    if (!grepl("/$", outpath))
      outpath <- paste0(outpath, "/")
    
    if (!dir.exists(outpath))
      dir.create(outpath, recursive = T)
  }
  
  for (ident in names(idents)) {
    Idents(data) <- as.character(idents[ident])
    marker <- FindAllMarkers(object = data)
    marker <- marker[marker$avg_log2FC > 0 & marker$p_val_adj < 0.05, ]
    marker <- marker[order(-marker$avg_log2FC, marker$p_val_adj), ] %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(order = dplyr::row_number()) %>% 
      dplyr::slice_head(n = 100) %>% 
      as.data.frame()
    
    marker <- marker[, c("cluster", "gene", "order", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")] %>% mutate(ident = ident)
    
    if (express == T) {
      exp <- as.data.frame(AverageExpression(object = data, group.by = idents[ident])$RNA)
      exp <- reshape2::melt(as.matrix(exp[rownames(exp) %in% unique(marker$gene), ]))
      exp <- data.frame(cluster = gsub("g", "", exp$Var2), gene = exp$Var1, expression = exp$value)
      marker <- merge(marker, exp, by = c("cluster", "gene"), all.x = T)
    }
    
    if (save == T)
      write.table(marker, file = paste0(outpath, "degene.", ident, ".txt"), quote = F, sep = "\t", row.names = F)
    
    if (return == T) {
      if (ident == names(idents)[1])
        allmarker <- marker
      else
        allmarker <- rbind(allmarker, marker)
    }
  }
  
  if (return == T)
    return(allmarker)
}

integrate <- function(data, reference, method = "cca") {

  common_feature <- intersect(rownames(data), rownames(reference))
  integratedata <- merge(x = data, y = reference, add.cell.ids = c("data", "ref"))
  if (sum(grepl("^counts\\.", Layers(integratedata))) > 1)
    integratedata[["RNA"]] <- JoinLayers(integratedata[["RNA"]])
  integratedata <- integratedata[common_feature, ]
  
  if (method == "cca")
    integratedata[["RNA"]] <- split(integratedata[["RNA"]], f = integratedata$orig.ident)
  
  integratedata <- NormalizeData(object = integratedata)
  integratedata <- FindVariableFeatures(object = integratedata, nfeatures = 3000)
  integratedata <- ScaleData(integratedata)
  integratedata <- RunPCA(integratedata)
  integratedata <- RunUMAP(integratedata, dims = seq_len(rna.dims), reduction = "pca", reduction.name = "umap.unintegrated")
  
  if (method == "cca") {
    integratedata <- IntegrateLayers(object = integratedata, method = CCAIntegration, orig.reduction = "pca", new.reduction = "rnapca", verbose = F)
    integratedata[["RNA"]] <- JoinLayers(integratedata[["RNA"]])
  }
  else if (method == "harmony") {
    integratedata <- RunHarmony(integratedata, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmonypca", max.iter = 100)
    integratedata@reductions$rnapca <- integratedata@reductions$harmonypca
  }
  
  # 下面这两行（FindNeighbors和FindClusters是为了整合library临时解除注释的，用完了注释回去）
  # integratedata <- FindNeighbors(integratedata, reduction = "rnapca", dims = seq_len(rna.dims))
  # integratedata <- FindClusters(integratedata, resolution = 0.8)

  integratedata <- RunUMAP(integratedata, reduction = "rnapca", dims = seq_len(rna.dims), reduction.name = "rnaumap", reduction.key = "rnaUMAP_")
  integratedata <- RunTSNE(integratedata, reduction = "rnapca", dims = seq_len(rna.dims), reduction.name = "rnatsne", reduction.key = "rnaTSNE_", check_duplicates = F)
  
  return(integratedata)
}

mapmarker <- function(data = NULL, idents = NULL, degene = NULL, reference, marker, outpath, resultname, savemarker = T) {
  if (!grepl("/$", outpath))
    outpath <- paste0(outpath, "/")
  
  if (!dir.exists(outpath))
    dir.create(outpath, recursive = T)
  
  if (!is.null(idents) & !is.null(data) & is.null(degene)) {
    Idents(data) <- idents
    marker <- FindAllMarkers(object = data)
    marker <- marker[marker$avg_log2FC > 0 & marker$p_val_adj < 0.05, ]
    marker <- marker[order(-marker$avg_log2FC, marker$p_val_adj), ] %>%
      group_by(cluster) %>%
      mutate(order_origin = row_number()) %>% 
      as.data.frame()
    
    if (savemarker == T) {
      marker <- marker[, c("cluster", "gene", "order_origin", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]
      write.table(marker, file = paste0(outpath, resultname, ".degene.txt"), quote = F, sep = "\t", row.names = F)
    }
    
    exp <- as.data.frame(AverageExpression(object = data, assays = "RNA", group.by = idents)$RNA)
    exp <- reshape2::melt(as.matrix(exp[rownames(exp) %in% unique(marker$gene), ]))
    exp <- data.frame(cluster = gsub("g", "", exp$Var2), gene = exp$Var1, expression = exp$value)
    
    degene <- merge(marker, exp, by = c("cluster", "gene"), all.x = T)
  }
  
  annodata <- merge(degene, marker, by = "gene", all.x = T)
  annodata <- annodata[!is.na(annodata$order), ]
  
  for (i in sort(unique(annodata$cluster))) {
    tmp1 <- annodata[annodata$cluster == i, c("order_origin", "expression", "avg_log2FC", "gene", "group")] %>% 
      mutate(x = "target")
    colnames(tmp1)[1] <- "order"
    tmp2 <- annodata[annodata$cluster == i, c("order", "expression", "avg_log2FC", "gene", "group")] %>% 
      mutate(x = "reference")
    plotdata <- rbind(tmp1, tmp2)
    n <- length(unique(plotdata$group))
    p <- ggplot(plotdata, aes(x = x, y = order)) +
      geom_line(aes(group = gene), color = 'grey40', lwd = 0.5)+
      geom_point(aes(color = expression, size = avg_log2FC))+
      theme_bw() +
      theme(axis.title.x = element_blank())
    if (length(unique(plotdata$group)) > 1)
      p <- p + facet_wrap(. ~ group, ncol = 15) + labs(title = paste0("cluster ", i, "; ", reference))
    else
      p <- p + labs(title = paste0("cluster: ", i, "; group: ", unique(plotdata$group), "; ", reference))
    ggsave(paste0(outpath, resultname, ".marker1." ,i, ".pdf"), p, width = n*5+1, height = 6, limitsize = F)
  }
  system(paste0("pdfunite ", paste0(outpath, resultname, ".marker1.*.pdf", collapse = " "), " ", paste0(outpath, resultname, ".marker1.pdf")))
  system(paste0("rm ", paste0(outpath, resultname, ".marker1.*.pdf", collapse = " ")))
  for (i in sort(unique(annodata$cluster))) {
    k <- 1
    plotdata <- annodata[annodata$cluster == i, ]
    n <- length(unique(plotdata$group))
    for (j in c("avg_log2FC", "pct.1", "expression")) {
      p <- ggplot(annodata[annodata$cluster == i, ], aes(x = cluster, y = !!as.name(j))) +
        geom_violin(aes(fill = group)) +
        geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
        theme_bw() +
        theme(legend.position = "none", 
              axis.title.x = element_blank())
      if (length(unique(plotdata$group)) > 1)
        p <- p + facet_wrap(. ~ group, ncol = 15) + labs(title = paste0("cluster ", i, "; ", reference))
      else
        p <- p + labs(title = paste0("cluster: ", i, "; group: ", unique(plotdata$group), "; ", reference))
      assign(paste0("p", k), p)
      k <- k + 1
    }
    ggsave(paste0(outpath, resultname, ".marker2.", i, ".pdf"), plot = wrap_plots(p1, p2, p3, ncol = 1, nrow = 3), width = n*5, height = 18, limitsize = F)
  }
  system(paste0("pdfunite ", paste0(outpath, resultname, ".marker2.*.pdf", collapse = " "), " ", paste0(outpath, resultname, ".marker2.pdf")))
  system(paste0("rm ", paste0(outpath, resultname, ".marker2.*.pdf", collapse = " ")))
}

clusterdeeper <- function(data, idents, outpath = ".", resultname, marker = NULL, savedata = F, plotcluster = T, plotmarker = F) {
  if (!grepl("/$", outpath))
    outpath <- paste0(outpath, "/")
  
  if (!dir.exists(outpath))
    dir.create(outpath, recursive = T)
  
  metadata <- as.data.frame(data@meta.data)
  metadata$cluster <- "unclassified"
  
  for (i in unique(data[[idents]][, 1])) {
    cells <- rownames(metadata[metadata[idents] == i, ])
    ncells <- length(cells)
    dims <- 30
    if (ncells < 50) {
      npcs <- ncells-1
      dims <- ncells-1
    }
    else
      npcs <- 50
    
    tmpdata <- subset(data, cells = cells)
    
    DefaultAssay(tmpdata) <- "RNA"
    tmpdata <- NormalizeData(object = tmpdata)
    tmpdata <- FindVariableFeatures(object = tmpdata, nfeatures = 3000)
    tmpdata <- ScaleData(object = tmpdata, features = rownames(tmpdata))
    tmpdata <- runSelectedPCA(tmpdata, npcs = npcs, features = VariableFeatures(tmpdata))
    tmpdata <- FindNeighbors(object = tmpdata, dims = 1:dims)
    tmpdata <- FindClusters(tmpdata, resolution = 0.8)
    
    tmpmeta <- as.data.frame(tmpdata@meta.data)
    metadata[rownames(tmpmeta), ]$cluster <- paste0(i, "-", tmpmeta$RNA_snn_res.0.8)
    
    if (length(unique(tmpdata@meta.data$seurat_clusters)) < 2)
      next
    
    tmpdata <- RunUMAP(tmpdata, dims = 1:dims, reduction.name = "pcaumap", reduction.key = "pcaUMAP_")
    tmpdata <- RunTSNE(tmpdata, dims = 1:dims, reduction.name = "pcatsne", perplexity = 5, reduction.key = "pcaTSNE_", check_duplicates = F)
    
    if (savedata == T)
      saveRDS(tmpdata, file = paste0(outpath, resultname, ".", i, ".rds"))
    
    if (plotcluster == T) {
      pdf(paste0(outpath, resultname, ".", i, ".cluster.pdf"), width = 15)
      p1 <- DimPlot(object = tmpdata, pt.size = 0.1, reduction = "pca", group.by = "RNA_snn_res.0.8", label = T, label.size = 5) + ggtitle("PCA") + NoLegend()
      p3 <- DimPlot(object = tmpdata, pt.size = 0.1, reduction = "pcaumap", group.by = "RNA_snn_res.0.8", label = T, label.size = 5) + ggtitle("UMAP") + NoLegend()
      p2 <- DimPlot(object = tmpdata, pt.size = 0.1, reduction = "pcatsne", group.by = "RNA_snn_res.0.8", label = T, label.size = 5) + ggtitle("TSNE") + NoLegend()
      print((p1 + p2 + p3))
      dev.off()
    }
  }
  
  return(metadata)
}

plotcluster <- function(data, idents, filename, plottheme = NULL, height = 8, width = 8) {
  if (!dir.exists(dirname(filename)))
    dir.create(dirname(filename))
  
  if (is.null(plottheme))
    plottheme <- theme(axis.text = element_blank(), axis.ticks = element_blank())
  
  pdf(file = filename, height = height, width = width)
  for (i in idents) {
    p1 <- DimPlot(data, reduction = "rnapca", label = T, label.size = 5, raster = F, group.by = i) + ggtitle("PCA") + NoLegend() + plottheme
    p2 <- DimPlot(data, reduction = "rnaumap", label = T, label.size = 5, raster = F, group.by = i) + ggtitle("UMAP") + NoLegend() + plottheme
    p3 <- DimPlot(data, reduction = "rnatsne", label = T, label.size = 5, raster = F, group.by = i) + ggtitle("t-SNE") + plottheme
    print(p1 + p2 + p3)
  }
  dev.off()
}

plotref <- function(data, idents, outpath = ".") {
  if (!grepl("/$", outpath))
    outpath <- paste0(outpath, "/")
  
  if (!dir.exists(outpath))
    dir.create(outpath, recursive = T)
  
  for (ident in names(idents)) {
    Idents(data) <- as.character(idents[ident])
    pdf(paste0(outpath, "reference.count.", ident, ".pdf"), width = 10)
    p <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, ncol = 2)
    print(p)
    dev.off()
    
    pdf(paste0(outpath, "reference.featurescatter.", ident, ".pdf"), width = 10)
    p <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(p)
    dev.off()
    
    cluster_theme <- theme(axis.text = element_blank(), axis.ticks = element_blank())
    
    pdf(paste0(outpath, "reference.cluster.", ident, ".pdf"), width = 24)
    p1 <- DimPlot(object = data, pt.size = 0.1, reduction = "rnapca", group.by = idents[ident]) + ggtitle("PCA") + NoLegend() + cluster_theme
    p2 <- DimPlot(object = data, pt.size = 0.1, reduction = "rnaumap", group.by = idents[ident]) + ggtitle("UMAP") + NoLegend() + cluster_theme
    p3 <- DimPlot(object = data, pt.size = 0.1, reduction = "rnatsne", group.by = idents[ident]) + ggtitle("t-SNE") + cluster_theme
    print(p1 + p2 + p3)
    dev.off()
  }
}