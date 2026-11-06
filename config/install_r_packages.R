options(repos = c(CRAN = Sys.getenv("CHAIR_CRAN_MIRROR", "https://cloud.r-project.org")))

# Bioconductor 3.20 UCSC sources use C semantics predating GCC 15's C23 default.
user.makevars <- Sys.getenv("R_MAKEVARS_USER", unset = path.expand("~/.R/Makevars"))
installer.makevars <- tempfile(fileext = ".Makevars")
writeLines(c(if (file.exists(user.makevars)) readLines(user.makevars) else character(),
             "CFLAGS += -std=gnu17"), installer.makevars)
Sys.setenv(R_MAKEVARS_USER = installer.makevars)
genome <- commandArgs(TRUE)
if (!length(genome)) genome <- "mm10"
if (!all(genome %in% c("mm10", "hg38"))) stop("Genome: mm10 and/or hg38")
install.packages(c("remotes", "BiocManager"))
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
install.version <- function(package, version) {
  if (!requireNamespace(package, quietly = TRUE) || packageVersion(package) != package_version(version)) {
    remotes::install_version(package, version, upgrade = "never")
  }
}
install.version("Matrix", "1.7-3")
install.version("SeuratObject", "5.1.0")
install.version("Seurat", "5.3.0")
bioc.core <- c("GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools", "S4Vectors", "BiocGenerics", "rtracklayer")
BiocManager::install(bioc.core, ask = FALSE, update = FALSE)
missing.core <- bioc.core[!vapply(bioc.core, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing.core)) stop("Bioconductor dependencies failed to install: ", paste(missing.core, collapse = ", "))
install.version("Signac", "1.14.0")
if (!requireNamespace("DoubletFinder", quietly = TRUE) || packageVersion("DoubletFinder") != package_version("2.0.4")) {
  remotes::install_github("chris-mcginnis-ucsf/DoubletFinder@03e9f37f891ef76a23cc55ea69f940c536ae8f9f", upgrade = "never")
}
# Use the Deriv release compatible with R 4.4 for ggpubr's dependencies.
install.version("Deriv", "4.2.0")
install.packages(c("harmony", "dplyr", "magrittr", "ggplot2", "patchwork", "stringr", "reshape2", "matrixStats", "data.table", "RColorBrewer", "pheatmap", "plotrix", "ggpubr", "hdf5r"))
annotations <- c("GenomicRanges", "GenomeInfoDb")
if ("mm10" %in% genome) annotations <- c(annotations, "EnsDb.Mmusculus.v79", "BSgenome.Mmusculus.UCSC.mm10")
if ("hg38" %in% genome) annotations <- c(annotations, "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install(annotations, ask = FALSE, update = FALSE)
required <- c("remotes", "BiocManager", "Matrix", "SeuratObject", "Seurat", "Signac", "DoubletFinder",
              "harmony", "dplyr", "magrittr", "ggplot2", "patchwork", "stringr", "reshape2", "matrixStats",
              "data.table", "RColorBrewer", "pheatmap", "plotrix", "ggpubr", "hdf5r", "Deriv", bioc.core, annotations)
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Packages failed to install or load: ", paste(missing, collapse = ", "))
expected <- c(Matrix = "1.7-3", SeuratObject = "5.1.0", Seurat = "5.3.0", Signac = "1.14.0", DoubletFinder = "2.0.4", Deriv = "4.2.0")
actual <- vapply(names(expected), function(package) as.character(packageVersion(package)), character(1))
versions.match <- vapply(names(expected), function(package) packageVersion(package) == package_version(expected[[package]]), logical(1))
if (!all(versions.match)) stop("Installed package versions differ from the installation recipe")
print(sessionInfo())
