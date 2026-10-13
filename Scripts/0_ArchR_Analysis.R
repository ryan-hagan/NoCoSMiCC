#### Script 1 - ArchR_Analysis.R ####

# Script for ArchR project creation, identification of epithelial cells and peak calling using ENCODE and HuBMAP colon scATAC datasets
# Adapted from Hickey et al. 2023 (doi:10.1038/s41586-023-05915-x)

#######################################
# Setup
#######################################
# load packages
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(pheatmap)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tibble)
library(FNN)
library(rtracklayer)

# set seed for reproducibility
set.seed(1)

# set genome and threads
addArchRGenome("hg38")
addArchRThreads(threads = 64) # adjust the number of cores

# set working directory (base directory should be set to where the NoCoSMiCC repo was cloned to)
base_dir <- "."
setwd(base_dir)

# set the path to Macs2 for peak calling
pathToMacs2 <- Sys.which("macs2")
# pathToMacs2 <- "/path/to/your/macs2"

# Define arrow/frag file directories
arrowdir <- file.path(base_dir, "data/ArrowFiles")
fragDir  <- file.path(base_dir, "data/scATAC")

#######################################
# Define paths to fragment files and nFrags thresholds
#######################################

# ENCODE samples
ENCODE_frags <- c(
  "ENCSR916RYB" = 1000, "ENCSR904WIW" = 1000, "ENCSR830FPR" = 1000,
  "ENCSR997YNO" = 2000, "ENCSR007QIO" = 2000, "ENCSR349XKD" = 2000,
  "ENCSR434SXE" = 2000, "ENCSR388NCA" = 2000, "ENCSR506YMX" = 2000,
  "ENCSR367GKP" = 2000
)

# HuBMAP samples
HuBMAP_frags <- c(
  # 1K threshold
  "B006-A-001" = 1000, "B006-A-101" = 1000, "B006-A-201-R2" = 1000,
  # 1.5K threshold
  "B001-A-302" = 1500, "B001-A-401" = 1500,
  "B001-A-406" = 1500, "B001-A-501" = 1500,
  # 2K threshold
  "B004-A-004" = 2000, "B008-A-001" = 2000,
  # 3K threshold
  "B001-A-301" = 3000, "B005-A-001" = 3000,
  "B005-A-002" = 3000, "B005-A-101" = 3000,
  "B005-A-201" = 3000, "B006-A-002" = 3000,
  "B006-A-201" = 3000, "B008-A-002" = 3000,
  "B008-A-101" = 3000, "B008-A-201" = 3000,
  "B010-A-001" = 3000, "B010-A-002" = 3000,
  "B010-A-101" = 3000, "B011-A-001" = 3000,
  "B011-A-002" = 3000, "B011-A-201" = 3000,
  "B012-A-001" = 3000, "B012-A-002" = 3000,
  "B012-A-101" = 3000,
  # 4K threshold
  "B004-A-004-R2" = 4000, "B004-A-204" = 4000,
  "B009-A-101" = 4000, "B010-A-201" = 4000,
  "B011-A-101" = 4000,
  # 5K threshold
  "B004-A-008" = 5000, "B009-A-001" = 5000,
  # 6K threshold
  "B012-A-201" = 6000
)

# define fragment file paths
HuBMAP_inputFiles <- setNames(
  file.path(fragDir, paste0(names(HuBMAP_frags), "_atac_fragments.tsv.gz")),
  names(HuBMAP_frags)
)
ENCODE_inputFiles <- setNames(
  file.path(fragDir, paste0(names(ENCODE_frags), ".fragments.tsv.gz")),
  names(ENCODE_frags)
)

all_frags <- c(HuBMAP_inputFiles, ENCODE_inputFiles)

#######################################
# Create arrow files
#######################################

dir.create(arrowdir, recursive = TRUE, showWarnings = FALSE)

all_thresholds <- unique(c(HuBMAP_frags, ENCODE_frags))

for(threshold in sort(all_thresholds)) {
  
  hubmap_samples <- names(HuBMAP_frags)[HuBMAP_frags == threshold]
  encode_samples <- names(ENCODE_frags)[ENCODE_frags == threshold]
  all_samples <- c(hubmap_samples, encode_samples)
  
  cat("\nCreating Arrow files with minFrags =", threshold, "\n")
  cat("Samples:", paste(all_samples, collapse = ", "), "\n")
  
  inputFiles <- all_frags[all_samples]
  
  createArrowFiles(
    inputFiles  = inputFiles,
    sampleNames = names(inputFiles),
    outputNames = file.path(arrowdir, names(inputFiles)),
    minFrags    = threshold,
    minTSS      = 5,
    addTileMat  = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE
  )
}

#######################################
# Calculate doublet scores for ENCODE samples
#######################################

ENCODE_Arrows <- file.path(arrowdir, paste0(names(ENCODE_frags), ".arrow"))
names(ENCODE_Arrows) <- names(ENCODE_frags)

doubletScores <- addDoubletScores(ENCODE_Arrows, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

#######################################
# Create an ArchR project for ENCODE samples
#######################################

ENCODEproj <- ArchRProject(
  ArrowFiles = ENCODE_Arrows,
  outputDirectory = "data/ArchR_projects/ArchR_ENCODE")

#######################################
# Filter doublets
#######################################

ENCODEproj <- filterDoublets(ENCODEproj, filterRatio = 1.2)

saveArchRProject(ArchRProj = ENCODEproj, 
                 outputDirectory = file.path(base_dir, "outputs/ArchR_projects/ArchR_ENCODE"))
#ENCODEproj <- loadArchRProject(file.path(base_dir,"outputs/ArchR_projects/ArchR_ENCODE"))

#######################################
# Plot quality control metrics
#######################################

p1 <- plotGroups(
  ArchRProj = ENCODEproj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.8,
  addBoxPlot = TRUE
)

p2 <- plotGroups(
  ArchRProj = ENCODEproj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

plotPDF(p1, p2,
  name = "ENCODEproj_QC",
  ArchRProj = ENCODEproj,
  addDOC = FALSE,
  width = 9, height = 6
)

#######################################
# Dim reduction, clustering and epithelial cluster identification
#######################################

# define marker genes
markerSets <- list(
  General_Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  T_cells = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B","TBX21", "IL7R", "CD4", "CD2"),
  B_cells = c("PAX5", "MS4A1", "CD19", "IGLL5", "VPREB3"),
  Myeloid = c("CD14", "CD68", "CSF1R", "LYZ", "ITGAM"),
  NK_cells = c("NCAM1", "KLRD1", "NKG7", "GNLY"),
  Fibroblasts = c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "FAP", "CBLN2", "SPOCK1", "ACSS3"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  SmoothMuscle = c("ACTA2", "MYH11", "TAGLN")
)

allSamples <- unique(ENCODEproj$Sample)

for(sample in allSamples) {
  
  cat("Processing sample:", sample, "\n")
  
  outDir <- file.path(base_dir, "outputs", paste0("ArchR_", sample))
  dir.create(outDir, showWarnings = FALSE)
  
  sampCells <- getCellNames(ENCODEproj)[ENCODEproj$Sample == sample]
  cat("Total cells:", length(sampCells), "\n")
  
  sampProj <- subsetArchRProject(
    ArchRProj = ENCODEproj,
    cells = sampCells,
    outputDirectory = outDir,
    force = TRUE
  )
  
  sampProj <- addIterativeLSI(
    ArchRProj = sampProj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(resolution = c(0.2), 
                         sampleCells = 10000, 
                         n.start = 10),
    varFeatures = 25000,
    dimsToUse = 1:30,
    seed = 1,
    force = TRUE
  )
  
  sampProj <- addClusters(
    input = sampProj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    seed = 1,
    force = TRUE
  )
  
  print(table(sampProj$Clusters))
  
  sampProj <- addUMAP(
    ArchRProj = sampProj,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    seed = 1,
    force = TRUE
  )
  
  # Overview UMAPs
  sample_umap <- plotEmbedding(
    ArchRProj = sampProj,
    colorBy = "cellColData",
    name = "Sample",
    embedding = "UMAP"
  )
  
  cluster_umap <- plotEmbedding(
    ArchRProj = sampProj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP"
  )
  
  # read GeneScoreMatrix
  gsMatrix <- getMatrixFromProject(sampProj, useMatrix = "GeneScoreMatrix")
  allMarkers <- unlist(markerSets)
  availMarkers <- allMarkers[allMarkers %in% rowData(gsMatrix)$name]
  
  gsRows <- which(rowData(gsMatrix)$name %in% availMarkers)
  gsMat  <- assay(gsMatrix)[gsRows, ]
  rownames(gsMat) <- rowData(gsMatrix)$name[gsRows]
  rm(gsMatrix); gc()
  
  clusters <- sampProj$Clusters
  
  meanScores <- sapply(sort(unique(clusters)), function(cl) {
    cells <- which(clusters == cl)
    if(length(cells) == 1) gsMat[, cells]
    else rowMeans(gsMat[, cells])
  })
  colnames(meanScores) <- sort(unique(clusters))
  
  # order markers
  orderedGenes <- unlist(markerSets)
  orderedGenes <- orderedGenes[orderedGenes %in% rownames(meanScores)]
  
  heatMatScaled <- t(scale(t(meanScores[orderedGenes, ])))
  
  geneSetAnnot <- data.frame(
    CellType = rep(names(markerSets), lengths(markerSets)),
    row.names = unlist(markerSets)
  )[orderedGenes, , drop = FALSE]
  
  p_heatmap <- pheatmap::pheatmap(
    heatMatScaled,
    annotation_row = geneSetAnnot,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    fontsize_row = 7,
    fontsize_col = 9,
    main = paste0(sample, " | Gene Activity Scores per Cluster"),
    silent = TRUE
  )
  
  # cluster scores table
  setScores <- sapply(names(markerSets), function(setName) {
    genes <- markerSets[[setName]][markerSets[[setName]] %in% rownames(gsMat)]
    if(length(genes) == 0) return(rep(NA, ncol(meanScores)))
    if(length(genes) == 1) return(meanScores[genes, ])
    colMeans(meanScores[genes, ])
  })
  
  setScoresDF <- round(as.data.frame(setScores), 3)
  setScoresDF$Cluster <- rownames(setScoresDF)
  
  cat("\nCluster scores for", sample, ":\n")
  print(setScoresDF)
  
  # save PDF
  pdf(file.path(outDir, paste0("MarkerAnalysis_", sample, ".pdf")),
      width = 15, height = 10)
  
  print(cowplot::plot_grid(
    sample_umap, cluster_umap,
    ncol = 2
  ))
  
  grid::grid.newpage()
  grid::grid.draw(p_heatmap$gtable)
  
  dev.off()
  cat("Saved PDF:", file.path(outDir, paste0("MarkerAnalysis_", sample, ".pdf")), "\n")
  
  saveArchRProject(sampProj)
  
  # tidy
  rm(sampProj, gsMat, meanScores, setScores, setScoresDF,
     p_heatmap, heatMatScaled, geneSetAnnot, orderedGenes)
  gc()
}

#######################################
# Define epithelial clusters
#######################################

ENCODE_epi_Clusters <- list(
  "ENCSR997YNO" = c("C1", "C2"),
  "ENCSR830FPR" = c("C1", "C2", "C3", "C4"),
  "ENCSR349XKD" = c("C1", "C2", "C3"),
  "ENCSR434SXE" = c("C1"),
  "ENCSR506YMX" = c("C1"),
  "ENCSR904WIW" = c("C1","C2", "C3", "C4","C5", "C6")
  # Samples with no epithelial clusters are omitted
)

allEpithelialCells <- c()

for(samp in names(ENCODE_epi_Clusters)) {
  
  cat("\nSubsetting epithelial cells from:", samp, "\n")
  
  # load the already-processed sample project
  sampProj <- loadArchRProject(file.path(base_dir, "outputs", paste0("ArchR_", samp)))
  
  # get cells from epithelial clusters
  epiClusters <- ENCODE_epi_Clusters[[samp]]
  epiCells <- getCellNames(sampProj)[sampProj$Clusters %in% epiClusters]
  
  cat("Epithelial cells:", length(epiCells), "\n")
  allEpithelialCells <- c(allEpithelialCells, epiCells)
}

cat("\nTotal epithelial cells across ENCODE samples:", length(allEpithelialCells), "\n")

# subset ENCODE project to epithelial cells
ENCODEProjEpi <- subsetArchRProject(
  ArchRProj = ENCODEproj,
  cells = allEpithelialCells,
  outputDirectory = "ArchR_Epithelial_ENCODE",
  force = TRUE
)

table(ENCODEProjEpi$Sample)

#######################################
# LSI, clustering on combined ENCODE epithelial cells
#######################################

ENCODEProjEpi <- addIterativeLSI(
  ArchRProj = ENCODEProjEpi,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  dimsToUse = 1:30,
  varFeatures = 25000,
  force = TRUE
)

ENCODEProjEpi <- addHarmony(
  ArchRProj = ENCODEProjEpi,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

ENCODEProjEpi <- addClusters(
  input = ENCODEProjEpi,
  reducedDims = "Harmony",
  resolution = 1.0,
  force = TRUE
)

ENCODEProjEpi <- addUMAP(
  ArchRProj = ENCODEProjEpi,
  reducedDims = "Harmony",
  name = "UMAP_Harmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

# plots
p1 <- plotEmbedding(
  ArchRProj = ENCODEProjEpi,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Harmony"
)

p2 <- plotEmbedding(
  ArchRProj = ENCODEProjEpi,
  colorBy = "GeneScoreMatrix",
  name = c("EPCAM", "KRT8", "KRT18", "MUC2"),
  embedding = "UMAP_Harmony"
)

p3 <- plotEmbedding(
  ArchRProj = ENCODEProjEpi,
  colorBy = "GeneScoreMatrix",
  name = c("CD3E", "CD14", "COL1A1", "FAP"),
  embedding = "UMAP_Harmony"
)

plotPDF(p1, p2, p3,
  name = "Epithelial_ENCODE_QC",
  ArchRProj = ENCODEProjEpi,
  addDOC = FALSE,
  width = 5, height = 5
)

#######################################
# Create combined epithelial project
#######################################

# define HuBMAP arrows
HuBMAP_Arrows <- file.path(arrowdir, paste0(names(HuBMAP_frags), ".arrow"))
names(HuBMAP_Arrows) <- names(HuBMAP_frags)

# create combined project
CombinedProj <- ArchRProject(
  ArrowFiles = c(ENCODE_Arrows, HuBMAP_Arrows),
  outputDirectory = "ArchR_Epithelial_Combined"
)

# read in the cell barcodes from multiome and non-multiome HuBMAP samples
HuBMAP_multiome_colon_epi_barcodes <- read.table(
  file.path(base_dir, "files/colon_epithelial_peak_matrix_cells.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)$V1

HuBMAP_Nonmultiome_colon_epi_barcodes <- read.table(
  file.path(base_dir, "files/non_multiome_colon_epithelial_peak_matrix_cells.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)$V1

# combine HuBMAP multiome and non-multiome epithelial barcodes
HuBMAPEpiCells <- union(HuBMAP_multiome_colon_epi_barcodes, HuBMAP_Nonmultiome_colon_epi_barcodes)

# Define all epithelial cells
AllEpiCells <- c(getCellNames(ENCODEProjEpi), HuBMAPEpiCells)

# Create combined epithelial project
CombinedEpiProj <- subsetArchRProject(
  ArchRProj = CombinedProj,
  cells = AllEpiCells,
  outputDirectory = "ArchR_Epithelial_Combined",
  force = TRUE
)

saveArchRProject(ArchRProj = CombinedEpiProj,
                 outputDirectory = file.path(base_dir, "outputs/ArchR_projects/ArchR_Epithelial_Combined"))
#CombinedEpiProj <- loadArchRProject(file.path(base_dir, "outputs/ArchR_projects/ArchR_Epithelial_Combined"))

#######################################
# Remove samples with < 100 epithelial cells
#######################################

table(CombinedEpiProj$Sample)

low_count_samples <- c("B004-A-004", "B005-A-101", "B005-A-201")

# get the cell names to keep
keep_cells <- getCellNames(CombinedEpiProj)[
  !(CombinedEpiProj$Sample %in% low_count_samples)
]

cat("Cells before filtering:", nCells(CombinedEpiProj), "\n")
cat("Cells to remove:", nCells(CombinedEpiProj) - length(keep_cells), "\n")
cat("Cells after filtering:", length(keep_cells), "\n")

# subset project
CombinedEpiProj <- subsetArchRProject(
  ArchRProj = CombinedEpiProj,
  cells = keep_cells,
  outputDirectory = "ArchR_Epithelial_Combined",
  force = TRUE
)

table(CombinedEpiProj$Sample)

saveArchRProject(ArchRProj = CombinedEpiProj,
                 outputDirectory = file.path(base_dir, "outputs/ArchR_projects/ArchR_Epithelial_Combined"))
#CombinedEpiProj <- loadArchRProject(file.path(base_dir, "outputs/ArchR_projects/ArchR_Epithelial_Combined"))

#######################################
# Add metadata
#######################################

donors <- c(
  # HuBMAP
  "B001-A-301" = "B001", "B001-A-302" = "B001", "B001-A-401" = "B001", "B001-A-406" = "B001", "B001-A-501" = "B001",
  "B004-A-004" = "B004", "B004-A-004-R2" = "B004", "B004-A-008" = "B004", "B004-A-204" = "B004",
  "B005-A-001" = "B005", "B005-A-002" = "B005", 
  "B006-A-001" = "B006", "B006-A-002" = "B006", "B006-A-101" = "B006", "B006-A-201" = "B006", "B006-A-201-R2" = "B006",
  "B008-A-001" = "B008", "B008-A-002" = "B008", "B008-A-101" = "B008", "B008-A-201" = "B008",
  "B009-A-001" = "B009", "B009-A-101" = "B009",
  "B010-A-001" = "B010", "B010-A-002" = "B010", "B010-A-101" = "B010", "B010-A-201" = "B010",
  "B011-A-001" = "B011", "B011-A-002" = "B011", "B011-A-101" = "B011", "B011-A-201" = "B011",
  "B012-A-001" = "B012", "B012-A-002" = "B012", "B012-A-101" = "B012", "B012-A-201" = "B012",
  # ENCODE
  "ENCSR830FPR" = "ENCSR830FPR", "ENCSR434SXE" = "ENCSR434SXE",
  "ENCSR997YNO" = "ENCSR997YNO", "ENCSR904WIW" = "ENCSR904WIW",
  "ENCSR506YMX" = "ENCSR506YMX", "ENCSR349XKD" = "ENCSR349XKD"
)

# add donor info to metadata
CombinedEpiProj <- addCellColData(
  ArchRProj = CombinedEpiProj,
  data      = donors[CombinedEpiProj$Sample],
  name      = "Donor",
  cells     = getCellNames(CombinedEpiProj),
  force     = TRUE
)

table(CombinedEpiProj$Donor)

#######################################
# Run LSI and Harmony on epithelial cells
#######################################

CombinedEpiProj <- addIterativeLSI(
  ArchRProj = CombinedEpiProj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  clusterParams = list(
        resolution = c(0.1, 0.2), 
        sampleCells = 20000, 
        n.start = 10), 
  sampleCellsPre = 25000,
  dimsToUse = 1:25,
  varFeatures = 15000,
  force = TRUE
)

CombinedEpiProj <- addHarmony(
  ArchRProj = CombinedEpiProj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Donor",
  force = TRUE
)

CombinedEpiProj <- addClusters(
  input = CombinedEpiProj,
  reducedDims = "Harmony",
  resolution = 1.0,
  force = TRUE
)

CombinedEpiProj <- addUMAP(
  ArchRProj = CombinedEpiProj,
  reducedDims = "Harmony",
  name = "UMAP_Harmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)
 
#######################################
# Plot combined epithelial UMAPs
#######################################

# add imputation weights for plotting
CombinedEpiProj <- addImputeWeights(ArchRProj = CombinedEpiProj, reducedDims = "Harmony")
 
p1 <- plotEmbedding(
  ArchRProj = CombinedEpiProj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Harmony"
)

p2 <- plotEmbedding(
  ArchRProj = CombinedEpiProj,
  colorBy = "GeneScoreMatrix",
  name = c("EPCAM", "KRT8", "KRT18", "MUC2"),
  embedding = "UMAP_Harmony"
)

p3 <- plotEmbedding(
  ArchRProj = CombinedEpiProj,
  colorBy = "GeneScoreMatrix",
  name = c("CD3E", "CD14", "COL1A1"),
  embedding = "UMAP_Harmony"
)

plotPDF(p1, p2, p3,
  name = "Epithelial_Combined_QC",
  ArchRProj = CombinedEpiProj,
  addDOC = FALSE,
  width = 5, height = 5
)

#######################################
# Annotate epithelial cell types using published labels
#######################################

# read in the previously published HuBMAP labels
multiome_labels <- read.table(
  file.path(base_dir, "files/scATAC_multiome_cell_types_epithelial_colon.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = ""
)

nonmultiome_labels <- read.table(
  file.path(base_dir, "files/scATAC_non_multiome_cell_types_epithelial_colon.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = ""
)

hubmap_labels <- bind_rows(
  multiome_labels    |> mutate(Source = "multiome"),
  nonmultiome_labels |> mutate(Source = "non_multiome")
) |> distinct(Cell, .keep_all = TRUE)

# add published cell labels to the combined project
cell_names <- getCellNames(CombinedEpiProj)

label_map <- tibble(Cell = cell_names) |>
  left_join(hubmap_labels, by = "Cell")

CombinedEpiProj <- addCellColData(
  ArchRProj = CombinedEpiProj,
  data  = label_map$CellType,
  cells = cell_names,
  name  = "CellType_HuBMAP",
  force = TRUE
)

# extract the harmony embedding
harmony_embed <- getEmbedding(
  ArchRProj = CombinedEpiProj,
  embedding = "UMAP_Harmony",
  returnDF  = TRUE
)

hubmap_idx <- which(!is.na(CombinedEpiProj$CellType_HuBMAP))
encode_idx  <- which(is.na(CombinedEpiProj$CellType_HuBMAP))

ref_embed   <- as.matrix(harmony_embed[hubmap_idx, ])
query_embed <- as.matrix(harmony_embed[encode_idx, ])
ref_labels  <- CombinedEpiProj$CellType_HuBMAP[hubmap_idx]

knn_result <- get.knnx(data = ref_embed, query = query_embed, k = 25)

encode_knn_labels <- apply(knn_result$nn.index, 1, function(idx) {
  neighbour_labels <- ref_labels[idx]
  names(sort(table(neighbour_labels), decreasing = TRUE))[1]
})

encode_knn_confidence <- apply(knn_result$nn.index, 1, function(idx) {
  neighbour_labels <- ref_labels[idx]
  max(table(neighbour_labels)) / length(neighbour_labels)
})

all_labels <- CombinedEpiProj$CellType_HuBMAP
all_labels[encode_idx] <- encode_knn_labels

confidence_all <- rep(NA_real_, length(cell_names))
confidence_all[encode_idx] <- encode_knn_confidence

CombinedEpiProj <- addCellColData(
  ArchRProj = CombinedEpiProj,
  data  = all_labels,
  cells = cell_names,
  name  = "CellTypeKNN",
  force = TRUE
)

CombinedEpiProj <- addCellColData(
  ArchRProj = CombinedEpiProj,
  data  = confidence_all,
  cells = cell_names,
  name  = "KNNConfidence",
  force = TRUE
)

cat("\n=== Final label distribution (all cells) ===\n")
print(sort(table(CombinedEpiProj$CellTypeKNN, useNA = "ifany"), decreasing = TRUE))

cat("\n=== ENCODE KNN label distribution ===\n")
print(sort(table(encode_knn_labels), decreasing = TRUE))

cat("\n=== ENCODE confidence summary per label ===\n")
tibble(label = encode_knn_labels, confidence = encode_knn_confidence) |>
  group_by(label) |>
  summarise(
    n = n(),
    median_conf = round(median(confidence), 3),
    pct_high    = round(mean(confidence > 0.6) * 100, 1)
  ) |>
  arrange(desc(n)) |>
  print()

# collapse cell types
epiCellType_map <- c(
  "Stem" = "Stem",
  "CyclingTA 1" = "Cycling TA",
  "CyclingTA 2" = "Cycling TA",
  "TA1" = "TA1",
  "TA2" = "TA2",
  "Immature Enterocytes"= "Immature Enterocytes",
  "Enterocytes" = "Enterocytes",
  "Best4+ Enterocytes" = "Best4+ Enterocytes",
  "Immature Goblet" = "Immature Goblet",
  "Goblet" = "Goblet",
  "Enteroendocrine" = "Enteroendocrine",
  "Enterochromaffin" = "Enteroendocrine",
  "EnteroendocrineUn" = "Enteroendocrine",
  "EnteroendocrineUn 1" = "Enteroendocrine",
  "NEUROG3high" = "Enteroendocrine",
  "L Cells" = "Enteroendocrine",
  "I Cells" = "Enteroendocrine",
  "D Cells" = "Enteroendocrine",
  "S Cells" = "Enteroendocrine",
  "Mo Cells" = "Enteroendocrine",
  "Tuft" = "Tuft"
)

epiCellType_labels <- epiCellType_map[CombinedEpiProj$CellTypeKNN]

CombinedEpiProj <- addCellColData(
  ArchRProj = CombinedEpiProj,
  data  = epiCellType_labels,
  cells = getCellNames(CombinedEpiProj),
  name  = "FinalCellType",
  force = TRUE
)

# sort
sort(table(CombinedEpiProj$FinalCellType), decreasing = TRUE)

p_ta <- plotEmbedding(
  ArchRProj = CombinedEpiProj,
  colorBy = "cellColData",
  name = "FinalCellType",
  embedding = "UMAP_Harmony"
)

plotPDF(p_ta,
        name = "CellType_Split",
        ArchRProj = CombinedEpiProj,
        addDOC = FALSE,
        width = 6, height = 6
)

# define colour palette for plotting
cell_type_colours <- c(
  "Stem" = "#E41A1C",
  "Cycling TA" = "#4A90D9",
  "TA2" = "#66C2A5",
  "TA1" = "#A6D854",
  "Goblet" = "#FC8D62",
  "Immature Goblet" = "#FFD92F",
  "Best4+ Enterocytes" = "#8DA0CB",
  "Immature Enterocytes" = "#B3B3E6",
  "Enterocytes" = "#984EA3",
  "Enteroendocrine" = "#00674F",
  "Tuft" = "#FF7F00" 
)

p_pub <- plotEmbedding(
  ArchRProj  = CombinedEpiProj,
  colorBy    = "cellColData",
  name       = "FinalCellType",
  embedding  = "UMAP_Harmony",
  pal        = cell_type_colours,
  labelMeans = FALSE,     
  size       = 0.2,        
  plotAs     = "points"
)

p_pub <- p_pub +
  theme_classic(base_size = 14) +
  theme(
    axis.line        = element_line(linewidth = 0.5),
    axis.ticks       = element_line(linewidth = 0.5),
    axis.text        = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.5, "cm"),
    plot.title       = element_blank()
  ) + guides(colour = guide_legend(
    override.aes = list(size = 3),
    ncol = 1
  )) + labs(
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    colour = "Cell Type"
  )

# save plot
ggsave(
  file.path(base_dir, "CombinedEpi_CellType_publication.png"),
  plot   = p_pub,
  width  = 7,
  height = 6,
  dpi    = 300
)

#######################################
# Call peaks
#######################################

# remove tuft cells (due to low number)
CombinedEpiProj <- CombinedEpiProj[CombinedEpiProj$FinalCellType != "Tuft", ]

CombinedEpiProj <- addGroupCoverages(
  ArchRProj = CombinedEpiProj,
  groupBy = "FinalCellType",
  force = TRUE
)

CombinedEpiProj <- addReproduciblePeakSet(
  ArchRProj = CombinedEpiProj,
  groupBy = "FinalCellType",
  pathToMacs2 = pathToMacs2,
  maxPeaks = 200000,
  force = TRUE
)

CombinedEpiProj <- addPeakMatrix(CombinedEpiProj)

# get the consensus peaks and export to BED
consensus_peaks <- getPeakSet(CombinedEpiProj)

dir.create("outputs/Peaks", recursive = TRUE, showWarnings = FALSE)
export(consensus_peaks, 
       con = file.path(base_dir, "outputs/Peaks/consensus_peaks.bed"),
       format = "BED")
