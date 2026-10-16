library(SpatialExperiment)
library(Matrix)
library(stringr)
library(scater)


#### download relevant data files from GSE176092 GEO repository
#### Day 1 and 14 perturbation data (WT, shCsrp3, Csrp3OE)


#### Create SpatialExperiment object ####

all_files <- list.files("./selected_GSE176092_RAW", 
                        full.names = FALSE)

sample_names <- unique(sub("_.*", "", all_files))
filepath <- "./selected_GSE176092_RAW"


spe_list <- lapply(sample_names, function(sid) {
    
    # get file names
    mtx_file <- 
        grep(paste0(sid, ".*matrix\\.mtx\\.gz$"), all_files, value = TRUE)
    barcode_file <- 
        grep(paste0(sid, ".*barcodes\\.tsv\\.gz$"), all_files, value = TRUE)
    feature_file <- 
        grep(paste0(sid, ".*features\\.tsv\\.gz$"), all_files, value = TRUE)
    pos_file <- 
        grep(paste0(sid, ".*tissue_positions_list\\.csv\\.gz$"), all_files, 
             value = TRUE)
    
    tissueName <- sub(paste0(sid, "_(.*?)_tissue_positions_list\\.csv\\.gz$"), 
                      "\\1", pos_file)
    expName <- str_split_i(tissueName, "_MI_", 1)
    timepoint <- str_split_i(tissueName, "_MI_", 2)
    
    # read files
    mat <- readMM(file.path(filepath, mtx_file))
    barcodes <- readLines(file.path(filepath, barcode_file))
    features <- read.delim(file.path(filepath, feature_file), header = FALSE)
    positions <- read.csv(file.path(filepath, pos_file), header = FALSE)
    colnames(positions) <- c("barcode", "in_tissue", "array_row", "array_col",
                             "pxl_row_in_img", "pxl_col_in_img")
    
    # assign barcode names
    spot_names <- paste0(sid, "_", barcodes)
    rownames(mat) <- features$V1
    colnames(mat) <- spot_names
    
    # match coordinates to barcodes
    rownames(positions) <- paste0(sid, "_", positions$barcode)
    # select only in_tissue positions
    positions <- positions[spot_names, ]
    
    # create metadata and spatialCoords
    col_data <- DataFrame(sample_id = sid, sample_group = tissueName, 
                          Condition = expName, Time_point = timepoint, positions, 
                          row.names = spot_names)
    spatial_coords <- as.matrix(positions[, c("pxl_row_in_img", "pxl_col_in_img")])
    colnames(spatial_coords) <- c('x', 'y')
    
    # create SpatialExperiment object
    spe <- SpatialExperiment(
        assays = list(counts = mat),
        colData = col_data,
        rowData = DataFrame(gene_id = features$V1, gene_name = features$V2),
        spatialCoords = spatial_coords)
    return(spe)
})

# combine all spe objects
spe <- do.call(cbind, spe_list)
spe$sample_group <- factor(
    spe$sample_group, 
    levels = c("WT_MI_day1_1", "WT_MI_day14_1", "WT_MI_day1_2", "WT_MI_day14_2",
               "WT_MI_day1_3", "WT_MI_day14_3", "shCsrp3_MI_day1", 
               "shCsrp3_MI_day14", "Csrp3OE_MI_day1", "Csrp3OE_MI_day14"))
saveRDS(spe, "./spe_Csrp3_perturb.Rds")



#### QC ####

# replace gene ensemble ids with gene names
rownames(spe) <- rowData(spe)[, 2] |> make.unique()

# QC - mitochondrial
mt_genes   <- grepl("^mt-", rownames(spe))
spe$percent_mito <- colSums(counts(spe)[mt_genes, ]) / spe$total_counts * 100
spe$percent_mito_15 <- spe$percent_mito < 15

# QC - ribosomal
ribo_genes <- grepl("^RPL|^RPS", rownames(spe))
spe$percent_ribo <- colSums(counts(spe)[ribo_genes, ]) / spe$total_counts * 100

# QC - library size
spe$total_counts <- colSums(counts(spe))
spe$total_counts_1000 <- spe$total_counts > 1000

# QC - gene detection
spe$n_genes <- colSums(counts(spe) > 0)
spe$n_genes_500 <- spe$n_genes > 500



#### Filtering ####

# filter only by library size and genes detected, similar to the approach in 
# original paper (based on their spatial plots)
spe_filtered <- spe[, spe$total_counts > 1000 &
                        spe$n_genes > 500]
saveRDS(spe_filtered, "./spe_Csrp3_perturb_QC.Rds")
