library(SpatialExperiment)
library(dplyr)
library(StabMap)
library(Matrix)
library(harmony)
library(stringr)

# load all the cells with zone & fibrotic program information 
spe <- readRDS("./cartana_spe_annotated.Rds")
# subset to CF cells
spe_fib <- spe[, spe$Major.groups == "Fibroblast"]

# load reference scRNA-seq dataset from Patrick et al. 
sce_ref <- readRDS("./SCE_Patrick_Forte.rds")

#### Integration ####

# exclude CM-specific genes from both datasets
ref_genes <- setdiff(rownames(sce_ref), c("Myh6", "Fabp3"))
query_genes <- setdiff(rownames(spe_fib), c("Myh6", "Fabp3"))

# extract expression matrices
ref_mat <- Matrix(assay(sce_ref[ref_genes,], "logcounts_forte"), sparse = TRUE)
query_mat <- Matrix(assay(spe_fib[query_genes,], "logcounts"), sparse = TRUE)

# scale expressions per cell
custom_scale <- function(mat) {
    colnames_mat <- colnames(mat)
    
    # Apply custom scaling column-wise
    scaled_mat <- t(apply(mat, 1, function(x) {
        if (var(x) == 0) {
            return(rep(0, length(x)))  # Return zeros if variance is 0
        } else {
            return(scale(x, center = TRUE, scale = TRUE))  # Scale if variance is non-zero
        }
    }))
    
    colnames(scaled_mat) <- colnames_mat
    return(scaled_mat)
}
ref_mat <- Matrix(t(custom_scale(t(ref_mat))), sparse = TRUE)
query_mat <- Matrix(t(custom_scale(t(query_mat))), sparse = TRUE)


# define dataset list
mat_list <- list(ref_mat = ref_mat, query_mat = query_mat)

# integrate query and reference data with StabMap
stab_embed <- stabMap(assay_list = mat_list, reference_list = names(mat_list),
                      plot = FALSE)

# generate UMAP of StabMAP embedding
stab_umap <- calculateUMAP(t(stab_embed)) 

# use harmony for batch correction of datasets
labels <- ifelse(rownames(stab_embed) %in% colnames(mat_list[["ref_mat"]]), 
                 names(mat_list)[1], names(mat_list)[2])
names(labels) <- rownames(stab_embed)
harmony_embed <- RunHarmony(stab_embed, labels, early_stop = early_stop, max_iter = max_iter)

# generate UMAP of harmony embedding
harmony_umap = calculateUMAP(t(harmony_embed))

result <- list(setNames(list(stab_embed, stab_umap), 
                        c("StabMap_PCA", "StabMap_UMAP")), 
               setNames(list(harmony_embed, harmony_umap), 
                        c("Harmony_PCA", "Harmony_UMAP")))
saveRDS(result, "./stabmap_integrate.Rds")




#### Classification ####

# use harmony-corrected stabmap embeddings for label transfer
embed_mat <- harmony_embed

# data frame of matching timepoint names in both datasets
condition_df <- data.frame(Cartana = levels(spe_fib$condition),
                           Patrick_forte = c("Healthy", "MI-day 1", "MI-day 3",
                                             "MI-day 5", "MI-day 7", "MI-day 14"),
                           row.names = levels(spe_fib$condition))

# classify CF states per timepoint
celllabels_list <- setNames(lapply(condition_df$Cartana, function(tp) {
    # query cells from timepoint for which labels to be predicted
    queryCells <- colnames(spe_fib[, spe_fib$condition == tp])
    
    # label transfer from matching timepoint
    ref_tp_exact <- condition_df[tp, "Patrick_forte"]
    sce_exact <- sce_ref[, sce_ref$Condition == ref_tp_exact]
    ctLabels_exact <- setNames(
        as.character(sce_exact$Clusters_updated_knn_basic_k25), 
        colnames(sce_exact))
    # classify cells
    ctLabels_exact_stab <- StabMap::classifyEmbedding(coords = embed_mat, 
                                                      labels = ctLabels_exact)
    
    out_df <- data.frame(
        Exact_transfer = ctLabels_exact_stab[queryCells, "predicted_labels"],
        row.names = queryCells)
    return(out_df)
}), condition_df$Cartana)




#### Update with predicted labels ####

celllabels_list_reord <- celllabels_list[as.character(unique(spe_fib$condition))]

cell_labels <- do.call(rbind, celllabels_list_reord)
cell_labels$condition <- str_split_i(rownames(cell_labels), "\\.", i = 1)
cell_labels$cell_IDs <- str_split_i(rownames(cell_labels), "\\.", i = 2)
rownames(cell_labels) <- cell_labels$cell_IDs

spe_fib$stabmapHarmonyDataset_lt_exact <- 
    factor(cell_labels[colnames(spe_fib), "Exact_transfer"], 
           levels = levels(spe_fib$patrick_label_transfer))
names(spe_fib$stabmapHarmonyDataset_lt_exact) <- colnames(spe_fib)

# adding predicted labels to cartana dataset with all cells
spe$stabmapHarmonyDataset_lt_exact <- 
    spe_fib$stabmapHarmonyDataset_lt_exact[colnames(spe)]
names(spe$stabmapHarmonyDataset_lt_exact) <- NULL
spe$stabmapHarmonyDataset_lt_exact <- 
    factor(spe$stabmapHarmonyDataset_lt_exact,
           levels = levels(spe$patrick_label_transfer))

saveRDS(spe, file = "./processed_data/cartana_spe_annotated_relabel.Rds")

