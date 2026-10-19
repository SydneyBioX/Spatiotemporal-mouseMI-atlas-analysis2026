library(SpatialExperiment)
library(scater)
library(scales)
library(dplyr)
library(UCell)
library(ggplot2)
library(tibble)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(BiocNeighbors)
library(stringr)
library(reshape2)

set.seed(2026)


# # load Yamada et al processed data with QC - see Data_Analysis_Scripts
spe_yam <- readRDS("./spe_Csrp3_perturb_QC.Rds")
# rename rownames by gene names in place of geneIDs
rownames(spe_yam) <- rowData(spe_yam)[, 2] |> make.unique()
# calculate log normalised counts
spe_yam <- logNormCounts(spe_yam)

# add timepoint column
spe_yam$Time <- ifelse(
  spe_yam$Time_point %in% c("day1_1", "day1_2", "day1_3", "day1"), "Day1",
  ifelse(spe_yam$Time_point %in% 
           c("day14_1", "day14_2", "day14_3", "day14"), "Day14", NA))

perturb_col <- setNames(c("#3ffbfd", "black", "#f828f8"),
                        c("shCsrp3", "WT", "Csrp3OE"))

#### Fig. 7A ####

# subset to Day1 samples
spe_day1 <- spe_yam[, spe_yam$Time == "Day1"]

markers_ct <- c("CMs", "strCMs", "CFs")

sig_df <- data.frame(x = spatialCoords(spe_day1)[,1],
                     y = spatialCoords(spe_day1)[,2],
                     samples = spe_day1$sample_group) %>%
    # sample factors
    mutate(samples = factor(samples, 
                            levels = c("shCsrp3_MI_day1", "WT_MI_day1_1", 
                                       "WT_MI_day1_2", "WT_MI_day1_3",
                                       "Csrp3OE_MI_day1", "shCsrp3_MI_day14", 
                                       "WT_MI_day14_1", "WT_MI_day14_2", 
                                       "WT_MI_day14_3", "Csrp3OE_MI_day14")))
# define sample labels
cond_labels <- setNames(
    c("Day1\nCsrp3 KO", "Day1\nWildtype_1", "Day1\nWildtype_2", "Day1\nWildtype_3", 
      "Day1\nCsrp3 OE", "Day14\nCsrp3 KO", "Day14\nWildtype_1", "Day14\nWildtype_2", 
      "Day14\nWildtype_3", "Day14\nCsrp3 OE"), 
    levels(sig_df$samples))

sig_plt_list <- setNames(lapply(markers_ct, function(gene) {
    if (gene == "CMs") {
        CM_genes = c("Myh6", "Gja1", "Tnnc1", "Cdh2")
        sig_df <- sig_df %>%
            mutate(gene_exp = colMeans(logcounts(spe_day1)[CM_genes, ]),
                   gene_exp_z = scale(gene_exp)) %>%
            arrange(gene_exp_z) 
    } else if (gene == "strCMs") {
        strCM_genes = c("Nppa", "Clu", "Cxcr6", "Csrp3", "Des", "Ankrd1")
        sig_df <- sig_df %>%
            mutate(gene_exp = colMeans(logcounts(spe_day1)[strCM_genes, ]),
                   gene_exp_z = scale(gene_exp)) %>%
            arrange(gene_exp_z) 
    } else if (gene == "CFs") {
        Fib_genes = c("Col1a1", "Postn", "Eln")
        sig_df <- sig_df %>%
            mutate(gene_exp = colMeans(logcounts(spe_day1)[Fib_genes, ]),
                   gene_exp_z = scale(gene_exp)) %>%
            arrange(gene_exp_z) 
    }
    
    plt <- ggplot(sig_df, aes(x = x, y = y, color = gene_exp_z)) +
        geom_point(size = 1) +
        labs(color = gene) +
        facet_grid(samples ~ "", switch = "y") +
        scale_color_gradientn(colors = c("#5a53a4", "#5ab0ab", "#fffbbb",
                                         "#f6794b", "#a10741"),
                              values = scales::rescale(c(-2, -1, 0, 1, 2)),
                              limits = c(-2, 2), 
                              oob = scales::squish) +
        theme_classic() +
        theme(plot.background = element_rect(fill = "black"),
              panel.background = element_rect(fill = "black"),
              panel.border = element_blank(),
              strip.background.x = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              legend.title = element_text(size = 12, colour = "white"),
              legend.background = element_rect(fill = "black"),
              legend.text = element_text(size = 10, colour = "white"))
    return(plt)
}), markers_ct)




#### Fig. 7B ####

# subset to Day1 samples
spe_day1 <- spe_yam[, spe_yam$Time == "Day1"]

markers_BZ <- c("Csrp3", "Ankrd1", "Des", "Nppa")

gene_df <- data.frame(x = spatialCoords(spe_day1)[,1],
                      y = spatialCoords(spe_day1)[,2],
                      samples = spe_day1$sample_group) %>%
    # sample factors
    mutate(samples = factor(samples, 
                            levels = c("shCsrp3_MI_day1", "WT_MI_day1_1", 
                                       "WT_MI_day1_2", "WT_MI_day1_3",
                                       "Csrp3OE_MI_day1", "shCsrp3_MI_day14", 
                                       "WT_MI_day14_1", "WT_MI_day14_2", 
                                       "WT_MI_day14_3", "Csrp3OE_MI_day14")))
# define sample labels
cond_labels <- setNames(
    c("Day1\nCsrp3 KO", "Day1\nWildtype_1", "Day1\nWildtype_2", "Day1\nWildtype_3", 
      "Day1\nCsrp3 OE", "Day14\nCsrp3 KO", "Day14\nWildtype_1", "Day14\nWildtype_2", 
      "Day14\nWildtype_3", "Day14\nCsrp3 OE"), 
    levels(gene_df$samples))

sig_plt_list <- setNames(lapply(markers_BZ, function(gene) {
    gene_df <- gene_df %>%
        mutate(gene_exp = as.numeric(logcounts(spe_day1)[gene,]),
               gene_exp_z = scale(gene_exp)) %>%
        arrange(gene_exp_z) 
    
    plt <- ggplot(gene_df, aes(x = x, y = y, color = gene_exp_z)) +
        geom_point(size = 1) +
        labs(color = gene) +
        facet_grid(samples ~ "", switch = "y") +
        scale_color_gradientn(colors = c("#5a53a4", "#5ab0ab", "#fffbbb",
                                         "#f6794b", "#a10741"),
                              values = scales::rescale(c(-2, -1, 0, 1, 2)),
                              limits = c(-2, 2), 
                              oob = scales::squish) +
        theme_classic() +
        theme(plot.background = element_rect(fill = "black"),
              panel.background = element_rect(fill = "black"),
              panel.border = element_blank(),
              strip.background.x = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              legend.title = element_text(size = 12, colour = "white"),
              legend.background = element_rect(fill = "black"),
              legend.text = element_text(size = 10, colour = "white"))
    return(plt)
}), markers_BZ)




#### Fig. 7C ####

# subset to Day1 samples
spe_day1 <- spe_yam[, spe_yam$Time == "Day1"]

markers_imm <- c("IR signature", "Immune", "Ccr2", "Neutrophils")

sig_df2 <- data.frame(x = spatialCoords(spe_day1)[,1],
                     y = spatialCoords(spe_day1)[,2],
                     samples = spe_day1$sample_group) %>%
    # sample factors
    mutate(samples = factor(samples, 
                            levels = c("shCsrp3_MI_day1", "WT_MI_day1_1", 
                                       "WT_MI_day1_2", "WT_MI_day1_3",
                                       "Csrp3OE_MI_day1", "shCsrp3_MI_day14", 
                                       "WT_MI_day14_1", "WT_MI_day14_2", 
                                       "WT_MI_day14_3", "Csrp3OE_MI_day14")))
# define sample labels
cond_labels <- setNames(
    c("Day1\nCsrp3 KO", "Day1\nWildtype_1", "Day1\nWildtype_2", "Day1\nWildtype_3", 
      "Day1\nCsrp3 OE", "Day14\nCsrp3 KO", "Day14\nWildtype_1", "Day14\nWildtype_2", 
      "Day14\nWildtype_3", "Day14\nCsrp3 OE"), 
    levels(sig_df2$samples))

sig_plt_list2 <- setNames(lapply(markers_imm, function(gene) {
    if (gene == "IR signature") {
        IR_genes = c("Mt2", "Ccl2", "Timp1", "Cxcl5", "Angptl4")
        sig_df2 <- sig_df2 %>%
            mutate(gene_exp = colMeans(logcounts(spe_day1)[IR_genes, ]),
                   gene_exp_z = scale(gene_exp)) %>%
            arrange(gene_exp_z) 
    } else if (gene == "Immune") {
        imm_genes = c("Ptprc", "Cd48", "Spn", "Itgam")
        sig_df2 <- sig_df2 %>%
            mutate(gene_exp = colMeans(logcounts(spe_day1)[imm_genes, ]),
                   gene_exp_z = scale(gene_exp)) %>%
            arrange(gene_exp_z) 
    } else if (gene == "Neutrophils") {
        neut_genes = c("S100a8", "S100a9", "Ly6g", "Csf3r", "Cxcr2", "Retnlg", 
                       "Mmp9", "Ngp", "Camp", "Ltf")
        sig_df2 <- sig_df2 %>%
            mutate(gene_exp = colMeans(logcounts(spe_day1)[neut_genes, ]),
                   gene_exp_z = scale(gene_exp)) %>%
            arrange(gene_exp_z) 
    } else {
        sig_df2 <- sig_df2 %>%
            mutate(gene_exp = as.numeric(logcounts(spe_day1)[gene,]),
                   gene_exp_z = scale(gene_exp)) %>%
            arrange(gene_exp_z) 
    }
    
    plt <- ggplot(sig_df2, aes(x = x, y = y, color = gene_exp_z)) +
        geom_point(size = 1) +
        labs(color = gene) +
        facet_grid(samples ~ "", switch = "y") +
        scale_color_gradientn(colors = c("#5a53a4", "#5ab0ab", "#fffbbb",
                                         "#f6794b", "#a10741"),
                              values = scales::rescale(c(-2, -1, 0, 1, 2)),
                              limits = c(-2, 2), 
                              oob = scales::squish) +
        theme_classic() +
        theme(plot.background = element_rect(fill = "black"),
              panel.background = element_rect(fill = "black"),
              panel.border = element_blank(),
              strip.background.x = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              legend.title = element_text(size = 12, colour = "white"),
              legend.background = element_rect(fill = "black"),
              legend.text = element_text(size = 10, colour = "white"))
    return(plt)
}), markers_imm)




#### Fig. 7D ####

# calculate UCell score for each spot to assess its signature enrichment
ct_gene_sigs <- list(
  IR = c("Mt2", "Ccl2", "Timp1", "Cxcl5", "Angptl4"),
  CF = c("Col1a1", "Postn", "Eln"),
  Immune = c("Ptprc", "Cd48", "Spn", "Itgam"),
  StressedCM = c("Nppa", "Clu", "Cxcr6", "Csrp3", "Des", "Ankrd1"),
  CM = c("Myh6", "Gja1", "Tnnc1", "Cdh2"))

spe_yam <- ScoreSignatures_UCell(spe_yam, features = gene_sigs, maxRank = 5000, 
                             ncores = 8)
colData(spe_yam) <- cbind(colData(spe_yam), t(assay(altExp(spe_yam), "UCell"))) 
spe_yam$Condition <-factor(spe_yam$Condition, levels = c("WT", "shCsrp3", "Csrp3OE"))

sig_enrich_plots <- lapply(names(gene_sigs), function(ct) {
  col <- paste0(ct, "_UCell")
  plt <- as.data.frame(colData(spe_yam)) %>%
    ggplot(aes(x = get(col), color = Condition, fill = Condition)) +
    geom_density(bw = 0.2, linewidth = 0.8, alpha = 0.1) +
    facet_wrap(~Time, nrow = 1) +
    xlim(0, 1) +
    ylim(0, 2) +
    labs(x = paste(ct, "signature enrichment"), y = "Spot density") +
    scale_fill_manual(values = perturb_col) +
    scale_color_manual(values = perturb_col) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
})




#### Fig. 7E ####

# list of cell type and subgroup marker genes
grp_gene_sigs <- list(
  # fibroblasts
  IR = c("Mt2", "Ccl2", "Timp1", "Cxcl5", "Angptl4"),
  all_CFs = c("Col1a1", "Postn", "Eln"),
  # immune cells
  Immune = c("Ptprc", "Cd48", "Spn", "Itgam"),
  Neutrophils = c("S100a8", "S100a9", "Csf3r", "Cxcr2", "Retnlg", "Mmp9", "Ngp", 
                  "Camp", "Ltf", "Ly6g"),
  Inflammatory_macs = c("Cd14", "Ccr2", "Fcgr1"),
  Reparative_macs = c("Cx3cr1", "C1qa", "C1qb", "Cd163", "Folr2", "Mgl2", "Csf1r", 
                      "Trem2"),
  # cardiomyocytes
  StrCM = c("Nppa", "Clu", "Cxcr6", "Csrp3", "Des", "Ankrd1"),
  CM = c("Myh6", "Gja1", "Tnnc1", "Cdh2"))

# calculate cell type and subgroup mean expression signatures
for (ct in names(grp_gene_sigs)) {
  spe_yam[[paste0(ct, "_sig")]] <-
    colMeans(logcounts(spe_yam)[grp_gene_sigs[[ct]], , drop = FALSE])
}

# calculate correlation between all combinations of cell groups in each sample
cor_res <- setNames(lapply(levels(spe_yam$sample_group), function(samp) {
  # subset data to the sample
  spe_sub <- spe_yam[, spe_yam$sample_group == samp]
  
  # initiate an empty data frame to collect relevant correlation data
  tmp_df <- data.frame(
    matrix(nrow = 0, ncol = 5,
           dimnames = list(c(), c("Sample", "Day", "Index_Celltype", 
                                  "Target_Celltype", "Correlation"))))
  for (index_ct in paste0(names(grp_gene_sigs), "_sig")) { # for index group
    for (target_ct in paste0(names(grp_gene_sigs), "_sig")) { # for target group
      rho <- suppressWarnings(
        cor(x = spe_sub[[index_ct]], y = spe_sub[[target_ct]],
            method = "spearman", use = "pairwise.complete.obs"))
      tmp_df <- rbind(tmp_df,
                      data.frame(Sample = samp,
                                 Day = unique(spe_sub$Time),
                                 Index_Celltype = index_ct,
                                 Target_Celltype = target_ct,
                                 Correlation = rho,
                                 stringsAsFactors = FALSE))
      print(paste(index_ct, target_ct))
    }
  }
  return(tmp_df)
}), levels(spe_yam$sample_group))
# combine correlation data
final_cor_df <- do.call(rbind, cor_res)

# plot correlation between IR and cell groups in Day 1
correl_mat <- final_cor_df[final_cor_df$Index_Celltype == "IR_sig" &
                             final_cor_df$Day == "Day1", ] %>%
  mutate(Target_Celltype =
           recode(Target_Celltype,
                  "StrCM_sig" = "Str CM", "all_CFs_sig" = "All CFs", 
                  "CM_sig" = "CM", "Immune_sig" = "Immune", 
                  "Neutrophils_sig" = "Neutrophils",
                  "Inflammatory_macs_sig" = "Inflammatory\nmacs",
                  "Reparative_macs_sig" = "Reparative\nmacs"),
         Sample = recode(
           Sample, "WT_MI_day1_1" = "WT 1",
           "WT_MI_day1_2" = "WT", "WT_MI_day1_3" = "WT 3",
           "shCsrp3_MI_day1" = "shCsrp3",
           "Csrp3OE_MI_day1" = "Csrp3 OE")) %>%
  select(Sample, Target_Celltype, Correlation) %>%
  pivot_wider(names_from = Sample, values_from = Correlation) %>%
  column_to_rownames("Target_Celltype") %>%
  as.matrix()

selected_ct <- c("CM", "All CFs", "Str CM", "Immune", "Neutrophils",
                 "Inflammatory\nmacs", "Reparative\nmacs")
selected_samp <- c("WT", "shCsrp3", "Csrp3 OE")

correl_mat <- correl_mat[selected_ct, selected_samp]

color_info <- colorRamp2(c(min(correl_mat), 0, max(correl_mat)),
                         c("#2166ac", "white", "#b2182b"))
Heatmap(correl_mat,
        name = "IR\ncorrelation",
        col = color_info,
        rect_gp = gpar(col = "white", lwd = 1),
        row_order = selected_ct,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 12),
        column_order = selected_samp,
        column_names_side = "top",
        column_names_rot = 0,
        column_names_centered = TRUE,
        column_names_gp = gpar(fontsize = 12),
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 10)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", correl_mat[i, j]), x, y,
                    gp = gpar(fontsize = 10))
        })




#### Fig. 7F & Extended Fig. 8 ####

# load processed ISS data
spe <- readRDS("./spe_processed_final.Rds")

# organise samples
spe$sample_id <- factor(spe$sample_id, 
                        levels = c("Sham_1", "Sham_2", "Day1_1", "Day1_2", 
                                   "Day3_1", "Day3_2", "Day5_1", "Day5_2", 
                                   "Day7_1", "Day7_2", "Day14_1", "Day14_2"))

# exclude F-WntX cells with false WntX signal (only keep Wif1+ F-WntX cells) 
# add Wif1 expression per cell
spe$Wif1_exp <- as.numeric(assay(spe, "logcounts_area")["Wif1",])
# add Wif1 detection status
spe$Wif1_detetct <- ifelse(spe$Wif1_exp == 0, "No", "Yes")
# relabel F-WntX to F-WntX (Wif1), if the cell has Wif1 expression, else NA
spe$fibStates <- ifelse(spe$stabmapHarmonyDataset_lt_exact == "F-WntX",
                        ifelse(spe$Wif1_detetct == "Yes", "F-WntX_Wif1",
                               NA),
                        as.character(spe$stabmapHarmonyDataset_lt_exact))
spe$fibStates <- factor(spe$fibStates, 
                        levels = c("F-SH", "F-SL", "F-Trans", "F-WntX_Wif1", 
                                   "F-Act", "IR", "F-CI", "F-Cyc", "MYO",
                                   "MFC", "F-IFNS"))

# add cell ids to fibroblast substates
names(spe$fibStates) <- colnames(spe)
# reorganise major cell group labels
spe$Major.groups <- factor(spe$Major.groups, 
                           levels = c("Stressed cardiomyocyte", "Cardiomyocyte",
                                      "Immune", "Endothelial", "Fibroblast"))
names(spe$Major.groups) <- colnames(spe)

# create cell labels that are a mix of major groups and fibroblast substates
spe$cell_annots <- ifelse(spe$Major.groups == "Fibroblast",
                          as.character(spe$fibStates), 
                          as.character(spe$Major.groups))
names(spe$cell_annots) <- colnames(spe)
celltype_labels <- unique(spe$cell_annots)

# remove F-WntX cells with no Wif1 expression
spe <- spe[, !is.na(spe$cell_annots)]


# DETECT 5-NEAREST NEIGHBOURS

# initiate matrices for collecting relevant information
nnCells <- matrix(nrow = 0, ncol = 5) # cell ids of 5-nn
for (s in unique(spe$sample_id)) {
  # subset to sample to look for neighbours in respective samples
  speX <- spe[, spe$sample_id == s]
  coords <- reducedDim(speX, "spatial") # get locations
  knn_res <- findKNN(coords, k = 5, num.threads = 8) # find 5-nn
  # extract 5-nn indices
  k_index <- knn_res$index
  # extract cell ids of 5-nn using their row indices in sample
  nn <- matrix(colnames(speX)[k_index], dim(k_index))
  rownames(nn) <- rownames(coords)
  nnCells <- rbind(nnCells, nn)
  print(s)
}


# FUNCTIONS FOR NEIGHBOURHOOD ANALYSES

findNeighCount <- function(spe, nnCells, samples, k = 5, celltype_labels, 
                           index_labels) {
  
  # function to identify the cell type or substate labels of 5-nn
  
  coloc_matList <- setNames(lapply(samples, function(s) {
    # subset sample
    speX <- spe[, spe$sample_id == s]
    # initialise neighbour count matrix 
    coloc_mat <- matrix(0, nrow = length(index_labels), 
                        ncol = length(celltype_labels), 
                        dimnames = list(index_labels, celltype_labels))
    # extract cell ids of substate fibroblast cells
    fib_cells <- colnames(speX[, speX$Major.groups == "Fibroblast"])
    
    for(i in fib_cells) {
      # identify samples
      source_grp <- speX$fibStates[i]
      # select cell neighbours
      neigh_ids <- nnCells[i, 1:k] 
      # identify neighbour celltype annotation
      neigh_grps <- speX$cell_annots[neigh_ids]
      # calculate the number of neighbouring cells belonging to each celltype
      coloc_count <- table(factor(neigh_grps, levels = celltype_labels))
      # convert total counts to binary detection calls 
      # detected = 1, not detected = 0
      coloc_binary <- setNames(as.integer(coloc_count > 0), 
                               names(coloc_count))
      # sum neighbour counts for each sample
      coloc_mat[source_grp, ] <- coloc_mat[source_grp, ] + coloc_binary
    }
    return(coloc_mat)
  }), samples)
  print(s)
  return(coloc_matList)
}


NN_meanlog_substate_heatmap <- function(spe, coloc_df_all, p_mat, 
                                        samples_ord, substate_order) {
  # function to generate heatmap of mean of log ratios of 
  # nearest neighbour proportions and the statistical inference
  
  substate <- substate_order[1] # main substate
  non_substate <- substate_order[2] # all other fibroblasts
  
  # remove sample and source label names
  coloc_df_all[, c("Sample", "Source_labels", "Condition")] <- NULL
  
  # calculate neighbor proportions
  coloc_prop <- as.matrix(coloc_df_all) / rowSums(as.matrix(coloc_df_all))
  
  # get mean log ratio of substate/non-substate cells
  
  coloc_df <- coloc_prop %>%
    as.data.frame() %>%
    # add column for substate and non-substate fibroblasts
    mutate(Source_label = str_split_i(rownames(coloc_prop), "---", i = 2),
           Sample = factor(str_split_i(rownames(coloc_prop), "---", i = 1),
                           levels = samples_ord)) %>%
    # melt dataframe
    pivot_longer(!c(Sample, Source_label), 
                 names_to = "Neighbouring_celltypes", 
                 values_to = "counts") %>%
    # elongate dataframe to compare between fibrotic classes
    pivot_wider(names_from = Source_label, values_from = counts) %>%
    # calculate and add log ratio of substate vs non-substate 
    mutate(logratio = log2((get(substate) + 1) / 
                             (get(non_substate) + 1))) %>%
    # select specific columns to elongate
    select(Sample, Neighbouring_celltypes, logratio) %>%
    # add Day and repeat columns
    mutate(Day = str_split_i(Sample, "_", i = 1),
           Repeat = str_split_i(Sample, "_", i = 2),
           Sample = NULL) %>%
    # select specific columns to elongate
    select(Day, Repeat, Neighbouring_celltypes, logratio) %>%
    # elongate data frame to compare between repeats
    pivot_wider(names_from = Repeat, values_from = logratio) %>%
    # calculate mean of logRatio in sample repeats
    mutate(MeanLogRatio = (`1` + `2`) / 2) %>%
    # select specific columns to elongate
    select(Day, Neighbouring_celltypes, MeanLogRatio) %>%
    # elongate dataframe to sample x neighbourhood with logratio values
    pivot_wider(names_from = Neighbouring_celltypes, 
                values_from = MeanLogRatio) %>%
    as.data.frame() # converting to data frame
  # add timepoint and fibroblast substate names as rownames
  rownames(coloc_df) <- coloc_df$Day 
  
  # exclude sample names and convert to numeric matrix
  coloc_mat <- as.matrix(coloc_df[, -c(1)])
  
  # count substate cells in each condition   
  substate_counts <- 
    table(spe$condition, spe$fibStates)[, substate]
  # create top annotation bar plot for the heatmap
  col_annot <- do.call(HeatmapAnnotation, setNames(
    list(anno_barplot(substate_counts, gp = gpar(fill = "black"))),
    paste0("Total\n", substate, " cells")))
  # define color range
  col_func <- colorRamp2(c(min(coloc_mat, na.rm = TRUE), 0, 
                           max(coloc_mat, na.rm = TRUE)), 
                         c("blue", "white", "red"))
  
  # Check p_mat rownames before using
  if (identical(rownames(p_mat), colnames(coloc_mat)) != TRUE) {
    p_mat <- p_mat[colnames(coloc_mat), ]
  } 
  
  # create heatmap
  if (substate == "F-WntX_Wif1") {
    substate = "F-WntX"
    non_substate = "nonF-WntX"
  }
  colocal_ht <- Heatmap(
    t(coloc_mat), 
    name = paste0("Logratio of\n", substate,"/", non_substate, 
                  "\nFibroblast cells"),
    col = col_func,
    na_col = "grey95",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    border_gp = grid::gpar(col = "black", lty = 2, lwd = 1),
    rect_gp = grid::gpar(col = "white", lwd = 1),
    column_title_gp = gpar(fontsize = 14),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12),
    column_names_side = "top",
    top_annotation = col_annot,
    # add adjusted p-values
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- p_mat[i, j]
      if (!is.na(p)) {
        if (p <= 0.001) {
          grid.text("***", x, y)
        } else if (p <= 0.01) {
          grid.text("**", x, y)
        } else if (p <= 0.05) {
          grid.text("*", x, y)
        }
      }
    })
  return(colocal_ht)
}


NN_substates_analysis <- function(spe, fib_state, samples_ord, 
                                  knn5_neigh_data) {
  
  # function to calculate mean log ratios and statistical significance 
  message("START ", fib_state, " ANALYSIS.")
  
  # PERFORM STATISTICAL ANALYSIS
  
  message("Assessing statistical significance.")
  
  # order of source cell labels
  substate_order <- c(fib_state, paste0("non", fib_state))
  # combine data from all timepoints
  knn5_df_list <- lapply(knn5_neigh_data, function(df) {
    df %>% as.data.frame() %>% 
      rownames_to_column("Source_labels")
  })
  knn5_df <- bind_rows(knn5_df_list) %>%
    # add sample name column
    mutate(Sample = factor(rep(names(knn5_neigh_data), each = 2), 
                           levels = samples_ord))
  
  # add source cell type and condition columns
  knn5_df <- knn5_df %>%
    mutate(Source_labels = factor(Source_labels, levels = substate_order),
           Condition = factor(str_split_i(Sample, "_", i = 1), 
                              levels = c("Sham", "Day1", "Day3", 
                                         "Day5", "Day7", "Day14"))) %>% 
    arrange(Sample, Source_labels)
  # rename the rows - add sample names
  rownames(knn5_df) <- paste0(knn5_df$Sample, "---", 
                              knn5_df$Source_labels)
  
  # merge each neighbour cell type counts from repeats in each condition
  knn5_df_merge <- reshape2::melt(knn5_df, variable.name = "Celltype", 
                                  value.name = "Counts") %>% # melt dataframe
    # add broad celltype level info
    mutate(CellGroup = ifelse(Celltype %in% levels(spe$Major.groups), 
                              as.character(Celltype), "Fibroblast")) %>%
    # group the data by condition and other variables
    group_by(Condition, Source_labels, Celltype, CellGroup) %>%
    # summarise counts in replicates
    summarise(tot_Counts = sum(Counts))
  
  cond_labels <- levels(knn5_df_merge$Condition) # unique conditions
  group_labels <- unique(knn5_df_merge$CellGroup) # unique cell groups
  fib_labels <- unique(as.character( # unique fibroblast substates
    knn5_df_merge[knn5_df_merge$CellGroup == "Fibroblast",]$Celltype))
  # define empty dataframe
  test_df <- data.frame(matrix( 
    nrow = 0, ncol = 5, dimnames = list(c(),
                                        c("Condition", "Celltype", 
                                          "Total_neighbours_detected",
                                          "X_statistic", "Pvalue"))))
  for (c in cond_labels) { # for each condition
    for (g in group_labels) { # for each cell group
      if (g != "Fibroblast") {
        # subset input data frame by condition
        sub_df <- knn5_df_merge[knn5_df_merge$Condition == c, ] %>%
          # create a label column for cell group and non-cell group class
          mutate(Cell_group_class = ifelse(CellGroup == g, g, 
                                           paste("non", g)))
        # build contingency table of sum of neighbour cell counts
        count_table <- as.matrix(
          xtabs(tot_Counts ~ Cell_group_class + Source_labels, 
                data = sub_df))
        # apply chi square test of significance
        chi_test <- chisq.test(count_table)
        # add total neighbours, chi square statistic, p-value to 
        # output data frame
        test_df <- rbind(
          test_df, 
          data.frame(Condition = c, Celltype = g, 
                     Total_neighbours_detected = sum(count_table),
                     X_statistic = chi_test$statistic, 
                     Pvalue = chi_test$p.value))
        show(paste("Finished", c, g))
      } else {
        for (f in fib_labels) { # for each fibroblast substate
          # subset input data frame by condition and fibroblast substate
          sub_fib_df <- knn5_df_merge[
            knn5_df_merge$Condition == c &
              knn5_df_merge$CellGroup == g, ] %>%
            # create label column for substate and non-substate classes
            mutate(Celltype_class = ifelse(Celltype == f, f, 
                                           paste("non", f)))
          # build contingency table of sum of neighbour cell counts
          count_table_fib <- as.matrix(
            xtabs(tot_Counts ~ Celltype_class + Source_labels, 
                  data = sub_fib_df))
          # apply chi square test of significance
          chi_test <- chisq.test(count_table_fib)
          # add total neighbours, chi square statistic, p-value to 
          # output data frame
          test_df <- 
            rbind(
              test_df,
              data.frame(
                Condition = c, Celltype = f,
                Total_neighbours_detected = sum(count_table_fib),
                X_statistic = chi_test$statistic,
                Pvalue = chi_test$p.value))
          show(paste("Finished", c, f))
        }
      }
    }
  }
  # calculate adjusted P-values
  test_df$adj_Pvalue <- p.adjust(test_df$Pvalue, method = "BH")
  rownames(test_df) <- NULL
  
  message("Statistical significance assessed.")
  
  
  # GENERATE MEAN LOG RATIO HEATMAP WITH STATISTICAL INFORMATION
  
  message("Generating heatmap.")
  
  # create adjusted p-value matrix
  test_df_mod <- test_df %>% 
    select(Condition, Celltype, adj_Pvalue) %>%
    pivot_wider(names_from = Condition, values_from = adj_Pvalue) %>%
    as.data.frame()
  rownames(test_df_mod) <- test_df_mod$Celltype
  p_mat <- test_df_mod %>% 
    select(-Celltype) %>%
    as.matrix()
  
  knn5_meanlogprop_ht <- 
    NN_meanlog_substate_heatmap(spe, knn5_df, p_mat, samples_ord,
                                substate_order)
  
  message("Heatmap generated.")
  
  message("FINISHED ", fib_state, " ANALYSIS.")
  
  return(knn5_meanlogprop_ht)
}



# PERFORM NEIGHBOURHOOD ANALYSES

# get matrix of cell type and substate associations per sample
knn5_neigh_data <- findNeighCount(
  spe = spe, nnCells = nnCells, samples = levels(spe$sample_id)[1:2], k = 5,
  celltype_labels = celltype_labels,
  index_labels = levels(spe$fibStates))

# assess each substates neighbours in each sample
ht_plt_list <- list()
for (fib in levels(spe$fibStates)) { 
  neigh_data <- setNames(lapply(names(knn5_neigh_data), function(s) {
    # extract neighbourhood association information for sample
    matX <- knn5_neigh_data[[s]]
    # generate neighbour association info for target substate vs non-target CFs
    mat <- rbind(matX[fib, ], colSums(matX[rownames(matX) != fib, ]))
    rownames(mat) <- c(fib, paste0("non", fib))
    return(mat)
  }), names(knn5_neigh_data))
  
  # assess and visualise substate neighbours
  ht_plt_list[[fib]] <- NN_substates_analysis(
    spe, fib_state = fib,
    samples_ord = levels(spe$sample_id), 
    knn5_neigh_data = neigh_data)
}




#### Extended Fig. 7B ####

markers_ss <- c("F-SH_sig", "F-SL_sig", "F-Trans_sig", "F-Act_sig", "MFC_sig") 

fib_df <- data.frame(x = spatialCoords(spe_yam)[,1],
                     y = spatialCoords(spe_yam)[,2],
                     samples = spe_yam$sample_group) %>%
  # sample factors
  mutate(samples = factor(samples, 
                          levels = c("WT_MI_day1_1", "WT_MI_day1_2", 
                                     "WT_MI_day1_3", "shCsrp3_MI_day1", 
                                     "Csrp3OE_MI_day1", "WT_MI_day14_1", 
                                     "WT_MI_day14_2", "WT_MI_day14_3", 
                                     "shCsrp3_MI_day14", "Csrp3OE_MI_day14")))
# define sample labels
cond_labels <- setNames(
  c("Day1\nWildtype_1", "Day1\nWildtype_2", "Day1\nWildtype_3", "Day1\nCsrp3 KO", 
    "Day1\nCsrp3 OE", "Day14\nWildtype_1", "Day14\nWildtype_2", 
    "Day14\nWildtype_3", "Day14\nCsrp3 KO", "Day14\nCsrp3 OE"), 
  levels(fib_df$samples))

sig_ss_plt_list <- setNames(lapply(markers_ss, function(gene) {
  if (gene == "MFC_sig") {
    MFC_genes = c("Fgl2", "Comp", "Dkk3", "Fmod", "Cthrc1",
                  "Col1a1", "Fn1", "Cfh", "Postn")
    fib_df <- fib_df %>%
      mutate(gene_exp = colMeans(logcounts(spe_yam)[MFC_genes, ]),
             gene_exp_z = scale(gene_exp)) %>%
      arrange(gene_exp_z)
  } else if (gene == "F-SH_sig") {
    FSH_genes = c("Ly6a", "Tmem100", "Efhd1", "Cd248", "Pi16", "Thy1")
    fib_df <- fib_df %>%
      mutate(gene_exp = colMeans(logcounts(spe_yam)[FSH_genes, ]),
             gene_exp_z = scale(gene_exp)) %>%
      arrange(gene_exp_z)
  } else if (gene == "F-SL_sig") {
    FSL_genes = c("Hsd11b1", "Cxcl14")
    fib_df <- fib_df %>%
      mutate(gene_exp = colMeans(logcounts(spe_yam)[FSL_genes, ]),
             gene_exp_z = scale(gene_exp)) %>%
      arrange(gene_exp_z)
  } else if (gene == "F-Act_sig") {
    FAct_genes = c("Cilp", "Meox1", "Col8a1", "Postn")
    fib_df <- fib_df %>%
      mutate(gene_exp = colMeans(logcounts(spe_yam)[FAct_genes, ]),
             gene_exp_z = scale(gene_exp)) %>%
      arrange(gene_exp_z)
  } else if (gene == "F-Trans_sig") {
    FTrans_genes = "Fgl2"
    fib_df <- fib_df %>%
      mutate(gene_exp = logcounts(spe_yam)[FTrans_genes, ],
             gene_exp_z = scale(gene_exp)) %>%
      arrange(gene_exp_z)
  }
  
  plt <- ggplot(fib_df, aes(x = x, y = y, color = gene_exp_z)) +
    geom_point(size = 1) +
    labs(color = gene) +
    facet_wrap(~samples, nrow = 2) +
    scale_color_gradientn(colors = c("#5a53a4", "#5ab0ab", "#fffbbb",
                                     "#f6794b", "#a10741"),
                          values = scales::rescale(c(-2, -1, 0, 1, 2)),
                          limits = c(-2, 2), 
                          oob = scales::squish) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.title = element_text(size = 12, colour = "white"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(size = 10, colour = "white"))
    return(plt)
}), markers_ss)




#### Extended Fig. 7C ####

# load processed data 
spe <- readRDS("./spe_processed_final.Rds")

# immune markers
specificImm_markers <- c(
  # Neutrophils
  "S100a8", "Fut4",
  # Inflammatory Macrophages
  "Cd14", "Ccr2", "Fcgr1",
  # Reparative Macrophages 
  "Cx3cr1", "C1qa", "C1qb", "Cd163", "Folr2", "Mgl2", "Csf1r", "Trem2")

gexp_marker <- assay(spe, "logcounts_area")[specificImm_markers, ] |> 
  as.matrix() |> t() |> as.data.frame()
# add cluster label to expression data frame
gexp_marker$cluster <- factor(as.character(spe$Immune_ct))
gexp_marker_melt <- reshape2::melt(gexp_marker, id.vars = "cluster", 
                                   variable.name = "gene", 
                                   value.name = "expression")
# summarize data
summary_data <- gexp_marker_melt %>%
  group_by(cluster, gene) %>%
  summarise(
    meanExpression = mean(expression),
    PerctageCells = sum(expression > 0) / n() * 100) %>%
  as.data.frame() %>%
  # scale gene expression
  group_by(gene) %>%
  mutate(exp_scaled = 
           (meanExpression - mean(meanExpression)) / sd(meanExpression),
         cluster = factor(cluster, 
                          levels = c("Neutrophil", "Inflammatory_Macrophage",
                                     "Reparative_Macrophage", "Other_Immune", 
                                     "Other")))
markerPlot <- summary_data %>%
  ggplot(aes(x = cluster, y = gene, size = PerctageCells, 
             color = exp_scaled)) +
  geom_point() +
  scale_size_continuous(limit = c(0, 100), range = c(0, 5)) + 
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") + 
  labs(color = "Scaled\nexpression", size = "Cells\ndetected") +
  scale_y_discrete(limits = rev) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.box = "vertical") 

