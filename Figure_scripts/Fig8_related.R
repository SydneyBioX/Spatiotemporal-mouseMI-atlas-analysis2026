library(SpatialExperiment)
library(dplyr)
library(scater)
library(ggplot2)
library(tidyr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(spdep)
library(stringr)

set.seed(2026)


# load processed data 
spe <- readRDS("./spe_processed_final.Rds")

# organise samples
spe$sample_id <- factor(spe$sample_id, 
                        levels = c("Sham_1", "Sham_2", "Day1_1", "Day1_2", 
                                   "Day3_1", "Day3_2", "Day5_1", "Day5_2", 
                                   "Day7_1", "Day7_2", "Day14_1", "Day14_2"))

# relabel immune zone to IZ - since these are just immune cells found in IZ 
spe$major_zone_ed <- ifelse(spe$major_zone == "Immune", "IZ",
                            as.character(spe$major_zone))

# subset to fibroblasts
spe_fib <- spe[, spe$Major.groups == "Fibroblast"]

# consider only MYO and MFC fibrotic programs
spe_fib$fibro_labels <-
    ifelse(spe_fib$stabmapHarmonyDataset_lt_exact %in% c("MYO", "MFC"),
           spe_fib$fibrotic_score_class, "Other CF")
spe_fib$fibro_labels <- factor(
    spe_fib$fibro_labels,
    levels = c("Anti-fibrotic", "Pro-fibrotic", "Other CF"))

fibro_cols <- setNames(c("#71b587", "#ca6331", "grey10"),
                       c("Anti-fibrotic", "Pro-fibrotic", "Other CF"))




#### Fig. 8A & Extended Fig. 9A ####

fibrotic_sp_plt <- plotReducedDim(spe_fib, dimred = "spatial", 
                                  color_by = "fibro_labels", point_size = 0.05, 
                                  point_alpha = 1) + 
    facet_grid(spe_fib$replicate ~ spe_fib$condition, switch = "y") +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    scale_color_manual(values = fibro_cols) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(size = 14, colour = "white"),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.position = "none")




#### Fig. 8B ####

fibrotic_bar_plt <- table(spe_fib$sample_id, spe_fib$fibro_labels) %>%
    as.data.frame() %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = fibro_cols) +
    labs(y = "Number of cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_blank())




#### Fig. 8C ####

# subset data to specific fibroblast states in specific timepoints 
spe_fib_sub <- 
    spe_fib[, (spe_fib$condition %in% c("Day5", "Day7") & 
                   spe_fib$stabmapHarmonyDataset_lt_exact == "MYO") | 
                (spe_fib$condition == "Day14" & 
                     spe_fib$stabmapHarmonyDataset_lt_exact == "MFC")]
spe_fib_sub$condition <- factor(as.character(spe_fib_sub$condition), 
                                levels = c("Day5", "Day7", "Day14"))

# fibrotic program markers in ISS data
anti_markers <- c("Sfrp2", "Eln")
pro_markers <- c("Cilp", "Scx", "Thbs4", "Crlf1", "Col15a1")

exp_df <- data.frame(
    Sample = spe_fib_sub$sample_id, # add sample names
    Condition = spe_fib_sub$condition, # add timepoints
    Zone = spe_fib_sub$major_zone_ed, # add zone names
    # fibrotic program in specific substate
    FibState = paste(spe_fib_sub$fibrotic_score_class, 
                     as.character(spe_fib_sub$stabmapHarmonyDataset_lt_exact)),
    # anti-fibrotic signature
    Anti_sig = scale(colMeans(
        assay(spe_fib_sub, "logcounts_area")[anti_markers, ])),
    # pro-fibrotic signature
    Pro_sig = scale(colMeans(
        assay(spe_fib_sub, "logcounts_area")[pro_markers, ]))) %>%
    mutate(FibState = factor(FibState,
                             levels = c("Anti-fibrotic MYO", "Pro-fibrotic MYO",
                                        "Anti-fibrotic MFC", "Pro-fibrotic MFC")))

hm_df <- exp_df %>%
    # scaled mean signature per sample
    group_by(Sample, Condition, FibState, Zone) %>%
    summarise(Anti_sig = mean(Anti_sig, na.rm = TRUE),
              Pro_sig = mean(Pro_sig, na.rm = TRUE),
              .groups = "drop") %>%
    # scaled mean signature per replicate
    group_by(Condition, FibState, Zone) %>%
    summarise(Anti_sig = mean(Anti_sig, na.rm = TRUE),
              Pro_sig = mean(Pro_sig, na.rm = TRUE),
              .groups = "drop") %>%
    # convert signatures into rows
    pivot_longer(cols = c(Anti_sig, Pro_sig),
                 names_to = "Signature",
                 values_to = "Expression") %>%
    # convert fibrotic groups and zones to columns
    pivot_wider(names_from = Zone, values_from = Expression)

# anti-fibrotic information
anti_mat <- hm_df[hm_df$Signature == "Anti_sig", ] %>%
    select(Condition, FibState, RZ, BZ, IZ) %>%
    unite(rowname, Condition, FibState, sep = " ") %>%
    column_to_rownames("rowname") %>%
    as.matrix()
# pro-fibrotic information
pro_mat <- hm_df[hm_df$Signature == "Pro_sig", ] %>%
    select(Condition, FibState, RZ, BZ, IZ) %>%
    unite(rowname, Condition, FibState, sep = " ") %>%
    column_to_rownames("rowname") %>%
    as.matrix()
min_val <- min(anti_mat, pro_mat, na.rm = T)
max_val <- max(anti_mat, pro_mat, na.rm = T)

col_func <- colorRamp2(c(min_val, 0, max_val), 
                       c("#2166ACFF", "#F7F7F7FF", "#B2182BFF")) 
anti_hm <- Heatmap(anti_mat, 
                   name = paste0("Scaled\nexpression"),
                   col = col_func,
                   na_col = "grey95",
                   rect_gp = grid::gpar(col = "white", lwd = 1),
                   column_title = "Anti-fibrotic\nsignature",
                   column_title_gp = gpar(fontsize = 14),
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 12),
                   # show_column_names = FALSE,
                   column_names_gp = gpar(fontsize = 12),
                   column_names_side = "bottom")
pro_hm <- Heatmap(pro_mat, 
                  name = paste0("Scaled\nexpression"),
                  col = col_func,
                  na_col = "grey95",
                  rect_gp = grid::gpar(col = "white", lwd = 1),
                  column_title = "Pro-fibrotic\nsignature",
                  column_title_gp = gpar(fontsize = 14),
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 12),
                  column_names_gp = gpar(fontsize = 12),
                  column_names_side = "bottom")
gene_zone_hm <- anti_hm + pro_hm




#### Fig. 8D & Extended Fig. 9C ####

# ------------ #

##### For Fig. 8D ####
# subset to MYO and MFC cells
spe_sub <- spe[, spe$stabmapHarmonyDataset_lt_exact %in% c("MYO", "MFC")]

##### For Extended Fig. 9C ####
# subset to fibroblast cells
spe_sub <- spe[, spe$Major.groups == "Fibroblast"]

# ------------ #

# binarise fibrotic program labels
# 1 = Anti-fibrotic; 0 = Pro-fibrotic
spe_sub$Fib_group <- as.numeric(spe_sub$fibrotic_score_class == "Anti-fibrotic")
samples <- levels(spe$sample_id)


# MORAN'S I RANDOMNESS TEST

# intialise data frame to collect relevant information
moran_rand_df <- data.frame(matrix(nrow = 0, ncol = 5))

# assess Moran's I across different neighbourhood sizes
for (i in c(6, 10, 20, 50, 100)) {
    moran_rand_out <- lapply(samples, function(s) {
        # subset to sample
        speX <- spe_sub[, spe_sub$sample_id == s]
        # add check - if cells in sample < 100 then Moran's I is NA
        if (ncol(speX) < 100) {
            df <- data.frame(Sample = s,
                             moran_I = NA,  
                             p_val = NA)
            print(paste(s, ncol(speX)))
            return(df)
        }
        # get cell locations
        coords <- reducedDim(speX, "spatial")
        # identify i neighbours
        knn_obj <- knearneigh(coords, k = i)
        # convert knn_obj to list
        nb <- knn2nb(knn_obj)
        # add spatial weights
        lw <- nb2listw(nb, style = "W")
        # perform Moran’s I permutation test to account for 
        # randomness and uncertainity
        moran <- moran.mc(speX$Fib_group, lw, nsim = 1000)
        df <- data.frame(Sample = s,
                         moran_I = moran$statistic,  
                         p_val = moran$p.value)
        print(paste(s, ncol(speX)))
        return(df)
    })
    # combine Moran's I output from all samples
    moran_rand <- do.call(rbind, moran_rand_out)
    # calculate BH-corrected P-values
    moran_rand$adj_pval <- p.adjust(moran_rand$p_val, method = "BH")
    moran_rand$k <- i # add neighbourhood size
    moran_rand_df <- rbind(moran_rand_df, moran_rand)
}

moran_rand_df <- moran_rand_df %>%
    mutate(k = factor(k, levels = c(6, 10, 20, 50, 100)),
           Sample = factor(Sample, levels = levels(spe$sample_id)),
           moran_I = ifelse(is.na(moran_I), 0, moran_I),
           adj_pval = ifelse(is.na(adj_pval), 1, adj_pval),
           Significant = factor(ifelse(adj_pval <= 0.05, "Yes", "No"), 
                                levels = c("Yes", "No")))


# MORAN'S I ROBUSTNESS TEST

# intialise data frame to collect relevant information
moran_rob_df <- data.frame(matrix(nrow = 0, ncol = 5))
# assess Moran's I across different neighbourhood sizes
for (i in c(6, 10, 20, 50, 100)) { # for each neighbourhood size 
    moran_rob_out <- lapply(samples, function(s) {
        # subset to sample
        speX <- spe_sub[, spe_sub$sample_id == s]
        # add check - if cells in sample < 100 then Moran's I is NA
        if (ncol(speX) < 100) {
            df <- data.frame(Sample = s,
                             moran_base = NA,
                             moran_boot_mean = NA,
                             moran_boot_sd = NA)
            print(paste(s, ncol(speX)))
            return(df)
        }
        # get cell locations
        coords <- reducedDim(speX, "spatial")
        # identify i neighbours
        knn_obj <- knearneigh(coords, k = i)
        # convert knn obj to list
        nb <- knn2nb(knn_obj)
        # add spatial weights
        lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
        # calculate base Moran's I without bootstrapping
        fib_bin <- speX$Fib_group
        moran_base <- moran(fib_bin, lw, n = ncol(speX), S0 = Szero(lw))$I
        # perform bootstrapped Moran’s I calculation for robustness without 
        # replacement
        boot_vals <- replicate(100, { # number of bootstraps
            # 50% sub-sampling
            idx <- sample(seq_along(fib_bin), size = floor(0.5 * length(fib_bin)), 
                          replace = F)
            x_boot <- fib_bin[sort(idx)]
            # re-index spatial structure - use the underlying spatial structure
            # while re-labelling the signal, and then assess how Moran's I varies
            nb_boot <- subset(nb, seq_along(fib_bin) %in% idx)
            lw_boot <- nb2listw(nb_boot, style = "W", zero.policy = TRUE)
            moran(x_boot, lw_boot, n = length(x_boot), S0 = Szero(lw_boot))$I })
        
        df <- data.frame(Sample = s,
                         moran_base = moran_base,
                         moran_boot_mean = mean(boot_vals),
                         moran_boot_sd = sd(boot_vals))
        print(paste(s, ncol(speX)))
        return(df)
    })
    # combine Moran's I output from all samples
    moran_rob <- do.call(rbind, moran_rob_out)
    moran_rob$k <- i # add neighbourhood size
    moran_rob_df <- rbind(moran_rob_df, moran_rob)
}

moran_rob_df <- moran_rob_df %>%
    mutate(k = factor(k, levels = c(6, 10, 20, 50, 100)),
           Sample = factor(Sample, levels = levels(spe$sample_id)),
           moran_base = ifelse(is.na(moran_base), 0, moran_base),
           moran_boot_mean = ifelse(is.na(moran_boot_mean), 0, moran_boot_mean),
           moran_boot_sd = ifelse(is.na(moran_boot_sd), 0, moran_boot_sd))


# VISUALISE MORAN'S I OUTPUT

moran_df <- left_join(moran_rand_df, moran_rob_df, by = c("Sample", "k"))
k_cols <- setNames(c("#355070", "#6d597a", "#b56576", "#e56b6f", "#eaac8b"),
                   c("6", "10", "20", "50", "100"))
morani_plot <- ggplot(moran_df, aes(x = Sample, y = moran_I, group = k)) +
    # add robustness output
    geom_ribbon(aes(ymin = moran_boot_mean - moran_boot_sd,
                    ymax = moran_boot_mean + moran_boot_sd, fill = k),
                alpha = 0.4) + 
    # add base Moran's I values
    geom_line(aes(colour = k), linewidth = 0.5) +
    # add permutation test output
    geom_point(aes(shape = Significant, fill = Significant), 
               size = 2.5, stroke = 0.5, color = "black",
               fill = ifelse(moran_df$Significant == "Yes", "red", "black")) + 
    scale_color_manual(values = k_cols) +
    scale_shape_manual(values = c("Yes" = 21, "No" = 22)) +
    scale_fill_manual(values = k_cols) +
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3))) +
    labs(y = "Moran's I", color = "Neighbourhood\nsize", 
         fill = "Neighbourhood\nsize") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))




#### Fig. 8E & Extended Fig. 9E ####


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

# add cell ids to fibrotic class labels
names(spe$fibrotic_score_class) <- colnames(spe) 
# add cell ids to fibroblast substates
names(spe$fibStates) <- colnames(spe) 
# add cell ids to broad cell groups
names(spe$Major.groups) <- colnames(spe) 
# add cell ids to immune cell groups
names(spe$Immune_ct) <- colnames(spe)

# creating fibrotic class labels
fibrotic_labels <- 
    ifelse(spe$fibStates %in% c("MYO", "MFC"), 
           paste0(spe$fibStates, "_", 
                  spe$fibrotic_score_class), NA)
myo_grps <- unique(fibrotic_labels[!is.na(fibrotic_labels)])

# creating neighbouring cell labels
cell_annots <- case_when(
    # Fibroblast neighbours
    spe$Major.groups == "Fibroblast" & spe$fibStates %in% c("MYO", "MFC") ~
        paste0(spe$fibStates, "_", spe$fibrotic_score_class),
    spe$Major.groups == "Fibroblast" & !(spe$fibStates %in% c("MYO", "MFC")) ~
        as.character(spe$fibStates),
    # Immune neighbours
    spe$Major.groups == "Immune" ~ spe$Immune_ct,
    # all other neighbours
    TRUE ~ spe$Major.groups)
names(cell_annots) <- colnames(spe)
celltype_labels <- c(unique(cell_annots[!is.na(cell_annots)]), # unique labels
                     c("MYO", "MFC", "Fibroblast", "Immune")) # add other broad groups

# getting CF substate labels - use to add additional CF column in final mat
fibStates <- as.character(unique(na.omit(spe$fibStates)))
# getting CF substate labels - use to add additional Immune column in final mat
immCelltypes <- as.character(unique(spe$Immune_ct[spe$Immune_ct != "Other"]))



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

# DETECT NEIGHBOURS

# function for identifying k nearest neighbours of fibrotic cells in each sample
findNeighCount <- function(spe, samples, k, cell_annots, myo_grps, 
                           celltype_labels, fibStates, immCelltypes) {
    
    # function to identify the cell type or substate labels of 5-nn
    
    coloc_matList <- setNames(lapply(samples, function(s) {
        # subset sample
        speX <- spe[, spe$sample_id == s] 
        # initialise neighbour count matrix 
        coloc_mat <- matrix(0, nrow = length(myo_grps),
                            ncol = length(celltype_labels), 
                            dimnames = list(myo_grps, celltype_labels))
        
        # extract cell ids of MYO and MFC fibrotic cells
        fibrotic_cells <- 
            speX[, speX$fibStates %in% c("MYO", "MFC")]$fibrotic_score_class |> 
            names()
        print(length(fibrotic_cells))
        
        for(i in fibrotic_cells) { # for each fibrotic cell in the sample
            
            # identify fibrotic label
            source_grp <- paste0(speX$fibStates[i], "_", 
                                 speX$fibrotic_score_class[i])
            
            # select nearest neighbours
            neigh_ids <- nnCells[i, 1:k]
            # identify neighbour cell type annotation
            neigh_grps <- cell_annots[neigh_ids]
            # add CF group label for every CF substate neighbour
            if (length(intersect(neigh_grps, fibStates)) != 0) {
                fib_num <- length(intersect(neigh_grps, fibStates))
                neigh_grps <- c(neigh_grps, rep("Fibroblast", fib_num))
            }
            # add MYO group label for every MYO substate neighbours
            if(length(intersect(str_split_i(neigh_grps, "_", i = 1), 
                                "MYO")) != 0) {
                myo_num <- length(intersect(str_split_i(neigh_grps, "_", i = 1), 
                                            "MYO"))
                neigh_grps <- c(neigh_grps, rep("MYO", myo_num))
            }
            # add MFC group label for every MFC substate neighbours
            if(length(intersect(str_split_i(neigh_grps, "_", i = 1), 
                                "MFC")) != 0) {
                mfc_num <- length(intersect(str_split_i(neigh_grps, "_", i = 1), 
                                            "MFC"))
                neigh_grps <- c(neigh_grps, rep("MFC", mfc_num))
            }
            # add Immune group label for every immune cell type neighbour
            if (length(intersect(neigh_grps, immCelltypes)) != 0) {
                immune_num <- length(intersect(neigh_grps, immCelltypes))
                neigh_grps <- c(neigh_grps, rep("Immune", immune_num))
            }
            
            # calculate the number of neighbouring cells belonging to each cell type
            coloc_count <- table(factor(neigh_grps, levels = celltype_labels))
            # convert total counts to binary detection calls 
            # detected = 1, not detected = 0
            coloc_binary <- setNames(as.integer(coloc_count > 0), names(coloc_count))
            # sum neighbour counts for each sample
            coloc_mat[source_grp, ] <- coloc_mat[source_grp, ] + coloc_binary
        }
        print(s)
        return(coloc_mat)
    }), samples)
    return(coloc_matList)
}

# detect neighbours in each sample
knn5_data <- findNeighCount(spe, samples = levels(spe$sample_id), k = 5, 
                            cell_annots, myo_grps, celltype_labels, fibStates, 
                            immCelltypes)

# combine data from all timepoints
knn5_df_list <- lapply(knn5_data, function(df) {
    df %>% as.data.frame() %>% 
        rownames_to_column("Source_labels")
})
knn5_df <- bind_rows(knn5_df_list) %>%
    # add sample name column
    mutate(Sample = factor(rep(names(knn5_data), each = 4), 
                           levels = levels(spe$sample_id))) %>%
    # add source celltype and condition columns
    mutate(Source_labels = factor(Source_labels),
           Condition = factor(str_split_i(Sample, "_", i = 1), 
                              levels = c("Sham", "Day1", "Day3", 
                                         "Day5", "Day7", "Day14"))) %>% 
    arrange(Sample, Source_labels)
# rename the rows - add sample names
rownames(knn5_df) <- paste0(knn5_df$Sample, "---", 
                            knn5_df$Source_labels)



# STATISTICAL ANALYSIS

perform_chiqtest <- function(knn5_df) {
    # merge each neighbour cell type counts from repeats in each condition
    knn5_df_merge <- reshape2::melt(knn5_df, variable.name = "Celltype", 
                                    value.name = "Counts") %>%
        # add celltype level info
        mutate(CF_type = str_split_i(Source_labels, "_", i = 1),
               Fibrotic_type = str_split_i(Source_labels, "_", i = 2)) %>%
        # group the data by condition and other variables
        group_by(Condition, Source_labels, Celltype, CF_type, Fibrotic_type) %>%
        # summarise counts in replicates
        summarise(tot_Counts = sum(Counts))
    
    index_labels <- unique(knn5_df_merge$CF_type)
    cond_labels <- levels(knn5_df_merge$Condition) # unique conditions
    group_labels <- levels(knn5_df_merge$Celltype) # unique neighbour celltypes
    test_df <- data.frame(matrix( # defining empty dataframe
        nrow = 0, ncol = 6, 
        dimnames = list(c(), c("Condition", "I_celltype", "N_celltype", 
                               "Total_neighbours_detected", "X_statistic", 
                               "Pvalue"))))
    for (s in index_labels){ # for each index cell type
        for (c in cond_labels) { # for each condition
            for (g in group_labels) { # for each neighbour cell group
                # subset input data frame by index celltype and condition
                print(paste(s, c, g))
                sub_df <- knn5_df_merge[knn5_df_merge$CF_type == s & 
                                            knn5_df_merge$Condition == c, ] %>%
                    # create a label column for cell group and non-cell group class
                    mutate(Cell_group_class = ifelse(Celltype == g, g, 
                                                     paste("non", g))) %>%
                    as.data.frame()
                print(dim(sub_df))
                # build contingency table of sum of neighbour cell counts
                count_table <- as.matrix(
                    xtabs(tot_Counts ~ Cell_group_class + Fibrotic_type, 
                          data = sub_df))
                if(sum(count_table) != 0) {
                    # apply chi square test of significance
                    chi_test <- chisq.test(count_table)
                    # add total neighbours, chi square statistic, p-value to 
                    # output data frame
                    test_df <- rbind(
                        test_df, 
                        data.frame(Condition = c, I_celltype = s, N_celltype = g, 
                                   Total_neighbours_detected = sum(count_table),
                                   X_statistic = chi_test$statistic, 
                                   Pvalue = chi_test$p.value))
                    show(paste("Finished", c, g))
                } else {
                    test_df <- rbind(
                        test_df, 
                        data.frame(Condition = c, I_celltype = s, N_celltype = g, 
                                   Total_neighbours_detected = sum(count_table),
                                   X_statistic = NA, 
                                   Pvalue = NA))
                    show(paste("Finished", c, g))
                }}}}
    rownames(test_df) <- NULL
    # calculate adjusted P-values
    test_df$adj_Pvalue <- p.adjust(test_df$Pvalue, method = "BH")
    # create adjusted p-value matrix
    test_df_mod <- test_df %>% 
        mutate(Sample = paste0(Condition, "-", I_celltype)) %>%
        select(Sample, N_celltype, adj_Pvalue) %>%
        pivot_wider(names_from = Sample, values_from = adj_Pvalue) %>%
        as.data.frame()
    rownames(test_df_mod) <- test_df_mod$N_celltype
    p_mat <- test_df_mod %>% 
        select(-N_celltype) %>%
        as.matrix()
    return(list(test_df = test_df,
                p_mat = p_mat))
}

# run test on broad cell types
knn5_df_ct <- knn5_df %>%
    select(Source_labels, Sample, Condition, unique(spe$Major.groups))
test_ct <- perform_chiqtest(knn5_df_ct)
# run test on CF substates
knn5_df_fib <- knn5_df %>%
    select(Source_labels, Sample, Condition, levels(spe$fibStates))
test_fib <- perform_chiqtest(knn5_df_fib)
# run test on MFC cells
knn5_df_mfc <- knn5_df %>%
    select(Source_labels, Sample, Condition, `MFC_Anti-fibrotic`, 
           `MFC_Pro-fibrotic`)
test_mfc <- perform_chiqtest(knn5_df_mfc)
# run test on MYO cells
knn5_df_myo <- knn5_df %>%
    select(Source_labels, Sample, Condition, `MYO_Anti-fibrotic`, 
           `MYO_Pro-fibrotic`)
test_myo <- perform_chiqtest(knn5_df_myo)
# run test on CF substates
knn5_df_immune <- knn5_df %>%
    select(Source_labels, Sample, Condition, immCelltypes)
test_immune <- perform_chiqtest(knn5_df_immune)

# gather p-values from all tests
p_mat <- rbind(test_ct$p_mat, test_fib$p_mat, test_mfc$p_mat, test_myo$p_mat, 
               test_immune$p_mat)


# CALCULATE MEAN LOG RATIO

# row arrangement for heatmap
rownames_arrange <- c("MYO_Anti-fibrotic", "MYO_Pro-fibrotic",
                      "MFC_Anti-fibrotic", "MFC_Pro-fibrotic",
                      levels(spe$fibStates),
                      "Neutrophil", "Inflammatory_Macrophage",
                      "Reparative_Macrophage", "Other_Immune",
                      "Cardiomyocyte", "Stressed cardiomyocyte",
                      "Immune", "Endothelial", "Fibroblast")

# reorganise p_mat rows
p_mat <- p_mat[rownames_arrange, ]

# remove sample and source label names
knn5_df[, c("Sample", "Source_labels", "Condition")] <- NULL

# calculate neighbor proportions
coloc_prop <- as.matrix(knn5_df) / rowSums(as.matrix(knn5_df))

# get mean log ratio of substate vs non-substate cells
coloc_mat <- coloc_prop %>%
    as.data.frame() %>% # convert to data frame
    # add column for substate and non-substate fibroblasts
    mutate(Source_label = str_split_i(rownames(coloc_prop), "---", i = 2),
           Sample = factor(str_split_i(rownames(coloc_prop), "---", i = 1),
                           levels = levels(spe$sample_id)),
           Fibrotic_group = str_split_i(Source_label, "_", i = 2),
           Fibroblast_Substate = str_split_i(Source_label, "_", i = 1),
           Source_label = NULL) %>%
    # melt dataframe
    pivot_longer(!c(Sample, Fibrotic_group, Fibroblast_Substate), 
                 names_to = "Neighbouring_celltypes", 
                 values_to = "counts") %>%
    # elongate dataframe to compare between fibrotic classes
    pivot_wider(names_from = Fibrotic_group, values_from = counts) %>%
    # calculate and add log ratio of substate vs non-substate 
    mutate(logratio = log2((`Anti-fibrotic` + 1) / 
                               (`Pro-fibrotic` + 1)),
           `Anti-fibrotic` = NULL, 
           `Pro-fibrotic` = NULL) %>%
    # add Day and repeat columns
    mutate(Day = str_split_i(Sample, "_", i = 1),
           Repeat = str_split_i(Sample, "_", i = 2),
           Sample = NULL) %>%
    # elongate data frame to compare between repeats
    pivot_wider(names_from = Repeat, values_from = logratio) %>%
    # order data by fibroblast substate
    arrange(Fibroblast_Substate) %>%
    # calculate mean of logRatio in sample repeats
    mutate(MeanLogRatio = (`1` + `2`) / 2,
           Sample = paste0(Day, "-", Fibroblast_Substate),
           `1` = NULL, `2` = NULL, Day = NULL, Fibroblast_Substate = NULL) %>%
    # elongate dataframe to sample x neighbourhood with mean logratio values
    pivot_wider(names_from = Sample, 
                values_from = MeanLogRatio) %>%
    as.data.frame() %>% # converting to data frame
    # add timepoint and fibroblast substate names as rownames
    column_to_rownames("Neighbouring_celltypes") %>%
    as.matrix()



# HEATMAP OF MEAN LOG RATIOS WITH STATISTICAL INFERENCE

# count substate cells in each condition   
antifib_count_df <- 
    table(spe$condition, spe$fibStates,
          spe$fibrotic_score_class)[,c("MFC", "MYO"), "Anti-fibrotic"] %>%
    as.data.frame()
antifib_count <- setNames(antifib_count_df$Freq, 
                          paste0(antifib_count_df$Var1, "-", 
                                 antifib_count_df$Var2))

# add top annotation
col_annot <- do.call(HeatmapAnnotation, setNames(
    list(anno_barplot(antifib_count, gp = gpar(fill = "black"))),
    paste0("Total\nanti-fibrotic\ncells")))

# define color range
col_func <- colorRamp2(c(min(coloc_mat, na.rm = TRUE), 0, 
                         max(coloc_mat, na.rm = TRUE)), 
                       c("#ca6331", "white", "#71b587"))

# match ratio and p-value matrix colnames
if (identical(colnames(p_mat), colnames(coloc_mat)) != TRUE) {
    print("Different colnames.")
    p_mat <- p_mat[colnames(coloc_mat), ]
} else {
    # keep only condition name
    colnames(p_mat) <- str_split_i(colnames(p_mat), "-", i = 1)
    colnames(coloc_mat) <- str_split_i(colnames(coloc_mat), "-", i = 1)
}
coloc_mat <- coloc_mat[rownames_arrange,]

# create heatmap
colocal_ht_all <- Heatmap(
    coloc_mat, 
    name = paste0("Log ratio of\nAnti-/Pro-\nfibrotic cells"),
    col = col_func,
    na_col = "grey95",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    border_gp = grid::gpar(col = "black", lty = 2, lwd = 1),
    rect_gp = grid::gpar(col = "white", lwd = 1),
    column_title_gp = gpar(fontsize = 14),
    row_title = NULL,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 12),
    row_split = factor(c(rep("MYO", 2), rep("MFC", 2), rep("CF substate", 11),
                         rep("Immune", 4), rep("Celltype", 5)),
                       levels = c("MYO", "MFC", "CF substate", "Immune", "Celltype")),
    column_names_gp = gpar(fontsize = 12),
    column_names_side = "top",
    column_split = factor(rep(c("MFC", "MYO"), each = 6), 
                          levels = c("MYO", "MFC")),
    
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



# BAR PLOTS of MEAN PROPORTIONS

knn5_new_df <- bind_rows(knn5_df_list) %>%
    # add sample name column
    mutate(Sample = rep(names(knn5_data), each = 4),
           CF_state = str_split_i(Source_labels, "_", i = 1),
           Fibrotic_group = str_split_i(Source_labels, "_", i = 2),
           Replicate = str_split_i(Sample, "_", i = 2),
           Day = str_split_i(Sample, "_", i = 1))

# cell group information
celltypes <- c("Cardiomyocyte", "Stressed cardiomyocyte", "Immune", 
               "Endothelial", "Fibroblast")
cellstates <- c("F-SH", "F-SL", "F-Trans", "F-WntX_Wif1", "F-Act", "IR", "F-CI", 
                "F-Cyc", "MYO", "MFC", "F-IFNS")
immunecells <- c("Neutrophil", "Inflammatory_Macrophage", 
                 "Reparative_Macrophage", "Other_Immune")
myocells <- c("MYO_Anti-fibrotic", "MYO_Pro-fibrotic")
mfccells <- c("MFC_Anti-fibrotic", "MFC_Pro-fibrotic")


df_prop <- knn5_new_df %>%
    pivot_longer(cols = -c(Source_labels, CF_state, Fibrotic_group, Sample,
                           Replicate, Day), names_to = "Neighbor",
                 values_to = "Count") %>%
    mutate(Cell_group = case_when(Neighbor %in% celltypes ~ "Cardiac cell types",
                                  Neighbor %in% cellstates ~ "Fibroblasts",
                                  Neighbor %in% immunecells ~ "Immune cells",
                                  Neighbor %in% myocells ~ "MYO",
                                  Neighbor %in% mfccells ~ "MFC",
                                  TRUE ~ NA)) %>%
    
    # calculate proportions first
    group_by(CF_state, Cell_group, Fibrotic_group, Day, Replicate) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    # average proportions over replicates
    group_by(CF_state, Cell_group, Fibrotic_group, Day, Neighbor) %>%
    summarize(Avg_prop = mean(Proportion),
              Tot_counts = mean(Count),
              .groups = "drop") %>%
    # invert proportions and counts for pro-fibrotic CFs for plotting
    mutate(Prop_invert = ifelse(Fibrotic_group == "Pro-fibrotic", -Avg_prop,
                                Avg_prop),
           Count_invert = ifelse(Fibrotic_group == "Pro-fibrotic", -Tot_counts,
                                 Tot_counts))


# plot all cell groups in select timepoints
df_subset <- df_prop[(df_prop$CF_state == "MYO" & 
                          df_prop$Day %in% c("Day5", "Day7")) | 
                         (df_prop$CF_state == "MFC" & 
                              df_prop$Day %in% c("Day14")), ] %>%
    mutate(Neighbor = factor(Neighbor, levels = c(myocells, mfccells, cellstates, 
                                                  immunecells, celltypes)),
           Cell_group = factor(Cell_group, 
                               levels = c("MYO", "MFC", "Fibroblasts", 
                                          "Immune cells", "Cardiac cell types")),
           Day = factor(Day, levels = c("Day5", "Day7", "Day14")),
           CF_state = factor(CF_state, levels = c("MYO", "MFC")))

main_plot <- ggplot(df_subset, aes(x = Neighbor, y = Avg_prop, 
                                   fill = Fibrotic_group)) +
    geom_col(position = "dodge", color = "black", linewidth = 0.2) +
    facet_grid(CF_state + Day ~ Cell_group, scales = "free_x", space = "free_x",
               switch = "y") +
    scale_fill_manual(values = c("Anti-fibrotic" = "#71b587", 
                                 "Pro-fibrotic" = "#ca6331")) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.1) +
    labs(y = "Mean proportion of neighbours",
         x = "Neighbours",
         fill = "Fibrotic\nprogram") +
    theme_classic() +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 12))




#### Extended Fig. 9B ####

fibrotic_bar_plt <- table(spe_fib$sample_id, spe_fib$fibro_labels) %>%
    as.data.frame() %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = fibro_cols) +
    labs(y = "Proportion of cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_blank())




#### Extended Fig. 9D ####

df_fibsclass <- data.frame(Sample = spe_fib$sample_id, 
                           Fibroblast = spe_fib$stabmapHarmonyDataset_lt_exact,
                           FibroticClass = spe_fib$fibrotic_score_class) %>%
    # identify all combinations of sample and fibroblasts 
    # add missing combinations and allocate NA value 
    # to corresponding FibroticClass
    complete(Sample, Fibroblast, fill = list(FibroticClass = NA))

summary_df <- df_fibsclass %>%
    group_by(Fibroblast, Sample) %>%
    # total cells of fibroblast subtype in sample
    # handle NA values
    summarise(n = sum(!is.na(FibroticClass)), 
              # total antifibrotic cells in fibroblast subtype and sample
              k = sum(FibroticClass == "Anti-fibrotic", na.rm = TRUE),
              # if no fibroblast in sample, percentage is zero
              perc = ifelse(n > 0, (k / n) * 100, 0),
              # if no fibroblast in sample, binomial sd is NA
              sd = ifelse(n > 0, sqrt(perc * (100 - perc) / n), NA),
              .groups = "drop")

fibrotic_barplt <- ggplot(summary_df, aes(x = Sample, y = perc)) +
    geom_col(fill = "#71b587") +
    geom_errorbar(aes(ymin = perc - sd, ymax = perc + sd), width = 0.2) +
    facet_wrap(~ Fibroblast, scales = "free_x", nrow = 1) +
    labs(x = "Sample", y = "% of Antifibrotic cells") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.background = element_blank(),
          strip.text = element_text(size = 14))

