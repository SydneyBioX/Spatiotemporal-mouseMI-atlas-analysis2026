library(SpatialExperiment)
library(dplyr)
library(scater)
library(ggplot2)
library(tidyr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(spdep)

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




#### Fig. 8A ####

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




#### Fig. 8D ####

# subset to MYO and MFC cells
spe_sub <- spe[, spe$stabmapHarmonyDataset_lt_exact %in% c("MYO", "MFC")]

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

