library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggridges)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(tibble)
library(tidyr)
library(scales)

set.seed(2026)


# load processed data 
spe <- readRDS("./spe_processed_final.Rds")

# organise samples
spe$sample_id <- factor(spe$sample_id, 
                        levels = c("Sham_1", "Sham_2", "Day1_1", "Day1_2", 
                                   "Day3_1", "Day3_2", "Day5_1", "Day5_2", 
                                   "Day7_1", "Day7_2", "Day14_1", "Day14_2"))

# relabel immune zone to IZ - since these are just immune cells found in IZ 
spe$major_zone_ed <- factor(ifelse(spe$major_zone == "Immune", "IZ", 
                                   as.character(spe$major_zone)),
                            levels = c("RZ", "BZ", "IZ"))

# subset to fibroblasts
spe_fib <- spe[, spe$Major.groups == "Fibroblast"]

fib_cols <- setNames(c("#4CD02D", "#E88321", "#DE331D", "#E9EA38", "#B7C1CD",
                       "#5CA6C9", "#F7CADE", "#479B74", "#ADF0E9", "#C850BA", 
                       "#4559B0"),
                     c("F-Trans", "F-SL", "F-Act", "F-SH", "MFC", "MYO", "F-IFNS",
                       "F-WntX (Wif1)", "F-CI", "IR", "F-Cyc"))



# subset data to fibroblasts in BZ and IZ for cell density plots
# exclude Sham samples, since they have no zones
spe_sub <- spe[, spe$Major.groups == "Fibroblast" & spe$major_zone != "RZ" &
                   spe$condition != "Sham"]

# combine two sets of distances - distance of BZ cells to IZ and 
# distance of IZ cells to BZ
# also, adjust distance of BZ cells to IZ so that BZ-IZ boundary is 0
spe_sub$dist_in_zone <- 
    ifelse(spe_sub$major_zone == "BZ", -spe_sub$distance_to_IZ, 
           spe_sub$distance_to_BZ)

# exclude F-WntX cells with false WntX signal (only keep Wif1+ F-WntX cells) 
# add Wif1 expression per cell
spe_sub$Wif1_exp <- as.numeric(assay(spe_sub, "logcounts_area")["Wif1", ])
# add Wif1 detection status
spe_sub$Wif1_detetct <- ifelse(spe_sub$Wif1_exp == 0, "No", "Yes")
# relabel F-WntX to F-WntX (Wif1), if the cell has Wif1 expression, else NA
spe_sub$fibStates <- ifelse(spe_sub$stabmapHarmonyDataset_lt_exact == "F-WntX",
                            ifelse(spe_sub$Wif1_detetct == "Yes", "F-WntX (Wif1)",
                                   NA),
                            as.character(spe_sub$stabmapHarmonyDataset_lt_exact))
spe_sub$fibStates <- factor(spe_sub$fibStates, 
                            levels = c("F-SH", "F-SL", "F-Trans", "F-WntX (Wif1)", 
                                       "F-Act", "IR", "F-CI", "F-Cyc", "MYO",
                                       "MFC", "F-IFNS"))
# remove F-WntX cells with no Wif1 expression
spe_sub <- spe_sub[, !is.na(spe_sub$fibStates)]

dense_df <- data.frame( # building data frame
    distance = spe_sub$dist_in_zone, # adding distance
    fibroblast = spe_sub$fibStates, # adding fib substate
    sample = as.character(spe_sub$sample_id)) %>% # adding sample id
    mutate(sample = factor(sample, # converting to factor
                           levels = rev(c("Day1_1", "Day1_2", "Day3_1", "Day3_2",
                                          "Day5_1", "Day5_2", "Day7_1", "Day7_2",
                                          "Day14_1", "Day14_2"))))

#### Fig. 6A ####

cell_density_plt <- dense_df %>%
    group_by(sample, fibroblast) %>%
    filter(n() > 100) %>% # keeping only subgroups with more than 100 cells
    do({ # calculating cell densities for subgroups to ensure comparability
        d <- density(.$distance, bw = 300) 
        data.frame(x = d$x, y = d$y * nrow(.)) 
    }) %>% 
    ungroup() %>%
    # plotting cell densities
    ggplot(aes(x = x, y = sample, height = y, color = fibroblast, 
               fill = fibroblast)) +
    geom_vline(xintercept = 0, col = "grey50", linetype = "dashed", linewidth = 0.2) +
    geom_density_ridges(stat = "identity", scale = 1, alpha = 0.7, linewidth = 0.2) +
    facet_wrap(~fibroblast, nrow = 1, drop = FALSE) +
    scale_color_manual(values = fib_cols) +
    scale_fill_manual(values = fib_cols) +
    xlim(-15000, 15000) +
    scale_x_continuous(breaks = c(-15000, -7500, 0, 7500, 15000)) +
    labs(x = "Euclidean spatial distance", y = "Cell density",
         color = NULL, fill = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 10, angle = 60, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 12, face = "bold"))




#### Extended Fig. 6B ####

cell_density_plt2 <- dense_df %>%
    group_by(sample, fibroblast) %>%
    do({ # calculating cell densities for subgroups to ensure comparability
        d <- density(.$distance, bw = 300) 
        data.frame(x = d$x, y = d$y * nrow(.)) 
    }) %>% 
    ungroup() %>%
    # plotting cell densities
    ggplot(aes(x = x, y = sample, height = y, color = fibroblast, 
               fill = fibroblast)) +
    geom_vline(xintercept = 0, col = "grey50", linetype = "dashed", linewidth = 0.2) +
    geom_density_ridges(stat = "identity", scale = 1, alpha = 0.5, linewidth = 0.2) +
    scale_color_manual(values = fib_cols) +
    scale_fill_manual(values = fib_cols) +
    labs(x = "Distance across border zone and injury zone", y = "Cell density",
         color = "Fibroblast\nsubstate", fill = "Fibroblast\nsubstate") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title = element_blank(),
          legend.title = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.placement = "none")




#### Fig. 6B ####

markers_ss <- c("Pi16", "Ly6a", "Fgl2", "Wif1", "Postn", "Mt1", "Timp1", "Cdc20", 
                "Ddah1", "Comp", "Cxcl10")

# extract log normalised expression of marker genes
gexp_marker <- assay(spe_fib, "logcounts_area")[markers_ss, ] |> 
    as.matrix() |> t() |> as.data.frame()
# add cluster label to expression data frame
gexp_marker$group <- spe_fib$stabmapHarmonyDataset_lt_exact
# add zone label to expression data frame
gexp_marker$zone <- spe_fib$major_zone_ed

summary_data <- reshape2::melt(gexp_marker, id.vars = c("group", "zone"), 
                               variable.name = "gene",
                               value.name = "expression") %>%
    mutate(gene = factor(gene, levels = markers_ss)) %>%
    # summarize data - mean gene expression in group and zone
    group_by(group, zone, gene) %>%
    summarise(
        meanExpression = mean(expression)) %>%
    as.data.frame() %>%
    # scale gene expression of each gene across all cells
    group_by(gene) %>%
    mutate(
        exp_scaled = 
            (meanExpression - mean(meanExpression)) / sd(meanExpression)) %>%
    as.data.frame() %>%
    arrange(gene, zone) %>%
    mutate(seg_name = paste(gene, zone)) %>%
    select(group, seg_name, exp_scaled) %>%
    pivot_wider(names_from = seg_name, values_from = exp_scaled) %>%
    column_to_rownames("group") %>%
    as.matrix()
colnames(summary_data) <- str_split_i(colnames(summary_data), " ", i = 2)

col_func <- colorRamp2(c(min(summary_data), 0, max(summary_data)), 
                       c("#2166ACFF", "#F7F7F7FF", "#B2182BFF")) 
gene_zone_hm <- Heatmap(summary_data, 
                        name = paste0("Scaled\nexpression"),
                        col = col_func,
                        na_col = "grey95",
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 12),
                        column_names_gp = gpar(fontsize = 12),
                        column_names_side = "bottom",
                        column_split = factor(rep(markers_ss, each = 3), 
                                              levels = markers_ss))



#### Fig. 6C ####

# excluding Sham cells from BZ and IZ as Sham shows no zonation
speX <- spe_fib[, !(spe_fib$condition == "Sham" & 
                        spe_fib$major_zone_ed %in% c("BZ", "IZ"))]

zn_ss_plot <- table(speX$sample_id, speX$stabmapHarmonyDataset_lt_exact,
                    speX$major_zone_ed) %>%
    as.data.frame() %>%
    # remove row combinations with no cells
    filter(rowSums(across(Freq) != 0) > 0) %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var2, group = Var2)) +
    geom_bar(stat = "identity") + 
    facet_wrap(~Var3, scales = "free_x") +
    scale_fill_manual(values = fib_cols) +
    labs(y = "Number of cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14, face = "bold"))




#### Extended Fig. 6A ####

gene_df <- data.frame(cells = colnames(spe_fib),
                      sample = spe_fib$sample_id,
                      x = reducedDim(spe_fib, "spatial_transform")[, 1],
                      y = reducedDim(spe_fib, "spatial_transform")[, 2])
# Vaibhao markers - project notes
gene_ss <- c("F-SH" = "Ly6a", "F-SL" = "Ndufa4l2", "F-Trans" = "Fgl2", 
             "F-WntX" = "Wif1", "F-Act" = "Postn", "IR" = "Cxcl5", 
             "F-CI" = "Top2a", "F-Cyc" = "Cdc20", "MYO" = "Col1a1", 
             "MFC" = "Comp", "F-IFNS" = "Cxcl10")

gene_plt_list <- setNames(lapply(gene_ss, function(gene) {
    gene_df <- gene_df %>%
        mutate(gene_exp = as.numeric(assay(spe_fib, "logcounts_area")[gene,]),
               gene_exp_z = scale(gene_exp)) %>%
        arrange(gene_exp_z)
    max_lim <- max(gene_df$gene_exp_z)
    min_lim <- min(gene_df$gene_exp_z)
    
    plt <- ggplot(gene_df, aes(x = x, y = y, color = gene_exp_z)) +
        geom_point(size = 0.05) +
        facet_wrap(~sample, nrow = 1) + 
        scale_color_viridis_c() +
        scale_color_gradientn(
            colors = c("#5a53a4", "#5ab0ab", "#fffbbb", "orange", "#a10741"),
            values = rescale(c(min_lim, min_lim/2, 0, max_lim/2, max_lim)),
            limits = c(min_lim, 5), oob = squish) +
        labs(color = gene) +
        theme_classic() +
        theme(plot.background = element_rect(fill = "black"),
              panel.background = element_rect(fill = "black"),
              panel.border = element_blank(),
              axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 20, colour = "white"),
              legend.background = element_rect(fill = "black"),
              legend.text = element_text(size = 14, colour = "white"),
              legend.title = element_text(size = 20, colour = "white"))
    
    return(plt)
}), names(gene_ss))




#### Extended Fig. 6C ####

# extract log normalised expression of IFNS marker gene
spe$IFNS_sig <- assay(spe, "logcounts_area")["Cxcl10", ]
# add IFNS detection status for F-IFNS substate, other fibroblasts, 
# and all other cells
spe$IFNS_detect <- case_when(
    spe$stabmapHarmonyDataset_lt_exact == "F-IFNS" & spe$IFNS_sig > 0 ~ "F-IFNS",
    spe$Major.groups == "Fibroblast" & 
        spe$stabmapHarmonyDataset_lt_exact != "F-IFNS" & 
        spe$IFNS_sig > 0 ~ "Other CF",
    spe$Major.groups != "Fibroblast" & spe$IFNS_sig > 0 ~ "Other cells",
    TRUE ~ "Background")
spe$IFNS_detect <- factor(spe$IFNS_detect, levels = c("Background", "Other CF", 
                                                      "F-IFNS", "Other cells"))

fifns_plt <- colData(spe) %>%
    as.data.frame() %>%
    arrange(IFNS_detect) %>%
    ggplot(aes(x = x_transform, y = y_transform, color = IFNS_detect, 
               size = IFNS_detect)) +
    geom_point() +
    scale_color_manual(
        name = NULL, 
        values = c("Background" = "gray97", "Other" = "deeppink2",
                   "Other CF" = "navy", "F-IFNS" = "green")) +
    
    scale_size_manual(name = NULL, 
                      values = c("Background" = 0.01, "Other" = 0.05, 
                                 "Other CF" = 0.05, "F-IFNS" = 0.05)) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    # facet by replicate along the rows and timepoint along the columns
    facet_grid(replicate ~ condition, switch = "y") +
    theme_void() +
    theme(strip.text = element_text(size = 12, face = "bold"),
          legend.position = "bottom")
