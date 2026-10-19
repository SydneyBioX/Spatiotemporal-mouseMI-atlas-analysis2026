library(SpatialExperiment)
library(scater)
library(dplyr)
library(reshape2)
library(tidyr)
library(tibble)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)

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

zone_cols <- setNames(c("#e63946", "#f39d12", "#8b49e2"), c("RZ", "BZ", "IZ"))

cardiacCells_cols <- setNames(c("#ff5732", "#18a085", "#1e618d", "#9b59b6", 
                                "#f29c13"), 
                              c("Cardiomyocyte" , "Endothelial", "Fibroblast",
                                "Immune", "Stressed cardiomyocyte"))

#### Panel A ####

allCells_sp_plot <- plotReducedDim(spe, dimred = "spatial", point_size = 0.5,
                                   color_by = "major_zone_ed",scattermore = TRUE,
                                   point_alpha = 1) +
    facet_wrap(~ spe$sample_id, scales = "free") +
    scale_color_manual(values = zone_cols) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 10, colour = "white"),
          legend.position = "none")



#### Panel B ####

spe <- runUMAP(spe, dimred = "PCA_area_cosine_sum", n_threads = 16,
               n_neighbors = 10, min_dist = 0.5, n_dimred = 30, 
               metric = "cosine", name = "UMAP_new")

# aggregate
umap_celltypes <- plotReducedDim(spe, dimred = "UMAP_new", 
                                 color_by = "major_zone_ed", point_alpha = 1,
                                 point_size = 0.1) +
    scale_color_manual(values = zone_cols) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = "none",
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"))

# per timepoint
umap_celltypes_list <- setNames(lapply(unique(spe$condition), function(cond) {
    speX <- spe[, spe$condition == cond]
    plt <- plotReducedDim(speX, dimred = "UMAP_new", 
                          color_by = "major_zone_ed", point_alpha = 1,
                          point_size = 0.5) +
        scale_color_manual(values = zone_cols) +
        theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              legend.position = "none",
              plot.background = element_rect(fill = "black"),
              panel.background = element_rect(fill = "black"),
              panel.border = element_blank())
    return(plt)
}), unique(spe$condition))



#### Panel C ####

markers_ct <- c("Gja1", "Tnnc1", "Egfl7", "Cdh5", "Col1a1", "Postn", "Cd74",
                "Tyrobp", "Nppa", "Clu")

# extract log normalised expression
gexp_marker <- assay(spe, "logcounts_area")[markers_ct, ] |> 
    as.matrix() |> t() |> as.data.frame()
# add broad cell type labels
gexp_marker$group <- spe$Major.groups
# add zone labels
gexp_marker$zone <- spe$major_zone_ed

summary_data <- reshape2::melt(gexp_marker, id.vars = c("group", "zone"), 
                               variable.name = "gene",
                               value.name = "expression") %>%
    mutate(zone = factor(zone, levels = c("RZ", "BZ", "IZ")),
           gene = factor(gene, levels = markers_ct)) %>%
    # summarize data - mean gene expression in group and zone
    group_by(group, zone, gene, .drop = T) %>%
    summarise(meanExpression = mean(expression)) %>%
    as.data.frame() %>%
    # scale expression of each gene
    group_by(gene) %>%
    mutate(exp_scaled = 
               (meanExpression - mean(meanExpression)) / sd(meanExpression)) %>%
    as.data.frame() %>%
    arrange(gene, zone) %>%
    mutate(seg_name = paste(gene, zone)) %>%
    select(group, seg_name, exp_scaled) %>%
    pivot_wider(names_from = seg_name, values_from = exp_scaled)
summary_mat <- summary_data %>%
    column_to_rownames("group") %>%
    as.matrix()
colnames(summary_mat) <- str_split_i(colnames(summary_mat), " ", i = 2)
col_func <- colorRamp2(c(min(summary_mat), 0, max(summary_mat)), 
                       c("#2166ACFF", "#F7F7F7FF", "#B2182BFF")) 
gene_zone_hm <- Heatmap(summary_mat, 
                        name = paste0("Scaled\nExpression"),
                        col = col_func,
                        na_col = "grey95",
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 12),
                        column_names_gp = gpar(fontsize = 12),
                        column_names_side = "bottom",
                        column_split = factor(rep(markers_ct, each = 3), 
                                              levels = markers_ct))




#### Panel D ####

bar_plot_all <- table(spe$sample_id, spe$major_zone_ed) %>%
    as.data.frame() %>%
    # excluding Sham, which shows no zonation
    filter(!(Var1 %in% c("Sham_1", "Sham_2"))) %>%
    mutate(Var1 = factor(Var1, levels = levels(spe$sample_id))) %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = zone_cols) +
    labs(y = "Number of cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_blank())



#### Panel E ####

# excluding Sham from BZ and IZ as Sham shows no zonation
speX <- spe[, !(spe$condition == "Sham" & spe$major_zone_ed %in% c("BZ", "IZ"))]

zn_ct_plot <- table(speX$sample_id, speX$Major.groups, speX$major_zone_ed) %>%
    as.data.frame() %>%
    # removing row combinations with no cells
    filter(rowSums(across(Freq) != 0) > 0) %>%
    mutate(Var1 = factor(Var1, levels = levels(spe$sample_id))) %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var2, group = Var2)) +
    geom_bar(stat = "identity") + 
    facet_wrap(~Var3, scales = "free_x") +
    scale_fill_manual(values = cardiacCells_cols) +
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
