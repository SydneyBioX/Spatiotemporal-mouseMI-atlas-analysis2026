library(MoleculeExperiment)
library(ggplot2)
library(SpatialExperiment)
library(dplyr)
library(scran)
library(reshape2)
library(scattermore)

set.seed(2026)


# load processed data 
spe <- readRDS("./spe_processed_final.Rds")

cardiacCells_cols <- setNames(c("#ff5732", "#18a085", "#1e618d", "#9b59b6", 
                                "#f29c13"), 
                              c("Cardiomyocyte", "Endothelial", "Fibroblast",
                                "Immune", "Stressed cardiomyocyte"))



#### Fig. 1B ####

me <- readRDS("./processed_data/ME_objects/me_Day7_2.Rds")

# plot nuclei
nuclei_plt <- ggplot_me() +
    geom_polygon_me(me, assayName = "nucleus", fill = "white") +
    coord_cartesian(xlim = c(10500, 11500), ylim = c(-13500, -12500)) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank())

# plot nuclei with cell boundaries
cell_plt <- ggplot_me() +
    coord_cartesian(xlim = c(10500, 11500), ylim = c(-13500, -12500)) +
    geom_polygon_me(me, assayName = "nucleus", fill = "blue") +
    geom_polygon_me(me, assayName = "cell", color = "yellow", fill = NA, 
                    linewidth = 1.5) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank())

# plot nuclei and cell boundaries with transcripts 
transcript_plt <- ggplot_me() +
    coord_cartesian(xlim = c(10500, 11500), ylim = c(-13500, -12500)) +
    geom_polygon_me(me, assayName = "nucleus", fill = "blue") +
    geom_point_me(me, size = 2, color = "magenta") +
    geom_polygon_me(me, assayName = "cell", color = "yellow", fill = NA, 
                    linewidth = 1) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank())



#### Fig. 1C ####
plot_1c <- table(spe$condition, spe$replicate) %>%
    reshape2::melt() %>%
    mutate(Var2 = as.character(Var2)) %>%
    ggplot(aes(x = Var1, y = value, group = Var2, fill = Var2, color = Var2)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("1" = "white", "2" = "grey50")) +
    scale_color_manual(values = c("black", "black")) +
    labs(x = NULL, y  = "Cell number") + 
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_blank())
    


#### Fig. 1D ####

# cardiac cell type markers
markers_ct <- c("Gja1", "Tnnc1", "Vegfa", "Cdh2", "Mt1", "Egfl7", "Cdh5", "Ly6a",
               "Kdr", "Cd34", "Col1a1", "Postn", "Fn1", "Eln", "Col15a1", "Apoe",
               "Cd74", "C1qa", "C1qb", "Tyrobp", "Nppa", "Clu", "Thbs4", "Cxcr6")

# extract normalised expression
gexp_marker <- assay(spe, "logcounts_area")[markers_ct, ] |> 
    as.matrix() |> t() |> as.data.frame()

# add cluster labels to expression data frame
gexp_marker$cluster <- factor(as.character(spe$Major.groups))
gexp_marker_melt <- reshape2::melt(gexp_marker, id.vars = "cluster", 
                                   variable.name = "gene", 
                                   value.name = "expression")

# summarize data
summary_data <- gexp_marker_melt %>%
    group_by(cluster, gene) %>%
    summarise(meanExpression = mean(expression),
              PerctageCells = sum(expression > 0) / n() * 100) %>%
    as.data.frame()

# scale gene expression
summary_data <- summary_data %>%
    group_by(gene) %>%
    mutate(exp_scaled = 
               (meanExpression - mean(meanExpression)) / sd(meanExpression))

markerPlot <- summary_data %>%
    ggplot(aes(x = gene, y = cluster, size = PerctageCells, 
               color = exp_scaled)) +
    geom_point() +
    scale_size_continuous(limit = c(0, 100), range = c(0, 5)) + 
    scale_color_gradient2(low = "gold", mid = "ivory", high = "red") +  
    labs(color = "Scaled\nexpression", size = "Cells\ndetected") +
    scale_y_discrete(limits = rev) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.box = "horizontal") 



#### Fig. 1E ####

# aggregate
umap_celltypes <- plotReducedDim(spe, dimred = "UMAP_new", 
                                 color_by = "Major.groups", point_alpha = 1,
                                 point_size = 0.1) +
    scale_color_manual(values = cardiacCells_cols) +
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
                          color_by = "Major.groups", point_alpha = 1,
                          point_size = 0.5) +
        scale_color_manual(values = cardiacCells_cols) +
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




#### Extended Fig. 1C ####

libsize_box <- colData(spe) %>%
    as.data.frame() %>%
    mutate(Major.groups = 
               factor(Major.groups, 
                      levels = c("Stressed cardiomyocyte", "Fibroblast", 
                                 "Cardiomyocyte", "Immune", "Endothelial"))) %>%
    ggplot(aes(x = Major.groups, y = total, fill = Major.groups)) +
    geom_boxplot() +
    scale_fill_manual(values = cardiacCells_cols) +
    labs(x = NULL, y = "Total transcripts\ndetected / cell") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))




#### Extended Fig. 1D ####

detect_box <- colData(spe) %>%
    as.data.frame() %>%
    mutate(Major.groups = 
               factor(Major.groups, 
                      levels = c("Fibroblast", "Stressed cardiomyocyte", 
                                 "Immune", "Endothelial", "Cardiomyocyte"))) %>%
    ggplot(aes(x = Major.groups, y = detected, fill = Major.groups)) +
    geom_boxplot() +
    scale_fill_manual(values = cardiacCells_cols) +
    labs(x = NULL, y = "Number of unique\ngenes detected / cell") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))




#### Extended Fig. 1E ####

me_files_list <- paste0("me_", c("Sham_1", "Sham_2", "Day1_1", "Day1_2", 
                               "Day3_1", "Day3_2", "Day5_1", "Day5_2", 
                               "Day7_1", "Day7_2", "Day14_1", "Day14_2"), ".Rds")

# plot transcripts
plots_all <- setNames(lapply(me_files_list, function(samp) {
    
    # load me object of sample
    me <- readRDS(me_list[[samp]])
    mols_df <- molecules(me, assayName = "detected", flatten = T)
    
    plt <- ggplot(mols_df, aes(x = x_location, y = y_location)) +
        geom_scattermore(pointsize = 1, alpha = 1, color = "grey20") +
        coord_fixed() +
        theme_void()
    
    return(plt)
}), me_files_list)
