library(MoleculeExperiment)
library(ggplot2)
library(SpatialExperiment)
library(dplyr)
library(scran)
library(reshape2)

set.seed(2026)

me <- readRDS("./processed_data/ME_objects/me_Day7_2.Rds")
spe <- readRDS("./cartana_spe_annotated_relabel.Rds")

cardiacCells_cols <- setNames(c("#ff5732", "#18a085", "#1e618d", "#9b59b6", 
                                "#f29c13"), 
                              c("Cardiomyocyte", "Endothelial", "Fibroblast",
                                "Immune", "Stressed cardiomyocyte"))



#### Panel B ####

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



#### Panel D ####

bar_plot <- table(spe$sample_id, 
                  spe$Major.groups) %>%
    as.data.frame() %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = cardiacCells_cols) +
    labs(y = "Number of cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_blank())



#### Panel E ####

# cardiac cell type markers
markers_v <- c("Gja1", "Tnnc1", "Vegfa", "Cdh2", "Mt1", "Egfl7", "Cdh5", "Ly6a",
               "Kdr", "Cd34", "Col1a1", "Postn", "Fn1", "Eln", "Col15a1", "Apoe",
               "Cd74", "C1qa", "C1qb", "Tyrobp", "Nppa", "Clu", "Thbs4", "Cxcr6")
# check markers are in dataset
top_markers <- intersect(markers_v, rownames(spe))


# extract normalised expression
gexp_marker <- assay(spe, "logcounts_area")[top_markers, ] |> 
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



#### Panel F ####

spe <- runUMAP(spe, dimred = "PCA_area_cosine_sum", n_threads = 16,
               n_neighbors = 10, min_dist = 0.5, n_dimred = 30, 
               metric = "cosine", name = "UMAP_new")

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
