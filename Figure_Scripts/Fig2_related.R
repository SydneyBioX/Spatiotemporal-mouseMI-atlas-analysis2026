library(SpatialExperiment)
library(scuttle)
library(scater)
library(ggrepel)
library(dplyr)

set.seed(2026)

# load processed data 
spe <- readRDS("./spe_processed_final.Rds")

day_cols <- setNames(c("#e3191c", "#377db8", "#4eaf4a", "#974ea3", 
                       "#ff7f00", "#ff01ff"),
                     c("Sham", "Day1", "Day3", "Day5", "Day7", "Day14"))

cardiacCells_cols <- setNames(c("#ff5732", "#18a085", "#1e618d", "#9b59b6", 
                                "#f29c13"), 
                              c("Cardiomyocyte", "Endothelial", "Fibroblast",
                                "Immune", "Stressed cardiomyocyte"))

#### Panel A ####

# sum cell expression in each sample
spe_pbulk <- aggregateAcrossCells(spe, ids = spe$sample_id)

# remove genes with zero expression
spe_pbulk <- spe_pbulk[rowSums(counts(spe_pbulk) != 0) > 0, ]

rm_col <- "sizeFactor"
colData(spe_pbulk) <- 
    colData(spe_pbulk)[, colnames(colData(spe_pbulk)) != rm_col]

# normalise expression
spe_pbulk <- logNormCounts(spe_pbulk)
# run PCA
spe_pbulk <- runPCA(spe_pbulk)

# aggregate
pseudobulk_plt <- data.frame(pc1 = reducedDim(spe_pbulk, "PCA")[, 1],
                             pc2 = reducedDim(spe_pbulk, "PCA")[, 2],
                             condition = spe_pbulk$condition,
                             sample = spe_pbulk$sample_id) %>%
    ggplot(aes(x = pc1, y = pc2, label = sample)) +
    geom_point(aes(color = condition), size = 3) +
    geom_label_repel(aes(label = sample, fill = condition), color = "ivory",
                     segment.color = NA, label.size = 0, size = 4,
                     label.padding = unit(0.4, "lines"), box.padding = 0.5, 
                     max.overlaps = Inf, point.padding = 0, fontface = "bold") +
    scale_colour_manual(values = day_cols) +
    scale_fill_manual(values = day_cols) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = "none")


# per cell type
pseudobulk_plt_ct <- setNames(lapply(unique(spe$Major.groups), function(ct) {
    # subset data
    speX <- spe[, spe$Major.groups == ct]
    # sum cell expression in a sample
    speX_pbulk <- aggregateAcrossCells(speX, ids = speX$sample_id)
    print(dim(speX_pbulk))
    # remove genes with zero expression across samples
    speX_pbulk <- speX_pbulk[rowSums(counts(speX_pbulk) != 0) > 0, ]
    print(dim(speX_pbulk))
    # remove sizeFactor
    rm_col <- "sizeFactor"
    colData(speX_pbulk) <- 
        colData(speX_pbulk)[, colnames(colData(speX_pbulk)) != rm_col]
    
    speX_pbulk <- logNormCounts(speX_pbulk)
    speX_pbulk <- runPCA(speX_pbulk)
    plt <- data.frame(pc1 = reducedDim(speX_pbulk, "PCA")[, 1],
                      pc2 = reducedDim(speX_pbulk, "PCA")[, 2],
                      condition = speX_pbulk$condition,
                      sample = speX_pbulk$sample_id) %>%
        ggplot(aes(x = pc1, y = pc2, label = sample)) +
        geom_point(aes(color = condition), size = 3) +
        geom_label_repel(aes(label = sample, fill = condition), color = "ivory",
                         segment.color = NA, label.size = 0, size = 4,
                         label.padding = unit(0.4, "lines"), box.padding = 0.5, 
                         max.overlaps = Inf, point.padding = 0, fontface = "bold") +
        scale_colour_manual(values = day_cols) +
        scale_fill_manual(values = day_cols) +
        theme_classic() +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              legend.position = "none")
    return(plt)
}), unique(spe$Major.groups))



#### Panel B ####

indiCells_sp_plot <- plotReducedDim(spe, dimred = "spatial", 
                                    color_by = "Major.groups", point_size = 2, 
                                    point_alpha = 1, scattermore = TRUE) +
    facet_grid(spe$Major.groups ~ spe$sample_id, switch = "y") + 
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = cardiacCells_cols) +
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
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(size = 10, colour = "white"),
          legend.title = element_blank())



#### Panel C ####

allCells_sp_plot <- plotReducedDim(spe, dimred = "spatial", point_size = 0.5,
                                   color_by = "Major.groups",scattermore = TRUE,
                                   point_alpha = 1) +
    facet_wrap(~ spe$sample_id, scales = "free") +
    scale_color_manual(values = cardiacCells_cols) +
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
