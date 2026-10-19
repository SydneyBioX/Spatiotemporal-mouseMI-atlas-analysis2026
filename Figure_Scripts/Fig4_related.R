library(SpatialExperiment)
library(ggplot2)
library(scattermore)
library(patchwork)
library(dplyr)
library(reshape2)
library(scater)

set.seed(2026)


# load Stabmap integration embedding with Harmony batch correction
# (Harmony_umap) and reference substate labels (ref_substate_labels)
load("./stabmap_integrate.RData")

# load processed data 
spe <- readRDS("./spe_processed_final.Rds")
spe$sample_id <- factor(spe$sample_id, 
                        levels = c("Sham_1", "Sham_2", "Day1_1", "Day1_2", 
                                   "Day3_1", "Day3_2", "Day5_1", "Day5_2", 
                                   "Day7_1", "Day7_2", "Day14_1", "Day14_2"))
# subset to fibroblasts
spe_fib <- spe[, spe$Major.groups == "Fibroblast"]

spe_fib <- runPCA(spe_fib, exprs_values = "logcounts_area",
                  name = "PCA_log_area")
spe_fib <- runUMAP(spe_fib, name = "UMAP_log_area")

# exclude F-WntX cells with false WntX signal (only keep Wif1+ F-WntX cells) 
# add Wif1 expression per cell
spe_fib$Wif1_exp <- as.numeric(assay(spe_fib, "logcounts_area")["Wif1",])
# add Wif1 detection status
spe_fib$Wif1_detetct <- ifelse(spe_fib$Wif1_exp == 0, "No", "Yes")
# relabel F-WntX to F-WntX (Wif1), if the cell has Wif1 expression, else NA
spe_fib$fibStates <- ifelse(spe_fib$stabmapHarmonyDataset_lt_exact == "F-WntX",
                            ifelse(spe_fib$Wif1_detetct == "Yes", "F-WntX",
                                   NA),
                            as.character(spe_fib$stabmapHarmonyDataset_lt_exact))
spe_fib$fibStates <- factor(spe_fib$fibStates, 
                            levels = c("F-SH", "F-SL", "F-Trans", "F-WntX", 
                                       "F-Act", "IR", "F-CI", "F-Cyc", "MYO",
                                       "MFC", "F-IFNS"))

fib_cols <- setNames(c("#4CD02D", "#E88321", "#DE331D", "#E9EA38", "#B7C1CD",
                       "#5CA6C9", "#F7CADE", "#479B74", "#ADF0E9", "#C850BA", 
                       "#4559B0"),
                     c("F-Trans", "F-SL", "F-Act", "F-SH", "MFC", "MYO", "F-IFNS",
                       "F-WntX", "F-CI", "IR", "F-Cyc"))



#### Fig. 4A ####

# data frame of integrated datasets
embed_df <- as.data.frame(harmony_umap) %>%
    mutate(substate = append(ref_substate_labels, 
                             as.character(spe_fib$stabmapHarmonyDataset_lt_exact)),
           dataset = append(rep("Patrick et al", length(ref_substate_labels)), 
                       rep("ISS data", ncol(spe_fib))))
    
plt_by_dataset <- embed_df %>%
    arrange(dataset) %>%
    ggplot(aes(x = UMAP1, y = UMAP2, colour = dataset)) +
    geom_point(size = 0.1, alpha = 1) +
    scale_color_manual(values = c("ISS data" = "#FA8072", 
                                  "Patrick et al" = "#ABD7EB")) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_classic() + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_text(),
          legend.text = element_text(size = 20, colour = "white"),
          legend.title = element_blank(),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          legend.background = element_rect(fill = "black"))

plt_by_substate <- ggplot(embed_df, aes(x = UMAP1, y = UMAP2, 
                                        colour = substate)) +
    geom_point(size = 0.1, alpha = 1) +
    scale_color_manual(values = fib_cols) + 
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_text(),
          legend.text = element_text(size = 20, colour = "white"),
          legend.title = element_blank(),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          legend.background = element_rect(fill = "black"))

plt_4a <- plt_by_dataset + plt_by_substate




#### Fig. 4B ####

substate_markers <- c("Postn", "Eln", "Col8a1", "Col1a1", "Cilp", "Timp1", 
                      "Vim", "Acta2", "Fn1", "Stmn1", "Cdc20", "Mki67", "Top2a", 
                      "Birc5", "Cxcl10", "Cxcl9", "Ddah1", "Cd40", "Tyrobp", 
                      "Fabp3", "Myh6", "Pi16", "Ly6a", "Tnnc1", "Gja1", "Nppa", 
                      "Fgl2", "Cfh", "C3", "Comp", "Fmod", "Ecrg4", "Fbln1", 
                      "Wif1", "Mt2", "Mt1", "Arg1", "Cxcl5", "Sfrp2", "Ccn5", 
                      "C1qa")

# extract normalised expression of marker genes
gexp_marker <- assay(spe_fib, "logcounts_area")[substate_markers, ] |> 
    as.matrix() |> t() |> as.data.frame()
# add substate label
gexp_marker$cluster <- factor(as.character(spe_fib$stabmapHarmonyDataset_lt_exact))
gexp_marker_melt <- reshape2::melt(gexp_marker, id.vars = "cluster", 
                                   variable.name = "gene", 
                                   value.name = "expression")

summary_data <- gexp_marker_melt %>%
    group_by(cluster, gene) %>%
    # mean gene expression in cluster and percentage of cluster cells 
    # in which gene is expressed
    summarise(
        meanExpression = mean(expression),
        PerctageCells = sum(expression > 0) / n() * 100) %>%
    as.data.frame() %>%
    # scale expression of each gene
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
          legend.position = "bottom")




#### Fig. 4C ####

# plot all samples together
umap_substates <- plotReducedDim(spe_fib, dimred = "UMAP_log_area", 
                                 color_by = "stabmapHarmonyDataset_lt_exact",
                                 point_alpha = 1, point_size = 0.1) +
    scale_color_manual(values = fib_cols) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.text = element_text(size = 30, colour = "white"),
          legend.title = element_blank(),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          legend.background = element_rect(fill = "black"))

# plot individual timepoints
umap_substates_samples <- 
    plotReducedDim(spe_fib, dimred = "UMAP_log_area", 
                   color_by = "stabmapHarmonyDataset_lt_exact",
                   point_alpha = 1, point_size = 0.1) +
    facet_wrap(~spe_fib$condition, ncol = 3) +
    scale_color_manual(values = fib_cols) +
    guides(color = guide_legend(override.aes = list(size = 15))) +
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          legend.text = element_text(size = 30, colour = "white"),
          legend.title = element_blank(),
          plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(size = 30, colour = "white"),
          legend.background = element_rect(fill = "black"))




#### Fig. 4D ####

substate_bar_plot <- table(spe_fib$sample_id, 
                           spe_fib$stabmapHarmonyDataset_lt_exact) %>%
    as.data.frame() %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = fib_cols) +
    labs(y = "Number of cells") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_blank())




#### Extended Fig. 4A ####

fib_table <- table(spe_fib$sample_id, spe_fib$fibStates)
fib_table <- rbind(fib_table, "Total cells" = colSums(fib_table))




#### Extended Fig. 4B ####

gene_df <- data.frame(cells = colnames(spe_fib),
                      sample = spe_fib$sample_id,
                      umap1 = reducedDim(spe_fib, "UMAP_log_area")[, 1],
                      umap2 = reducedDim(spe_fib, "UMAP_log_area")[, 2])

gene_ss <- c("F-SH" = "Ly6a", "F-SL" = "Ndufa4l2", "F-Trans" = "Fgl2", 
             "F-WntX" = "Wif1", "F-Act" = "Postn", "IR" = "Cxcl5", 
             "F-CI" = "Top2a", "F-Cyc" = "Cdc20", "MYO" = "Col1a1", 
             "MFC" = "Comp", "F-IFNS" = "Cxcl10")

gene_umap_plt_list <- setNames(lapply(gene_ss, function(gene) {
    gene_df <- gene_df %>%
        mutate(gene_exp = as.numeric(assay(spe_fib, "logcounts_area")[gene,]),
               gene_exp_z = scale(gene_exp)) %>%
        arrange(gene_exp_z)
    max_lim <- max(gene_df$gene_exp_z)
    min_lim <- min(gene_df$gene_exp_z)
    
    plt <- ggplot(gene_df, aes(x = umap1, y = umap2, color = gene_exp_z)) +
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

