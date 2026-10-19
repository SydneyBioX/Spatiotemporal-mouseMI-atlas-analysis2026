library(SpatialExperiment)
library(scater)

set.seed(2026)


# load processed data 
spe <- readRDS("./spe_processed_final.Rds")
spe$sample_id <- factor(spe$sample_id, 
                        levels = c("Sham_1", "Sham_2", "Day1_1", "Day1_2", 
                                   "Day3_1", "Day3_2", "Day5_1", "Day5_2", 
                                   "Day7_1", "Day7_2", "Day14_1", "Day14_2"))
# subset to fibroblasts
spe_fib <- spe[, spe$Major.groups == "Fibroblast"]

# account for cells with false F-WntX signal
spe_fib$Wif1_exp <- as.numeric(assay(spe_fib, "logcounts_area")["Wif1",])
# add Wif1 detection status
spe_fib$Wif1_detetct <- ifelse(spe_fib$Wif1_exp == 0, "No", "Yes")
# relabel F-WntX to F-WntX (Wif1), if the cell has Wif1 expression, else NA
spe_fib$fibStates <- 
    ifelse(spe_fib$stabmapHarmonyDataset_lt_exact == "F-WntX",
           ifelse(spe_fib$Wif1_detetct == "Yes", "F-WntX (Wif1)", NA),
           as.character(spe_fib$stabmapHarmonyDataset_lt_exact))
spe_fib$fibStates <- factor(spe_fib$fibStates, 
                            levels = c("F-SH", "F-SL", "F-Trans", "F-WntX (Wif1)", 
                                       "F-Act", "IR", "F-CI", "F-Cyc", "MYO",
                                       "MFC", "F-IFNS"))
# remove F-WntX cells with no Wif1 expression
spe_fib <- spe_fib[, !is.na(spe_fib$fibStates)]

fib_cols <- setNames(c("#4CD02D", "#E88321", "#DE331D", "#E9EA38", "#B7C1CD",
                       "#5CA6C9", "#F7CADE", "#479B74", "#ADF0E9", "#C850BA", 
                       "#4559B0"),
                     c("F-Trans", "F-SL", "F-Act", "F-SH", "MFC", "MYO", "F-IFNS",
                       "F-WntX", "F-CI", "IR", "F-Cyc"))



#### Panel A ####

indi_substate_plot <- plotReducedDim(spe_fib, dimred = "spatial", point_size = 3, 
                                     color_by = "fibStates", point_alpha = 1, 
                                     scattermore = TRUE) +
    facet_grid(spe_fib$fibStates ~ spe_fib$sample_id, switch = "y") + 
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = fib_cols) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 10, colour = "white"),
          strip.background = element_rect(fill = "black"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(size = 10, colour = "white"),
          legend.title = element_blank())



#### Panel B ####

all_substate_plot <- plotReducedDim(spe_fib, dimred = "spatial", point_size = 3, 
                                    color_by = "fibStates", point_alpha = 1, 
                                    scattermore = TRUE) +
    facet_wrap(~ spe_fib$sample_id, nrow = 2, scales = "free") + 
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = fib_cols) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 10, colour = "white"),
          strip.background = element_rect(fill = "black"),
          legend.background = element_rect(fill = "black"),
          legend.text = element_text(size = 10, colour = "white"),
          legend.title = element_blank())
