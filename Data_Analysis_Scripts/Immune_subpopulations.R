library(SpatialExperiment)
library(dplyr)

set.seed(2026)


# load all the cells with zone, fibrotic program, and CF label information 
spe <- readRDS("./processed_data/cartana_spe_annotated_relabel.Rds")



#### Immune signatures ####

# marker genes of specific immune subpopulations
# neutrophils
Neutrophils <- c("S100a8", "Fut4")
# inflammatory macrophages
InflamMacro <- c("Cd14", "Ccr2", "Fcgr1")
# reparative macrophages
RepairMacro <- c("Cx3cr1", "C1qa", "C1qb", "Cd163", "Folr2", "Mgl2", "Csf1r", 
                 "Trem2")

# calculate marker gene signatures
spe$Neutrophil_score <- 
    colMeans(as.matrix(assay(spe, "logcounts_area")[Neutrophils, ]))
spe$InflamMacro_score <- 
    colMeans(as.matrix(assay(spe, "logcounts_area")[InflamMacro, ]))
spe$RepairMacro_score <- 
    colMeans(as.matrix(assay(spe, "logcounts_area")[RepairMacro, ]))



#### Cluster immune signatures ####

# kmeans clustering of individual immune subpopulation signatures to identify 
# immune vs background cells

neutro_km <- kmeans(spe$Neutrophil_score, centers = 2)
inflam_km <- kmeans(spe$InflamMacro_score, centers = 2)
repair_km <- kmeans(spe$RepairMacro_score, centers = 2)

# combine information
spe$Neutrophil_group <- 
    ifelse(neutro_km$cluster == which.max(neutro_km$centers) & 
               spe$Major.groups == "Immune", "Neutrophil", "Other")
spe$InflamMacro_group <- 
    ifelse(inflam_km$cluster == which.max(inflam_km$centers) & 
               spe$Major.groups == "Immune", "Inflammatory_macrophage", "Other")
spe$RepairMacro_group <- 
    ifelse(repair_km$cluster == which.max(repair_km$centers) & 
               spe$Major.groups == "Immune", "Reparative_macrophage", "Other")




#### Allocate immune labels ####

# Immune label = (Neutrophil) | (Inflammatory macrophage) | (Reparative macrophage)?
is_neutro <- spe$Neutrophil_group == "Neutrophil"
is_inflam <- spe$InflamMacro_group == "Inflammatory_macrophage"
is_repair <- spe$RepairMacro_group == "Reparative_macrophage"

# Scale normalised signature for comparability between subpopulations
# Allocate initial Immune label = (Neutrophil) | (Macrophage) based on higher z-score
spe$Neutrophil_score_scaled <- as.vector(scale(spe$Neutrophil_score))
spe$Macro_score_scaled_max <- as.vector(pmax(scale(spe$InflamMacro_score), 
                                             scale(spe$RepairMacro_score)))
spe$specificImmune_lineage <- case_when(
    is_neutro & !is_inflam & !is_repair ~ "Neutrophil",
    is_neutro & (is_inflam | is_repair) & 
        (spe$Neutrophil_score_scaled > spe$Macro_score_scaled_max) ~ "Neutrophil",
    
    !is_neutro & (is_inflam | is_repair) ~ "Macrophage",
    is_neutro & (is_inflam | is_repair) & 
        (spe$Neutrophil_score_scaled <= spe$Macro_score_scaled_max) ~ "Macrophage",
    
    TRUE ~ "Other")

# For initial Immune label = Macrophage, allocate (inflammatory) | (reparative) 
# based on Ccr2 expression as tie-breaker
# Ccr2 expression threshold = 0.01, based on its overall expression distribution 
spe$Ccr2_exp <- assay(spe, "logcounts_area")["Ccr2", ]

spe$specificImmune_group <- case_when(
    spe$specificImmune_lineage == "Neutrophil" ~ "Neutrophil",
    
    spe$specificImmune_lineage == "Macrophage" & is_inflam & !is_repair ~
        "Inflammatory_Macrophage",
    spe$specificImmune_lineage == "Macrophage" & !is_inflam & is_repair ~
        "Reparative_Macrophage",
    
    spe$specificImmune_lineage == "Macrophage" & is_inflam & is_repair & 
        (spe$Ccr2_exp >= 0.01) ~ "Inflammatory_Macrophage",
    spe$specificImmune_lineage == "Macrophage" & is_inflam & is_repair & 
        (spe$Ccr2_exp < 0.01) ~ "Reparative_Macrophage",
    
    TRUE ~ "Other")


# final immune subpopulation labels
spe$Immune_ct <- case_when(
    spe$Major.groups == "Immune" & spe$specificImmune_group != "Other" ~
        spe$specificImmune_group, 
    spe$Major.groups == "Immune" & spe$specificImmune_group == "Other" ~
        "Other_Immune",
    TRUE ~ "Other")
spe$Immune_ct <- factor(spe$Immune_ct, 
                        levels = c("Neutrophil", "Inflammatory_Macrophage",
                                   "Reparative_Macrophage", "Other_Immune", 
                                   "Other"))

saveRDS(spe, "./processed_data/spe_processed_final.Rds")

