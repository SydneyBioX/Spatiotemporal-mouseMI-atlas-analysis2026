library(SingleCellExperiment)

set.seed(2024)


# load all the cells with zone information
spe <- readRDS("./cartana_spe_annotated.Rds")

# subset to fibroblasts
spe_fib <- spe[, spe$Major.groups == "Fibroblast"]

# extract CF cells metadata
fibro <- colData(spe_fib)

# reassign the labels to the full data
spe$patrick_label_transfer <- NA
colData(spe)[rownames(fibro),"patrick_label_transfer"] <- 
    fibro$patrick_label_transfer

spe$patrick_label_transfer <- factor(spe$patrick_label_transfer,
                                     levels = c("F-SH", "F-SL", "F-Trans",
                                                "F-WntX", "F-Act", "IR",
                                                "F-CI", "F-Cyc", "MYO", "MFC",
                                                "F-IFNS"))

# gene signatures of fibrotic programs
anti <- c("Ccn5", "Gucy1a1", "Sfrp2", "Fbln1")
pro <- c("Tgfb1", "Thbs4", "Crlf1", "Col15a1")

spe$fibrotic_score_anti <- colMeans(assay(spe, "logcounts_area")[anti,])
spe$fibrotic_score_pro <- colMeans(assay(spe, "logcounts_area")[pro,])

spe_sub <- spe[c(anti, pro), !is.na(spe$patrick_label_transfer)]

clusts <- kmeans(colData(spe_sub)[,c("fibrotic_score_anti", "fibrotic_score_pro")], 
                 centers = 2)

# identify anti-fibrotic cluster
anti_clust <- which.max(tapply(spe_sub$fibrotic_score_anti, clusts$cluster, mean))

spe_sub$fibrotic_score_class <- "Pro-fibrotic"
spe_sub$fibrotic_score_class[clusts$cluster == anti_clust] <- "Anti-fibrotic"

spe$fibrotic_score_class <- NA
colData(spe)[colnames(spe_sub), "fibrotic_score_class"] <- 
    spe_sub$fibrotic_score_class

saveRDS(spe, file = "./cartana_spe_annotated.Rds")
