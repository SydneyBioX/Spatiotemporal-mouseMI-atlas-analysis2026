library(Seurat)
library(ranger)
library(e1071)
library(caret)
library(ROCR)
library(gbm)
library(dplyr)
source("classifier_train.R")

### From Seurat object, define marker genes
seurat.file <- "TIP_54cl_panel_genes.RData"
load(seurat.file)

p.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO",
              "EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "Cyc",
              "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-WntX", "NKC", "Glial")

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
             "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
             "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")

TSNEPlot(tip.aggr, pt.size = 0.5, label = TRUE, cols = col.set)

curated.list.label <- "curated60_updated2"

exp.markers <- as.matrix(GetAssayData(tip.aggr, assay = "RNA", slot = "counts"))
all.genes <- rownames(exp.markers)
length(all.genes)

exp.means <- Seurat::AverageExpression(tip.aggr, slot = "counts")
exp.means <- exp.means$RNA

## Plot a few example genes
hist(exp.markers["Col1a1", which(Idents(tip.aggr) == "MYO")], main = "MYO: Col1a1 counts", xlab = "Cell gene counts")

div.by <- 10

exp.means.downsampled <- exp.means / div.by
sampled.counts <- rpois(table(Idents(tip.aggr))["MYO"], exp.means.downsampled["Col1a1", "MYO"])
hist(sampled.counts, main = paste0("MYO: Col1a1 Poisson sample, lambda = mean / ", div.by), xlab = "Cell gene counts")

## Generate downsampled counts across all genes/clusters
exp.means.downsampled <- exp.means / div.by

set.seed(1)
mat.downsampled <- exp.markers
clusters <- names(table(Idents(tip.aggr)))
for (gene in all.genes) {
  for (cl in clusters) {
    sampled.counts <- rpois(table(Idents(tip.aggr))[cl], exp.means.downsampled[gene, cl])
    mat.downsampled[gene, which(Idents(tip.aggr) == cl)] <- sampled.counts
  }
}

hist(colSums(mat.downsampled),  main = paste0("Poisson, lambda = mean / ", div.by), xlab = "Total cell gene counts")

sum(colSums(mat.downsampled) == 0)
sum(colSums(mat.downsampled) < 5)
sum(colSums(mat.downsampled) < 2)
sum(colSums(mat.downsampled) > 60)
mean(colSums(mat.downsampled))
median(colSums(mat.downsampled))

#exp.markers <- exp.markers[1:20, 1:10]
summed.counts <- colSums(exp.markers)

## remove cells with 1 or fewer counts
class.labels <- as.character(Idents(tip.aggr))
cells.keep <- which(colSums(mat.downsampled) > 1)

class.labels <- class.labels[cells.keep]
mat.downsampled <- mat.downsampled[, cells.keep]

### Define lineages
fibroblasts <- c("F-Act", "F-SH", "F-SL", "MYO", "F-WntX")
enothelial <- c("EC1", "EC2", "EC3", "Mural")
myeloid <- c("M1M\u03A6", "M2M\u03A6", "M1Mo", "DC", "MAC-IFNIC", "MAC8", "MAC6", "MAC7", "MAC-TR")
lymphoid <- c("TC1-Cd8", "TC2-Cd4", "BC", "NKC")
glial <- "Glial"

class.lineage.labels <- class.labels
## remove cycling cells

non.cycling <- which(class.lineage.labels != "Cyc")
class.lineage.labels <- class.lineage.labels[non.cycling]
mat.downsampled <- mat.downsampled[, non.cycling]


### Define lineages
fibroblasts <- c("F-Act", "F-SH", "F-SL", "MYO", "F-WntX")
endothelial <- c("EC1", "EC2", "EC3", "Mural")
myeloid <- c("M1M\u03A6", "M2M\u03A6", "M1Mo", "DC", "MAC-IFNIC", "MAC8", "MAC6", "MAC7", "MAC-TR")
lymphoid <- c("TC1-Cd8", "TC2-Cd4", "BC", "NKC")
glial <- "Glial"


#################################################################################
### Train a hierarchical classifer system to detect the individual cell types ###
#################################################################################

set.seed(1)
folds <- createFolds(class.lineage.labels, k = 10, list = TRUE, returnTrain = FALSE)
#Iterate through the training/test folds

label.set <- class.lineage.labels
foldNum = 0
feature.matrix <- t(as.matrix(mat.downsampled))
feature.matrix <- as.data.frame(feature.matrix)
all.predictions.continuous <- c()
all.predictions <- c()
all.test.labels <- c()
for (thisFold in folds) {
  test <- feature.matrix[thisFold, ]
  train <- feature.matrix[-thisFold, ]
  
  train.labels <- label.set[-thisFold]
  
  classifier.list <- TrainAllClassifiers(train, train.labels)
  
  #glial.rf <- classifier.list[[1]]
  lineage.rf <- classifier.list[[1]]
  fibroblast.rf <- classifier.list[[2]]
  ec.rf <- classifier.list[[3]]
  myeloid.rf <- classifier.list[[4]]
  lymphoid.rf <- classifier.list[[5]]
  
  #pred.test <- predict(rf.classifier, data.frame(test))
  #pred.glial <- predict(glial.rf, test)
  pred.lineage <- predict(lineage.rf, test)
  pred.fibroblasts <- predict(fibroblast.rf, test)
  pred.ec <- predict(ec.rf, test)
  pred.myeloid <- predict(myeloid.rf, test)
  pred.lymphoid <- predict(lymphoid.rf, test)
  
  #pred.test <- as.character(pred.glial$predictions)
  pred.test <- as.character(pred.lineage$predictions)
  pred.test[which(pred.test == "Fibroblast")] <- as.character(pred.fibroblasts$predictions)[which(pred.test == "Fibroblast")]
  pred.test[which(pred.test == "Endothelial")] <- as.character(pred.ec$predictions)[which(pred.test == "Endothelial")]
  pred.test[which(pred.test == "Myeloid")] <- as.character(pred.myeloid$predictions)[which(pred.test == "Myeloid")]
  pred.test[which(pred.test == "Lymphoid")] <- as.character(pred.lymphoid$predictions)[which(pred.test == "Lymphoid")]
  
  test.labels <- label.set[thisFold]
  
  all.predictions <- append(all.predictions, as.character(pred.test))
  
  #all.predictions.continuous <- rbind(all.predictions.continuous, pred.cont)
  all.test.labels <- c(all.test.labels, test.labels)
  foldNum = foldNum + 1
  print(paste0(foldNum, " folds completed"))
}

## Plot the prediction accuracy derived from the confusion matrix for fibroblasts
cm.fibro <- confusionMatrix(data = droplevels(as.factor(all.predictions[which(all.predictions %in% fibroblasts & all.test.labels %in% fibroblasts)])),
                            reference = as.factor(all.test.labels[which(all.predictions %in% fibroblasts & all.test.labels %in% fibroblasts)]))

class.accuracy.matrix.fibro <- cm.fibro$byClass
rownames(class.accuracy.matrix.fibro) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix.fibro))
ggData <- data.frame(
  CellType = rep(rownames(class.accuracy.matrix.fibro), 4),
  Metric = c(rep("Precision", 5), rep("Recall", 5), 
             rep("Balanced Accuracy", 5), rep("Specificity", 5)),
  Value = c(class.accuracy.matrix.fibro[, "Precision"], class.accuracy.matrix.fibro[, "Recall"],
            class.accuracy.matrix.fibro[, "Balanced Accuracy"], class.accuracy.matrix.fibro[, "Specificity"])
)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
#plot.label <- paste0("figures_final/panelLFC1_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_fibroblasts_predictions_accuracy_bar_9x4inch.pdf")
plot.label <- paste0("figures_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_fibroblasts_predictions_accuracy_bar_9x4inch.pdf")

ggsave(pl, filename = plot.label, width = 9, height = 4, units = "in", device = cairo_pdf)

file.label <- paste0("results_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_fibroblasts_predictions_accuracy.csv")
write.csv(ggData, file = file.label, row.names = FALSE)

## Plot the prediction accuracy derived from the confusion matrix for myeloid cells
cm.myl <- confusionMatrix(data = droplevels(as.factor(all.predictions[which(all.predictions %in% myeloid & all.test.labels %in% myeloid)])),
                          reference = as.factor(all.test.labels[which(all.predictions %in% myeloid & all.test.labels %in% myeloid)]))

class.accuracy.matrix.myl <- cm.myl$byClass
rownames(class.accuracy.matrix.myl) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix.myl))
ggData <- data.frame(
  CellType = rep(rownames(class.accuracy.matrix.myl), 4),
  Metric = c(rep("Precision", 9), rep("Recall", 9), 
             rep("Balanced Accuracy", 9), rep("Specificity", 9)),
  Value = c(class.accuracy.matrix.myl[, "Precision"], class.accuracy.matrix.myl[, "Recall"],
            class.accuracy.matrix.myl[, "Balanced Accuracy"], class.accuracy.matrix.myl[, "Specificity"])
)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
plot.label <- paste0("figures_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_myeloid_predictions_accuracy_bar_11x4inch.pdf")
ggsave(pl, filename = plot.label, width = 11, height = 4, units = "in", device = cairo_pdf)

file.label <- paste0("results_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_myeloid_predictions_accuracy.csv")
write.csv(ggData, file = file.label, row.names = FALSE)

## Plot the prediction accuracy derived from the confusion matrix for lymphoid cells
cm.lymphoid <- confusionMatrix(data = droplevels(as.factor(all.predictions[which(all.predictions %in% lymphoid & all.test.labels %in% lymphoid)])),
                               reference = as.factor(all.test.labels[which(all.predictions %in% lymphoid & all.test.labels %in% lymphoid)]))

class.accuracy.matrix.lymphoid <- cm.lymphoid$byClass
rownames(class.accuracy.matrix.lymphoid) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix.lymphoid))
ggData <- data.frame(
  CellType = rep(rownames(class.accuracy.matrix.lymphoid), 4),
  Metric = c(rep("Precision", 4), rep("Recall", 4), 
             rep("Balanced Accuracy", 4), rep("Specificity", 4)),
  Value = c(class.accuracy.matrix.lymphoid[, "Precision"], class.accuracy.matrix.lymphoid[, "Recall"],
            class.accuracy.matrix.lymphoid[, "Balanced Accuracy"], class.accuracy.matrix.lymphoid[, "Specificity"])
)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
plot.label <- paste0("figures_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_lymphoid_predictions_accuracy_bar_8.5x4inch.pdf")
ggsave(pl, filename = plot.label, width = 8.5, height = 4, units = "in", device = cairo_pdf)

file.label <- paste0("results_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_lymphoid_predictions_accuracy.csv")
write.csv(ggData, file = file.label, row.names = FALSE)

## Plot the prediction accuracy derived from the confusion matrix for endothelial cells
cm.ec <- confusionMatrix(data = droplevels(as.factor(all.predictions[which(all.predictions %in% endothelial & all.test.labels %in% endothelial)])),
                         reference = as.factor(all.test.labels[which(all.predictions %in% endothelial & all.test.labels %in% endothelial)]))

class.accuracy.matrix.ec <- cm.ec$byClass
rownames(class.accuracy.matrix.ec) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix.ec))
ggData <- data.frame(
  CellType = rep(rownames(class.accuracy.matrix.ec), 4),
  Metric = c(rep("Precision", 4), rep("Recall", 4), 
             rep("Balanced Accuracy", 4), rep("Specificity", 4)),
  Value = c(class.accuracy.matrix.ec[, "Precision"], class.accuracy.matrix.ec[, "Recall"],
            class.accuracy.matrix.ec[, "Balanced Accuracy"], class.accuracy.matrix.ec[, "Specificity"])
)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
plot.label <- paste0("figures_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_endothelial_predictions_accuracy_bar_8.5x4inch.pdf")
ggsave(pl, filename = plot.label, width = 8.5, height = 4, units = "in", device = cairo_pdf)

file.label <- paste0("results_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_endothelial_predictions_accuracy.csv")
write.csv(ggData, file = file.label, row.names = FALSE)

### Finally record the accuracy of lineage predictions
all.test.labels[which(all.test.labels %in% fibroblasts)] <- "Fibroblast"
all.test.labels[which(all.test.labels %in% enothelial)] <- "Endothelial"
all.test.labels[which(all.test.labels %in% myeloid)] <- "Myeloid"
all.test.labels[which(all.test.labels %in% lymphoid)] <- "Lymphoid"
all.test.labels[which(all.test.labels == glial)] <- "Glial"

all.predictions[which(all.predictions %in% fibroblasts)] <- "Fibroblast"
all.predictions[which(all.predictions %in% enothelial)] <- "Endothelial"
all.predictions[which(all.predictions %in% myeloid)] <- "Myeloid"
all.predictions[which(all.predictions %in% lymphoid)] <- "Lymphoid"
all.predictions[which(all.predictions == glial)] <- "Glial"

cm <- confusionMatrix( data=as.factor(all.predictions), reference=as.factor(all.test.labels) )
cm$byClass
class.accuracy.matrix <- cm$byClass
rownames(class.accuracy.matrix) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix))

ggData <- data.frame(
  CellType = rep(rownames(class.accuracy.matrix), 4),
  Metric = c(rep("Precision", 5), rep("Recall", 5), 
             rep("Balanced Accuracy", 5), rep("Specificity", 5)),
  Value = c(class.accuracy.matrix[, "Precision"], class.accuracy.matrix[, "Recall"],
            class.accuracy.matrix[, "Balanced Accuracy"], class.accuracy.matrix[, "Specificity"])
)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
plot.label <- paste0("figures_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_lineage_predictions_accuracy_bar_9x4inch.pdf")
ggsave(pl, filename = plot.label, width = 9, height = 4, units = "in", device = cairo_pdf)

file.label <- paste0("results_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_lineage_predictions_accuracy.csv")
write.csv(ggData, file = file.label, row.names = FALSE)

