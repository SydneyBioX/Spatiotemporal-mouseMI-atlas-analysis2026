library(Seurat)
library(ranger)
library(e1071)
library(caret)
library(ROCR)
library(gbm)

### From Seurat object, define marker genes
load( "GFP-day3_panel_genes.RData" )

TSNEPlot(gfp.day3, pt.size = 0.5, label = TRUE, cols = col.set)

exp.markers <- as.matrix(GetAssayData(gfp.day3, assay = "RNA", slot = "counts"))
all.genes <- rownames(exp.markers)

#exp.means <- Seurat::AverageExpression(aggr, features = all.genes, use.counts = TRUE)
exp.means <- Seurat::AverageExpression(gfp.day3, slot = "counts")
exp.means <- exp.means$RNA

## Generate downsampled counts across all genes/clusters
div.by <- 10
exp.means.downsampled <- exp.means / div.by

set.seed(1)
mat.downsampled <- exp.markers
clusters <- names(table(Idents(gfp.day3)))
for (gene in all.genes) {
  for (cl in clusters) {
    sampled.counts <- rpois(table(Idents(gfp.day3))[cl], exp.means.downsampled[gene, cl])
    mat.downsampled[gene, which(Idents(gfp.day3) == cl)] <- sampled.counts
  }
}

hist(colSums(mat.downsampled),  main = paste0("Poisson, lambda = mean / ",div.by), xlab = "Total cell gene counts")

sum(colSums(mat.downsampled) == 0)
sum(colSums(mat.downsampled) < 5)
sum(colSums(mat.downsampled) < 2)
sum(colSums(mat.downsampled) > 60)
mean(colSums(mat.downsampled))
median(colSums(mat.downsampled))

#exp.markers <- exp.markers[1:20, 1:10]
summed.counts <- colSums(exp.markers)

## remove cells with 1 or fewer counts
class.labels <- as.character(Idents(gfp.day3))
cells.keep <- which(colSums(mat.downsampled) > 1)

class.labels <- class.labels[cells.keep]
mat.downsampled <- mat.downsampled[, cells.keep]

class.lineage.labels <- class.labels
## remove cycling cells

#non.cycling <- which(class.lineage.labels != "Cyc")
#class.lineage.labels <- class.lineage.labels[non.cycling]
#mat.downsampled <- mat.downsampled[, non.cycling]


############################################################################
### Train a Random Forest classifier to detect fibroblast subpopulations ###
############################################################################

set.seed(1)
folds <- createFolds(class.labels, k = 10, list = TRUE, returnTrain = FALSE)
#Iterate through the training/test folds

label.set <- class.labels
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
  
  rf.classifier <- ranger(x = train,
                          y = as.factor(train.labels),
                          num.trees = 1000,
                          mtry = 25)
  print("Fibroblast classifier trained")
  
  #pred.test <- predict(rf.classifier, data.frame(test))
  pred.fibroblasts <- predict(rf.classifier, test)
  
  pred.test <- as.character(pred.fibroblasts$predictions)
  
  test.labels <- label.set[thisFold]
  
  all.predictions <- append(all.predictions, as.character(pred.test))
  
  all.test.labels <- c(all.test.labels, test.labels)
  foldNum = foldNum + 1
  print(paste0(foldNum, " folds completed"))
}

cm <- confusionMatrix( data=as.factor(all.predictions), reference=as.factor(all.test.labels) )
cm
cm$byClass
class.accuracy.matrix <- cm$byClass
rownames(class.accuracy.matrix) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix))

## Plot the prediction accuracy derived from the confusion matrix
ggData <- data.frame(
  CellType = rep(rownames(class.accuracy.matrix), 4),
  Metric = c(rep("Precision", 7), rep("Recall", 7), 
             rep("Balanced Accuracy", 7), rep("Specificity", 7)),
  Value = c(class.accuracy.matrix[, "Precision"], class.accuracy.matrix[, "Recall"],
            class.accuracy.matrix[, "Balanced Accuracy"], class.accuracy.matrix[, "Specificity"])
)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
plot.label <- paste0("figures_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_GFP-day3_predictions_accuracy_bar_11x4inch.pdf")
ggsave(pl, filename = plot.label, width = 11, height = 4, units = "in", device = cairo_pdf)

file.label <- paste0("results_final/panel_all_",curated.list.label, "_RF_poissonDivide",div.by,"Sim_GFP-day3_predictions_accuracy.csv")
write.csv(ggData, file = file.label, row.names = FALSE)



