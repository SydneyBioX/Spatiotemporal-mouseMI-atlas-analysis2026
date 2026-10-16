library(SingleCellExperiment)
library(scater)
library(scran)
library(BiocNeighbors)
library(StabMap)

set.seed(2024)

# load all the cells with broad cell type information
spe <- readRDS("./SCE_PatrickAno.rds")

kval <- 200 # number of neighbours to consider

# define zones based on cell type composition of 200 nearest neighbors
propwhich_list <- sapply(levels(spe$sample_id), function(i) {
    
    spe_sub <- spe[,spe$sample_id == i]
    spe_sub_sub <- spe_sub[,spe_sub$Major.groups %in% c("Cardiomyocyte", 
                                                        "Stressed cardiomyocyte", 
                                                        "Fibroblast")]
    
    ptab <- c(prop.table(table(spe_sub_sub$Major.groups)))
    exptab <- ptab*kval
    
    mat <- StabMap:::queryNamedKNN(reducedDim(spe_sub_sub, "spatial"), 
                                   reducedDim(spe_sub, "spatial"), k = kval)
    mat_names <- StabMap:::vectorSubset(setNames(spe_sub_sub$Major.groups, 
                                                colnames(spe_sub_sub)), mat)
    prop <- cbind(
        CM = rowMeans(1*(mat_names == "Cardiomyocyte")),
        Fibroblast = rowMeans(1*(mat_names == "Fibroblast")),
        SCM = rowMeans(1*(mat_names == "Stressed cardiomyocyte"))
    )
    nobs <- prop*kval
    
    resids <- t(t(t(t(nobs) - exptab)) / sqrt(exptab))
    
    propwhich <- c("RZ", "IZ", "BZ")[apply(resids, 1, which.max)]
    
    return(setNames(propwhich, colnames(spe_sub)))
}, simplify = FALSE)
names(propwhich_list) <- NULL
propwhich_named <- unlist(propwhich_list)

spe$major_zone <- as.character(propwhich_named[colnames(spe)])


# calculate distances to BZ
first_ind_list_BZ <- sapply(levels(spe$sample_id), function(i) {
    
    kval <- 5 # number of nearest neighbors
    
    spe_sub <- spe[,spe$sample_id == i]
    spe_sub_sub <- spe_sub[,spe_sub$major_zone %in% c("BZ")]
    
    out <- BiocNeighbors::queryKNN(
        reducedDim(spe_sub_sub, "spatial"), reducedDim(spe_sub, "spatial"), 
        k = kval, get.distance = TRUE,
        BNPARAM=ExhaustiveParam()
    )$distance
    
    dist_select <- out[,kval]
    
    return(setNames(dist_select, colnames(spe_sub)))
}, simplify = FALSE)
names(first_ind_list_BZ) <- NULL
spe$distance_to_BZ <- unlist(first_ind_list_BZ)[colnames(spe)]


# calculate distances to RZ
first_ind_list_RZ <- sapply(levels(spe$sample_id), function(i) {
    
    kval <- 5 # number of nearest neighbors
    
    spe_sub <- spe[,spe$sample_id == i]
    spe_sub_sub <- spe_sub[,spe_sub$major_zone %in% c("RZ")]
    
    out <- BiocNeighbors::queryKNN(
        reducedDim(spe_sub_sub, "spatial"), reducedDim(spe_sub, "spatial"), 
        k = kval, get.distance = TRUE,
        BNPARAM=ExhaustiveParam()
    )$distance
    
    dist_select <- out[,kval]
    
    return(setNames(dist_select, colnames(spe_sub)))
}, simplify = FALSE)
names(first_ind_list_RZ) <- NULL
spe$distance_to_RZ <- unlist(first_ind_list_RZ)[colnames(spe)]


# calculate distances to IZ
first_ind_list_IZ <- sapply(levels(spe$sample_id), function(i) {
    
    kval <- 5 # number of nearest neighbors
    
    spe_sub <- spe[,spe$sample_id == i]
    spe_sub_sub <- spe_sub[,spe_sub$major_zone %in% c("IZ")]
    
    out <- BiocNeighbors::queryKNN(
        reducedDim(spe_sub_sub, "spatial"), reducedDim(spe_sub, "spatial"), k = kval, get.distance = TRUE,
        BNPARAM=ExhaustiveParam()
    )$distance
    
    dist_select <- out[,kval]
    
    return(setNames(dist_select, colnames(spe_sub)))
}, simplify = FALSE)
names(first_ind_list_IZ) <- NULL
spe$distance_to_IZ <- unlist(first_ind_list_IZ)[colnames(spe)]

saveRDS(spe, file = "./cartana_spe_annotated.Rds")
