#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(rTASSEL)

# input parameter of number of traits in the phenotype data
col_num <- as.integer(args[1])

options(java.parameters = c("-Xmx2g", "-Xms1g"))
rTASSEL::startLogger(fullPath = getwd(), fileName = "RTassel.log")

tasGenoPheno <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = "genotype.hmp.txt",
  phenoPathDFOrObj = "phenotype.txt"
)

feature_names <- names(rTASSEL::getPhenotypeDF(tasObj = tasGenoPheno))
# Skip first column and last 5 PC1-5 columns to get a vector of trait names
feature_names <- feature_names[2:(length(feature_names)-5)]
tasKin <- readTasselDistanceMatrix("kinship.txt")

formula <- as.formula(paste(feature_names[col_num], "~ ."))
tasMLM <- rTASSEL::assocModelFitter(
  tasObj = tasGenoPheno,
  formula = formula,
  fitMarkers = TRUE,
  kinship = tasKin,
  fastAssociation = FALSE
)

stats <- tasMLM$MLM_Stats[!is.na(tasMLM$MLM_Stats$p), ]
stats <- stats[stats$p < 0.001, ]
stats <- stats[order(stats$p), ]
write.table(stats, paste("results/pMLM_", feature_names[col_num], ".tsv", sep=""), sep="\t", quote=FALSE)

effects <- tasMLM$MLM_Effects
effects <- effects[effects$Marker %in% stats$Marker, ]
write.table(effects, paste("effects/pMLM_", feature_names[col_num], ".tsv", sep=""), sep="\t", quote=FALSE)

png(paste("plots/pMLM_", feature_names[col_num], ".png", sep=""))
manhattanPlot(
  assocStats = tasMLM$MLM_Stats,
  trait      = feature_names[col_num],
  threshold  = 5
)
dev.off()