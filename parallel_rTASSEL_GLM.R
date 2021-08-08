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

formula <- as.formula(paste(feature_names[col_num], "~ ."))

tasGLM <- rTASSEL::assocModelFitter(
  tasObj = tasGenoPheno,
  formula = formula,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = FALSE
)

write.table(tasGLM$GLM_Stats[order(tasGLM$GLM_Stats$p), ], paste("results/parallelMLM_", feature_names[col_num], ".tsv", sep=""), sep="\t", quote=FALSE)

png(paste("plots/parallelMLM_", feature_names[col_num], ".png", sep=""))
manhattanPlot(
  assocStats = tasGLM$GLM_Stats,
  trait      = feature_names[col_num],
  threshold  = 5
)
dev.off()