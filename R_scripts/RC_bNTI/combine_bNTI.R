library(ape)     # read.tree
library(picante) # match.phylo.data
library(foreach)
library(abind)

# Set up the paths to the files
basePath <- "path/to/assembly/data"
rawPath <- file.path(basePath, "10xHT/HTU")
csvFile <- file.path(rawPath, "HTU_otutable.csv")
phyloFile <- file.path(rawPath, "HTU.tre")
calcPath <- file.path(rawPath, "bNTI")
bNTIfiles <- dir(calcPath, pattern="idiv_102_.*", full.names=T)

# Combine calculated randomized betaMNTD
rand.weighted.bMNTD.comp <- foreach (b = bNTIfiles, .combine='abind', .multicombine=TRUE) %do% {
   readRDS(b)
}

dim(rand.weighted.bMNTD.comp)

# Read in data
otutable <- read.csv(csvFile, sep=";", dec = ",", row.names=1, header=TRUE)
phylo <- read.tree(phyloFile)
match.phylo.otutable <- match.phylo.data(phylo, t(otutable))

# Calculate empirical betaMNTD
beta.mntd.weighted <- as.matrix(comdistnt(t(match.phylo.otutable$data), cophenetic(match.phylo.otutable$phy), abundance.weighted=T))
write.csv(beta.mntd.weighted, file.path(calcPath, 'betaMNTD_weighted.csv'), quote=F)

# Use randomized values to weight bNTI
weighted.bNTI = matrix(c(-Inf),nrow=ncol(match.phylo.otutable$data), ncol=ncol(match.phylo.otutable$data))
dimnames(weighted.bNTI) <- dimnames(beta.mntd.weighted)
for (columns in 1:(ncol(match.phylo.otutable$data)-1)) {
   for (rows in (columns+1):ncol(match.phylo.otutable$data)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,]
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
      rm("rand.vals")
   }
}

# Store results
write.csv(weighted.bNTI, file.path(calcPath, "weighted_bNTI.csv"), quote=F);

pdf(file.path(calcPath, "weighted_bNTI_Histogram.pdf"))
hist(weighted.bNTI)
dev.off()

