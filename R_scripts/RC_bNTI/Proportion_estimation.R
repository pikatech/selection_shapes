###### Donuts plot #####
#
# x:      numeric vector for each slice
# group:  vector identifying the group for each slice
# labels: vector of labels for individual slices
# col:    colors for each group
# radius: radius for inner and outer pie (usually in [0,1])

donuts <- function(x, group = 1, labels = NA, col = NULL, radius = c(.7, 1)) {
  group <- rep_len(group, length(x))
  ug  <- unique(group)
  tbl <- table(group)[order(ug)]
  
  col <- if (is.null(col))
    seq_along(ug) else rep_len(col, length(ug))
  col.main <- Map(rep, col[seq_along(tbl)], tbl)
  col.sub  <- lapply(col.main, function(x) {
    al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    Vectorize(adjustcolor)(x, alpha.f = al)
  })
  
  plot.new()
  
  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L],
      col = unlist(col.sub), labels = labels)
  
  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L],
      col = unlist(col.main), labels = NA)
}

################### Proportion estimation ################################
path <- "path/to/assembly/data"
setwd(path)
species <- "pnk75"

weighted.bNTI <- read.table(file.path(path, species, "bNTI/weighted_bNTI.csv"),
                                      sep=",", dec=".", row.names=1, header=TRUE)

dim(weighted.bNTI)
weighted.bNTI1 <- weighted.bNTI[lower.tri(weighted.bNTI)] # the lower triangular matrix
hist(weighted.bNTI1) # check it's output
length(weighted.bNTI1) # the value should be (n*(n-1))/2

##################################
rc <- read.table(file.path(path, species, "RC/rc_bray.csv"),
                           sep=",", dec=".", row.names=1, header=TRUE)
dim(rc)

# Keep only the half (upper) triangular matrix
rc <- rc[lower.tri(rc)]
hist(rc)
length(rc) # the value should be (n*(n-1))/2

##################################
#1) Homogenizing dispersal
homo_disp <- length(which(abs(weighted.bNTI1) < 2 & rc < -0.95)) / length(rc)

# 2) undominated
undo <- length(which(abs(weighted.bNTI1) <= 2 & abs(rc) <= 0.95)) / length(rc)

# 3) dispersal limitation
disp_limi <- length(which(abs(weighted.bNTI1) < 2 & rc > 0.95)) / length(rc)

# 4) variable selection
vari_sele <- length(which(weighted.bNTI1 > 2)) / length(rc)

# 5) Homogeneous selection
homo_sele <- length(which(weighted.bNTI1 < -2)) / length(rc)

# Check if the calculation is correct
Proportion <- c(homo_disp, disp_limi, undo, vari_sele, homo_sele) 
sum(Proportion) # this value should be 1

Mechanism <- c("HD",  "DL", "U", "VS", "HS")
Process <- c("Total stochasticity", "Total stochasticity", "Total stochasticity",
             "Total selection", "Total selection")

pdf(paste0("plot_", species, ".pdf"), width=6, height=4)
donuts(Proportion, Process, sprintf('%s: %2.2f%%', Mechanism, Proportion*100), col = c("blue", "red"))
title(main = species)
dev.off()

