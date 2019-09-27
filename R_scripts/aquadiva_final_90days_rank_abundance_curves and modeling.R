dataPath <- "/home/tech/Dokumente/AquaDiva/manuscript/data"
setwd(dataPath)

### load libraries
library(phyloseq)
library(reshape2)
library(vegan) # import adonis
library(readr)
library(gridExtra)
library(ggplot2)
library(goeveg)
library(magrittr)
library(ggpubr)
library(RADanalysis)
library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions

# some settings
set_colors_mcluster <- c("#8CFF8C", "cyan", "#FBC58C", "#DFE07F", "#D0A5fE", "#A2D9F7")

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Load data
load(file.path(dataPath, "data_rarefy.phyloseq"))

#############################################################################################
rads <- vegan_otu(data_rarefy)

### raw curves for visualization
line_cols <- set_colors_mcluster
sample_classes <- sample_data(data_rarefy)$mCluster

plot(1, xlim = c(1,10000), ylim = c(1,1000), col = "white", log = "xy", axes = FALSE,
     xlab = "Rank", ylab = "Abundance", main = "")
sfsmisc::eaxis(side = 1, at = c(1,10,100,1000, 10000))
sfsmisc::eaxis(side = 2)
for (i in 1:nrow(rads)) {
  temp <- sort(rads[i,], decreasing = TRUE)
  temp <- temp[temp>0]
  lines(x = temp, lwd = 2, col = line_cols[sample_classes[i]])
}
legend("topright", bty = "n", legend = c("mCluster 1","mCluster 2","mCluster 3",
                                         "mCluster 4", "mCluster 5", "mCluster 6"),
       col = line_cols, lwd = 2)

### Normalization ###
richness <- sapply(X = 1:nrow(rads), function(i) length(which(rads[i,] > 0)))
boxplot(richness, horizontal = TRUE, xlab = "Richness")
quantile(richness)

########################################################
##########   Normalization using MaxRank  ##############
########################################################
nrads <- RADnormalization_matrix(input = rads,max_rank = min(richness), average_over = 100,
                                 sample_in_row = TRUE,verbose = FALSE)
nrads <- nrads$norm_matrix

maxrank <- min(richness)

plot(1e10,xlim = c(1,maxrank),ylim = c(1e-4,1),log="xy",
     xlab = "Rank",ylab = "Relative abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,10,100,1000))
sfsmisc::eaxis(side = 2,at = c(1e-4,1e-3,1e-2,1e-1,1),las = 0)

# plot confidence intervals of representative nrads 
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 1"),
                        plot = TRUE,confidence = 0.95,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.3),border = NA)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 2"),
                        plot = TRUE,confidence = 0.95,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.3),border = NA)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 3"),
                        plot = TRUE,confidence = 0.95,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.3),border = NA)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 4"),
                        plot = TRUE,confidence = 0.95,with_conf = TRUE,
                        col = scales::alpha(line_cols[4],0.3),border = NA)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 5"),
                        plot = TRUE,confidence = 0.95,with_conf = TRUE,
                        col = scales::alpha(line_cols[5],0.3),border = NA)
a <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 6"),
                        plot = TRUE,confidence = 0.95,with_conf = TRUE,
                        col = scales::alpha(line_cols[6],0.3),border = NA)
# plot representative nrads
a1 <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 1"),
                         plot = TRUE,with_conf = FALSE,
                         col = scales::alpha(line_cols[1],1),lwd = 4)
a2 <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 2"),
                         plot = TRUE,with_conf = FALSE,
                         col = scales::alpha(line_cols[2],1),lwd = 4)
a3 <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 3"),
                         plot = TRUE,with_conf = FALSE,
                         col = scales::alpha(line_cols[3],1),lwd = 4)
a4 <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 4"),
                         plot = TRUE,with_conf = FALSE,
                         col = scales::alpha(line_cols[4],1),lwd = 4)
a5 <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 5"),
                         plot = TRUE,with_conf = FALSE,
                         col = scales::alpha(line_cols[5],1),lwd = 4)
a6 <- representative_RAD(norm_rad = nrads,sample_ids = which(sample_classes == "mCluster 6"),
                         plot = TRUE,with_conf = FALSE,
                         col = scales::alpha(line_cols[6],1),lwd = 4)

legend("bottomleft",bty = "n",legend = c("mCluster 1","mCluster 2","mCluster 3", "mCluster 4", "mCluster 5", "mCluster 6"),
       col = line_cols,lwd = 4)


dev.copy(png, file="rank_abundance_normalized.png", width=1000, height=800, res=144)
dev.off()

####################################################################################################################

#######################################
##############  MDS plot  #############
#######################################

d <- dist(x = nrads, method = "manhattan")
# ordination using classical multi-dimensional scaling
mds <- cmdscale(d = d, k = 2, eig = TRUE)
mds$eig[1]/sum(mds$eig)
mds$eig[2]/sum(mds$eig)
# plot the points 
sample_classes <- sample_data(data_rarefy)$mCluster
line_cols <- set_colors_mcluster
plot(mds$points, xlab = paste0("MDS1 (", round(mds$eig[1]/sum(mds$eig)*100,2), "% explained variance)"),
                 ylab = paste0("MDS2 (", round(mds$eig[2]/sum(mds$eig)*100,2), "% explained variance)"),
                 pch = 19, cex =0.5, col = line_cols[sample_classes], main = NULL)

# add the representative points with error bar to the previous plot
aa1 <- representative_point(input = mds$points,ids = which(sample_classes == "mCluster 1"),
                            col = scales::alpha(line_cols[1],0.5),
                            plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
aa2 <- representative_point(input = mds$points,ids = which(sample_classes == "mCluster 2"),
                            col = scales::alpha(line_cols[2],0.5),
                            plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
aa3 <- representative_point(input = mds$points,ids = which(sample_classes == "mCluster 3"),
                            col = scales::alpha(line_cols[3],0.5),
                            plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
aa4 <- representative_point(input = mds$points,ids = which(sample_classes == "mCluster 4"),
                            col = scales::alpha(line_cols[4],0.5),
                            plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
aa5 <- representative_point(input = mds$points,ids = which(sample_classes == "mCluster 5"),
                            col = scales::alpha(line_cols[5],0.5),
                            plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
aa6 <- representative_point(input = mds$points,ids = which(sample_classes == "mCluster 6"),
                            col = scales::alpha(line_cols[6],0.5),
                            plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)

dev.copy(pdf, file="rank_abundance_mds.pdf", width=5, height=3.5)
dev.off()

#### Test for significance between mClusters
adonis(d ~ sample_classes, perm=999)

#################################################################################
###################   Model fitting  for each mCluster   ########################
#################################################################################

### Test using radfit to see which model suits our data best
# Use a1, a2, a3, a4, a5 and a6, respectively
mod <- radfit(a1$average, family=gaussian)
mod

## We chose to use zipf model, because it has less deviance and lower AIC values
mod1 <- rad.zipf(a1$average, family=gaussian)
mod2 <- rad.zipf(a2$average, family=gaussian)
mod3 <- rad.zipf(a3$average, family=gaussian)
mod4 <- rad.zipf(a4$average, family=gaussian)
mod5 <- rad.zipf(a5$average, family=gaussian)
mod6 <- rad.zipf(a6$average, family=gaussian)

## Check the results to report (p1, gamma, Deviance, AIC, BIC values)
mod1
mod2
mod3
mod4
mod5
mod6

## For plotting the rank abundance patterns and fitted values ###
ddt1 <- data.frame(Rank = 1:maxrank, fitted.abundance = mod1$fitted.values, Abundance = a1$average, mCluster = "mCluster 1")
ddt2 <- data.frame(Rank = 1:maxrank, fitted.abundance = mod2$fitted.values, Abundance = a2$average, mCluster = "mCluster 2")
ddt3 <- data.frame(Rank = 1:maxrank, fitted.abundance = mod3$fitted.values, Abundance = a3$average, mCluster = "mCluster 3")
ddt4 <- data.frame(Rank = 1:maxrank, fitted.abundance = mod4$fitted.values, Abundance = a4$average, mCluster = "mCluster 4")
ddt5 <- data.frame(Rank = 1:maxrank, fitted.abundance = mod5$fitted.values, Abundance = a5$average, mCluster = "mCluster 5")
ddt6 <- data.frame(Rank = 1:maxrank, fitted.abundance = mod6$fitted.values, Abundance = a6$average, mCluster = "mCluster 6")

ddt <- rbind(ddt1, ddt2, ddt3, ddt4, ddt5, ddt6)

m1 <- ggscatter(ddt, x = "Rank", y = "Abundance", shape=1, size=2, facet.by = "mCluster") +
  geom_line(data=ddt, aes(x=Rank, y=fitted.abundance), color='red', size=1) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14),
        strip.text.x = element_text(size = 14, colour = "black",face="bold", angle = 0),
        legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) + 
  ylab("Relative abundance")

m1

ggsave(m1, file="rank_abundance_comparion_zipf_gaussian_normalized.pdf", width=12, height=6)
dev.copy(png, file="rank_abundance_comparion_zipf_gaussian_normalized.png", width=1500, height=900, res=144)
dev.off()

