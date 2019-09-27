dataPath <- "/home/tech/Dokumente/AquaDiva/manuscript/data"
setwd(dataPath)

# Load libraries
library(phyloseq)
library(Biostrings)
library(reshape2)
library(vegan)
library(ape) # to read tree/fasta file
library(readr)
library(gridExtra)
library(ggplot2)
library(ggalluvial)
library(NbClust)
library(ggrepel)

# some settings
set_colors_mcluster <- c("#8CFF8C", "cyan", "#FBC58C", "#DFE07F", "#D0A5fE", "#008DD2")
set_colors_well <- c(
  "green1", "green4" , # Well H1 
  # "#A2A980", "wheat", # Well H2
  "cyan2", "cyan4", # well H3
  "#008DD2","gold","lightgoldenrod4",# well H4
  "#FBC58C", "mediumpurple3", "plum" # Well H5
) 

# Load data
load(file.path(dataPath, "data_rarefy.phyloseq"))

#########################################################################################
################################  hierarchical clustering  ##############################
#########################################################################################
set.seed(123)

# Calculate the dissimilarity between samples using Bray-Curtis distance
d1 <- phyloseq::distance(data_rarefy, "bray")

# To find the optimal number of clusters
res.nb <- NbClust(diss=as.dist(d1), distance=NULL,
                  min.nc = 2, max.nc = 10, 
                  method = "ward.D2", index ="silhouette") 
res.nb # print the results
nc = res.nb$Best.nc[1] # best number of clusters 

# Plot
hc1 <- hclust(d1, method="ward.D2")
clusnc = cutree(hc1, nc)
clusnc # check the clusters

# replace mCluster numbers, because we will set H41 to an extra cluster (mCluster 6) and keep the number of other clusters

clusnc[clusnc == 3] <- 10
clusnc[clusnc == 5] <- 3
clusnc[clusnc == 6] <- 5
clusnc[clusnc == 10] <- 6
clusnc

### Now plot the dendrogram and color the tips (mClusters) with different colors
par(mar=c(3, 0, 0, 0), xpd=TRUE) # bottom,left, top, right
plot(as.phylo(hc1), tip.color = set_colors_mcluster[clusnc],direction = "downwards",
     cex = 1)

# add legend to the bottom
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

add_legend("bottom", legend=c("mCluster 1","mCluster 2","mCluster 3","mCluster 4","mCluster 5","mCluster 6"), 
           fill=set_colors_mcluster ,inset=.02,
           horiz=TRUE, cex=1.1)

dev.copy(pdf, "hierarchical_clustering_per_sample.pdf", width = 14, height = 5)
dev.off()

##################################################################
#         PCoA and PERMANOVA using Bray-Curtis distance         #
##################################################################
set.seed(123)

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(data_rarefy))

# Adonis test (Number of permutations: 999)
ht_mcluster <- adonis(d1 ~ mCluster, data = sampledf) #  ***
# This output tells us that our adonis test is significant 
# so we can reject the null hypothesis that our three countries have the same centroid.
ht_mcluster

# Export the PERMANOVA results
write.table(data.frame(ht_mcluster$aov.tab), "permanova_mcluster_norm.csv", sep=";", dec=",", col.names=NA)

# Ordination (PCoA)
set.seed(123)

ordu = ordinate(data_rarefy, method = "PCoA", distance = "bray")

p_pcoa <- plot_ordination(physeq = data_rarefy, ordination = ordu) + 
  theme_bw() +
  theme(strip.text = element_text(size = 14, colour = "black", angle = 0), 
        legend.title = element_text(size = 14, colour = "black", face="bold", angle = 0),
        axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_manual(values=set_colors_well, name="Well") +
  stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=mCluster)) + 
  geom_point(size=3, aes(color = Well)) +
  ggtitle("PCoA") +
  xlab(paste0("PCoA 1 (", round(ordu$values[1,2]*100,1),"%)")) + ylab(paste0("PCoA 2 (", round(ordu$values[2,2]*100,1),"%)")) + 
  scale_fill_manual(values=set_colors_mcluster) +
  annotate("label", label = paste0("adonis(Bray-Curtis distance ~ mCluster)\n", "p < ", ht_mcluster$aov.tab$"Pr(>F)"[1], ", R-squared = ", 
                                   round(ht_mcluster$aov.tab$R2[1],2)), x = -0.45, y = -0.32, size = 4.5, 
           colour = "black", hjust = 0) +
  guides(color = guide_legend(override.aes = list(fill=NA)),
         fill= guide_legend(override.aes = list(linetype = 0)))

p_pcoa
dev.copy(pdf, "pcoa_rarefied_data.pdf", width = 8, height = 6)
dev.off()

