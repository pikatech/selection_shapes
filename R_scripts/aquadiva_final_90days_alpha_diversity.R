
# load libraries
library(phyloseq)
library(Biostrings)
library(reshape2)
library(ape) # to read tree/fasta file
library(readr)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(picante)
library(FSA) # dunn test
library(rcompanion) #To have R convert this table to a compact letter display for us

# some settings
set_colors_mcluster <- c("#8CFF8C", "cyan", "#FBC58C", "#DFE07F", "#D0A5fE", "#A2D9F7")

set_colors_well <- c(
  "green1", "green4" , # Well H1 
  # "#A2A980", "wheat", # Well H2
  "cyan1", "cyan4", # well H3
  "#A2D9F7","lightgoldenrod1","lightgoldenrod4",# well H4
  "#FBC58C", "mediumpurple3", "plum" # Well H5
)

# Load data
dataPath <- "/home/tech/Dokumente/AquaDiva/manuscript/data"
load(file.path(dataPath, "aquadiva_DNA_90days.phyloseq")
load(file.path(dataPath, "data_rarefy.phyloseq"))

####################################################################
##################### rarefaction curve in phyloseq ################
####################################################################
#https://github.com/joey711/phyloseq/issues/143

setwd(dataPath)
psdata <- my_data90
min(sample_sums(psdata))
set.seed(123)

calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require(plyr) # ldply
  require(reshape2) # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'value')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(psdata, c('Observed', 'Shannon',"InvSimpson", "Chao1"), rep(c(1:20*5, 2:75 * 400), each = 10))
summary(rarefaction_curve_data)
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(value), Alpha_diversity_sd = sd(value))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')

View(rarefaction_curve_data_summary_verbose)


ggplot(data = subset(rarefaction_curve_data_summary_verbose, Measure !="se.chao1"),
       mapping = aes(
         x = Depth,
         y = Alpha_diversity_mean,
         colour = Well
       )) + 
  geom_smooth(method = "loess") +
  geom_vline(xintercept = 3000, linetype = "solid", size=1) +
  facet_wrap(facets = ~ Measure,
             scales = 'free_y') +
  theme_bw() +
  theme(strip.text = element_text(size = 18, colour = "black", angle = 0), 
        legend.title = element_text(size = 18, colour = "black", face="bold", angle = 0),
        axis.title=element_text(size=18), axis.text=element_text(size=18), legend.text=element_text(size=18)) +
  ylab("Alpha diversity") + 
  scale_color_manual(values=set_colors_well)

ggsave(file.path(dataPath, "rarefaction_well.pdf"), width = 20, height = 10)
#################################################################################


#######################################################################
############################  rarefy Data  ############################
#######################################################################

# According to the rarefaction analysis (shannon values stablized after 3000), 
# we decide to nomalize the library size at 3000 sequences/sample for downstream data analysis

set.seed(123)
data_rarefy <- rarefy_even_depth(my_data90, sample.size = 3000, verbose = FALSE, replace = TRUE)

### integrate the phylogenetic tree, which was output from mothur, for phylogenetic diversity estimation #####
tree <- read.tree("data_rarefy.phylip.tre") 
phytree <- phy_tree(tree)

data_rarefy <- merge_phyloseq(data_rarefy, phytree)
save(data_rarefy, file ="data_rarefy.phyloseq")


#######################################################################
#####################   alpha diversity indices #######################
#######################################################################

# 1) Estimate Shannon index
alphadiv <- estimate_richness(data_rarefy, measures="Shannon")
sample_data(data_rarefy)$Shannon <- alphadiv$Shannon

# 2)  Estimate Phylogenetic diversity 

estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}
pd <- estimate_pd(data_rarefy)
sample_data(data_rarefy)$PD <- pd$PD

### Plot the diversity indices
summary <- data.frame(sample_data(data_rarefy))

p <- ggboxplot(summary, x="Well", y = "Shannon", fill = "mCluster",
               palette = set_colors_mcluster,
               add = "jitter" 
) +
  stat_compare_means(label.x = 1, label.y = 8) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) 

p

ggsave(file.path(dataPath, "shannon_mcluster.pdf"), width = 7, height = 5)


###############################################################################

p1 <- ggboxplot(summary, x="Well", y = "PD", fill = "mCluster",
                palette = set_colors_mcluster,
                add = "jitter"
) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) +
  stat_compare_means(label.x = 1, label.y = 8) 

p1

ggsave(file.path(dataPath, "PD_well_boxplot.eps"), width = 7, height = 5)

# To calculate mean and SD of the alpha diverity indices
res_shannon <- aggregate(Shannon ~ mCluster, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res_shannon

res_pd <- aggregate(PD ~ mCluster, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res_pd


############################ Dunn test for multiple comparisons of alpha diversity indices ####################
# Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations

# 1) Shannon 
PX <-dunnTest(Shannon~mCluster,summary, method="bh")
PX <-PX$res

cldList(comparison = PX$Comparison,
        p.value    = PX$P.adj, # or PT$P.unadj
        threshold  = 0.05)
# To have R convert this table to a compact letter display for us

PX2<-dunnTest(Shannon~Well,summary, method="bh")
PX2<-PX2$res

cldList(comparison = PX2$Comparison,
        p.value    = PX2$P.adj, # or PT$P.unadj
        threshold  = 0.05)


# 2) PD 
PT<-dunnTest(PD~mCluster,summary, method="bh")
PT<-PT$res

cldList(comparison = PT$Comparison,
        p.value    = PT$P.adj, # or PT$P.unadj
        threshold  = 0.05)
# To have R convert this table to a compact letter display for us

PT2<-dunnTest(PD~Well,summary, method="bh")
PT2<-PT2$res

cldList(comparison = PT2$Comparison,
        p.value    = PT2$P.adj, # or PT$P.unadj
        threshold  = 0.05)

