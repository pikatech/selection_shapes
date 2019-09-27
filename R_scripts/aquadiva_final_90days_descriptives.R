# Loading libraries
library(phyloseq)
library(reshape2)
library(ggplot2)
library(ggalluvial)

# some settings
set_colors_cluster <- c("#8CFF8C", "#8FC5FF", "#FBC58C", "#DFE07F", "#D0A5fE")
set_colors_bcluster <- c("#8CFF8C", "cyan", "pink", "#DFE07F", "#D0A5fE")
set_colors_mcluster <- c("#8CFF8C", "cyan", "#FBC58C", "#DFE07F", "#D0A5fE", "#A2D9F7")
set_colors_well <- c(
  "springgreen", "green" , # Well H1 
  # "#A2A980", "wheat", # Well H2
  "steelblue2", "blue", # well H3
  "#A2D9F7","lightgoldenrod2","peachpuff",# well H4
  "deeppink1", "mediumpurple3", "plum" # Well H5
) 
color2 <- c("Alphaproteobacteria" = "purple4", "Gammaproteobacteria"="mediumpurple2", "Deltaproteobacteria" = "mediumorchid2", "Proteobacteria_unclassified"="magenta1",
            "Actinobacteria"="brown",
            "Bacteria_unclassified"="gray90",
            "Bacteroidia"="orange1", "Ignavibacteria" = "yellow",
            "Dehalococcoidia"="lightpink",
            "Elusimicrobia" = "burlywood",
            "Nitrospira"="blue", "Thermodesulfovibrionia"="royalblue", 
            "Saccharimonadia" = "green4", "Parcubacteria"="green2", "Gracilibacteria"="darkolivegreen4", "ABY1"="lightgreen",
            "Omnitrophicaeota_cl"="darkorange1", 
            "Brocadiae"="firebrick2",
            "Others"="gray",
            "Chlamydiae" = "black"
) 

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Load data

dataPath <- "/home/tech/Dokumente/AquaDiva/manuscript/data"
setwd(dataPath)
load(file.path(dataPath, "data_rarefy.phyloseq"))

#################################################################################

#######  Bar plot: specific well to characterize the temporal patterns  #########

#################################################################################
# convert the values to relative abundance
my_data_prop <- transform_sample_counts(data_rarefy, function(x) 100*x/sum(x))

# select class as rank
RANK="Class"
my_data_prop_rank <- tax_glom(my_data_prop, taxrank=RANK)

sample_data(my_data_prop_rank)$Filter <- "Filter2" 

## Now we need to select the abundant classes (relative abundance > 1%)
my_data_prop_abun = filter_taxa(my_data_prop_rank, function(x) mean(x) > 1, TRUE)

## Calculate the abundance sum of other rare taxa (relative abundance <= 1%), We want to calculate a new phyla called "Others"
others <- tax_glom(my_data_prop_abun, "Kingdom" )
otu_table(others) <- 100 - otu_table(others)
tax_table(others)@.Data[,2:6] <- "Others"# Define all taxanomic levels as "Others"
taxa_names(others) <- "Others" # Define the taxa name as "Others"

# Combine the abundant taxa with the sum of rare taxa
OTU1 <- otu_table(rbind(otu_table(others), otu_table(my_data_prop_abun) ), taxa_are_rows = TRUE)
TAX1 <- tax_table(rbind(tax_table(others), tax_table(my_data_prop_abun) )) 
METADATA1 <- sample_data(my_data_prop_abun)

final <- phyloseq(OTU1, TAX1, METADATA1)

# However, we need to show the result as an alluvia plot using ggplot2, so we have to export the phyloseq object as data frame
my_data_rank_melt <- psmelt (final)

# order the rank by abundance
my_data_rank_melt <- transform(my_data_rank_melt, Phylum=reorder(Class, order(Abundance, decreasing=TRUE)))
# write.table(my_data_rank_melt, "~/aquadiva/redo/phylum_over1percent_per_well.csv", sep=";", dec=",", col.names=NA)

# define the x and y axis variables
is_alluvia_form(my_data_rank_melt, silent = TRUE)
my_data_rank_melt$Class<- as.factor(my_data_rank_melt$Class)

#############################################################################

# To plot
ggplot(my_data_rank_melt,
       aes(x = Samp_Cam, stratum = Class, alluvium = OTU, 
           y=Abundance,
           fill = Class, label = Class)) +
  scale_fill_manual(values=color2) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  facet_wrap(Well~.)+
  theme_bw()+
  # geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "bottom") +
  ylab("Relative abundance (%)") +
  xlab("Sampling campaign") +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0.5))+
  ggtitle("Bacterial community changes over time")

ggsave(file.path(dataPath, "alluvial_class_90days_per_well.pdf"), width = 16, height = 12)



