# Load data
dataPath <- "/home/tech/Dokumente/AquaDiva/manuscript/data"
load(file.path(dataPath, "data_rarefy.phyloseq"))

# load libraries
library(phyloseq)
library(reshape2)
library(FSA)
library(rcompanion) # To have R convert this table to a compact letter display for us
library(dplyr)
library(ggplot2)
library(ggpubr) # ggboxplot
library(magrittr)
library(Hmisc)

##################################################################################
## Rare OTUs were defined as the relative abundance < 0.01% in a dataset (here, in all the samples)
my_data_rare_all_samples <- data_rarefy %>%
  transform_sample_counts(function(x) 100*x/sum(x))%>% 
  filter_taxa( function(x) mean(x) < 0.01, TRUE)

rareotuall <- data.frame(otu_table(my_data_rare_all_samples)@.Data)

##################################################################################
# Subset the dataset to each mCluster before identifying the rare species
## Rare OTUs were defined as the relative abundance < 0.01% in a dataset (subset by each mCluster)
my_data_prop <- data_rarefy %>%
  transform_sample_counts(function(x) 100*x/sum(x))

list_of_data_frames = list()
for (i in 1:6) {
  mcl <- paste("mCluster", i)

  my_data_rare <- data_rarefy %>%
    subset_samples(mCluster == mcl) %>% # mcl
    filter_taxa( function(x) sum(x) > 0, TRUE) %>%
    transform_sample_counts(function(x) 100*x/sum(x))%>% 
    filter_taxa( function(x) mean(x) < 0.01, TRUE)

  my_data_prop_rare <- subset_taxa(my_data_prop, rownames(tax_table(my_data_prop)) %in% taxa_names(my_data_rare)) 
  my_data_prop_rare <- subset_samples(my_data_prop_rare, mCluster == mcl)

  rareotu <- data.frame(otu_table(my_data_prop_rare)@.Data)
  rareotu$OTUid <- rownames(rareotu)

  list_of_data_frames = c(list_of_data_frames, list(rareotu))
}

df_rare = Reduce(function(...) merge(..., all=T), list_of_data_frames)
row.names(df_rare) <- df_rare[,1]
df_rare <- df_rare[,-1]
df_rare[is.na(df_rare)] <- 0

## Put it into phyloseq object
OTU <- otu_table(as.matrix(df_rare), taxa_are_rows = TRUE)
rare_data <- phyloseq(OTU,tax_table(data_rarefy),sample_data(data_rarefy) )
rare_data <- filter_taxa(rare_data, function(x) sum(x) > 0, TRUE) 
sample_data(rare_data) <- sample_data(data_rarefy)

sample_sums(rare_data)
save(rare_data, file ="rare_data_mcl.phyloseq")

######################################################################################
# load("rare_data_mcl.phyloseq")
  
# Calculate the average relative abundance of rare OTUs in each dataset (mCluster)
my_data_rare_total <- rare_data %>%
  tax_glom("Kingdom") %>%
  psmelt()

me <- aggregate(Abundance ~ mCluster, my_data_rare_total, function(x) c(mean = mean(x), sd = sd(x)))
me <- round(me, 1)

write.table(me, "rare_ave_total_per_mCluster.csv", sep=";", dec=",", col.names = NA)


p2 <- ggboxplot(my_data_rare_total,  x = "mCluster", y = "Abundance",
                fill = "#C5C6C6", color="black",
                add = "none"
) +
  geom_hline(yintercept = 25.9, linetype = "dashed", color = "red") +
  stat_compare_means(aes(label =paste0("Kruskal-Wallis: ", ..p.signif..)), label.x = 3) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  theme(text=element_text(size=14)) +
  ylab("Relative Abundance (%)")
p2

#### To see if there is abundance difference of the rare species between mCLuster
library(FSA)
PT4 <- dunnTest(Abundance~mCluster, my_data_rare_total, method="bh")
PT4 <- PT4$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT4$Comparison,
        p.value    = PT4$P.adj, # or PT$P.unadj
        threshold  = 0.05)

