# Load data
dataPath <- "/home/tech/Dokumente/AquaDiva/manuscript/data"
load(file.path(dataPath, "data_rarefy.phyloseq"))

# load libraries
library(phyloseq)
library(reshape2)
library(readr)
require(gridExtra)
library(ggplot2)
library(ggalluvial)
library(FSA)
library(rcompanion) #To have R convert this table to a compact letter display for us
library(broom)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Hmisc)
library(gplots)

# some settings
set_colors_mcluster <- c("#8CFF8C", "cyan", "#FBC58C", "#DFE07F", "#D0A5fE", "#A2D9F7")
set_colors_well <- c(
  "green1", "green4" , # Well H1 
  # "#A2A980", "wheat", # Well H2
  "cyan1", "cyan4", # well H3
  "#A2D9F7","lightgoldenrod1","lightgoldenrod4",# well H4
  "#FBC58C", "mediumpurple3", "plum" # Well H5
) 

set_colors <- c("Alphaproteobacteria" = "purple4", "Gammaproteobacteria"="mediumpurple2", "Deltaproteobacteria" = "mediumorchid2", "Proteobacteria_unclassified"="magenta1",
                "Actinobacteria"="brown","Thermoleophilia"="brown3","Acidimicrobiia"="brown4","Actinobacteria_unclassified"="brown2",
                "Bacteria_unclassified"="gray90",
                "Bacteroidia"="orange1", "Ignavibacteria" = "yellow",
                "Dehalococcoidia"="lightpink", "KD4-96"="hotpink", "Anaerolineae" = "pink1",  "JG30-KF-CM66"="pink2", # Chloroflexi
                "Elusimicrobia" = "burlywood","Elusimicrobia_unclassified" ="yellow4",
                "Nitrospira"="blue", "Thermodesulfovibrionia"="royalblue", "HDB-SIOI1093"="blue3", # Nitrospirae
                "Saccharimonadia" = "green4", "Parcubacteria"="green2", "Gracilibacteria"="darkolivegreen4", "ABY1"="lightgreen",
                "Patescibacteria_unclassified" ="greenyellow", "Berkelbacteria"="#5CA595",
                "Omnitrophicaeota_cl"="darkorange2", "Omnitrophia"="chocolate",
                "Brocadiae"="firebrick2","OM190"="brown1", "Phycisphaerae"="orangered3",
                "Others"="gray",
                "Gemmatimonadetes" = "black",
                "Chlamydiae"="bisque",
                "Kiritimatiellae"="chocolate2",
                "Lineage_llc"="blue2", "Verrucomicrobiae"="red2", "Deinococci" = "yellow2",
                "MVP-15" ="darkturquoise",
                "Blastocatellia_(Subgroup_4)"="cyan2","Holophagae"="cornflowerblue", "Subgroup_6"="cyan1", # acidobacteria
                "Planctomycetacia" = "deepskyblue4",
                "vadinHA49"="mediumpurple4",
                "Nitrospinia"="violetred4",
                "Clostridia"="gray", # Firmicutes
                "NC10"="yellow2", # Rokubacteria 
                "Margulisbacteria_cl"="yellow4",
                "Sericytochromatia"="lightcyan", 
                "Babeliae" = "grey70",
                "Cyanobacteria_unclassified" = "#00FFBF"            
)

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


setwd(dataPath)
#####################################################################################
##############################         core species   ###############################
#####################################################################################

# Definition: OTUs that occur only over 80% of the samples

my_data_core <- data_rarefy %>%
 filter_taxa(function(x) sum(x > 0) >= (0.8*length(x)), TRUE)


# Check the relative abundance of those core OTUs
my_data_prop <- data_rarefy %>%
  transform_sample_counts(function(x) 100*x/sum(x))  
save(my_data_core, file ="my_data_core_26.phyloseq")


###############################################################
##############     Bar plot --- For only core OTUs    #########
###############################################################

my_data_prop_core <- subset_taxa(my_data_prop, rownames(tax_table(my_data_prop)) %in% taxa_names(my_data_core)) 
### Merge the core OTUs at the class level
RANK="Class"# relative abundance
my_data_prop_rank <- tax_glom(my_data_prop_core, taxrank=RANK)

### To calculate the mean abundance of the core OTUs (at the class level) in each well
data_rank_melt <- psmelt(my_data_prop_rank)
# order the rank by abundance
data_rank_melt <- transform(data_rank_melt, Phylum=reorder(Phylum, order(Abundance, decreasing=TRUE)))

test <- psmelt(my_data_prop_core)
mean_core <- aggregate(Abundance ~ Well+Class, test, function(x) c(mean = mean(x), sd = sd(x)))
mean_core 
write.table(mean_core, "coreOTU_abundance_AVE_class.csv", sep=";", dec=",", col.names = NA)

#### Statistical test to see the abudnance difference of the core OTUs between different wells
PT<-dunnTest(Abundance~Well,data_rank_melt, method="bh")
PT<-PT$res

cldList(comparison = PT$Comparison,
        p.value    = PT$P.adj, # or PT$P.unadj
        threshold  = 0.05)
####################################################################
data_rank_melt_select  <- aggregate(Abundance ~ Well+Class+Phylum, data_rank_melt, function(x) c(mean = mean(x), sd = sd(x)))
data_rank_melt_select$Class<- as.factor(data_rank_melt_select$Class)

data_rank_melt_select$Abundance.mean <- data_rank_melt_select$Abundance[,1]

####### Combine average and per well and plot 
data_rank_melt_ave1 <-  aggregate(Abundance ~ Class+ Phylum, data_rank_melt, function(x) c(mean = mean(x), sd = sd(x)))
data_rank_melt_ave1$Well <- "ave."

data_rank_melt_ave1<-data_rank_melt_ave1[, c(4,1:3)]
data_rank_melt_ave1$Abundance.mean <- data_rank_melt_ave1$Abundance[,1]

# Combine the datasets and make a new data fram for plotting
df <- rbind(data_rank_melt_ave1,data_rank_melt_select)

################################################################################################
################  Plot the relative abundance of the 26 core OTUs at the class level ###########
################################################################################################

ggplot() + geom_bar(aes(y = Abundance.mean, x = Well, 
                        fill = Class, 
                       # fill = Phylum, 
                        group=Phylum), data = df,
                    stat="identity", position = "stack", lwd = 0,
                    width = 0.6, colour =  "black")  +    
  scale_fill_manual(values=set_colors) +
  theme_bw() +
  theme(legend.position = "right") +
  ylab("Abundance") +
  xlab("Well") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

ggsave("bar_26core_absolute_abundance_class_Phylum.pdf", width = 8, height = 4)

### To calculate the total abundance sum of core OTUs in each well showed in this plot
hh <- data.frame(sample_sums(my_data_prop_rank)) # calculate the sum of core per sample
names(hh)<- "Abundance"
hh$Well <- sample_data(my_data_prop_rank)$Well
mean_core1 <- aggregate(Abundance ~ Well, hh, function(x) c(mean = mean(x), sd = sd(x)))
write.table(mean_core1, "total_core_abundance_per_well.csv", sep=";", dec=",", col.names = NA)


###########################################################################################
######                             Statistical testing                               ######
###########################################################################################

dat <- my_data_prop_core %>%
  tax_glom(taxrank=RANK) %>%
  psmelt() # the lower the taxonomic level, the slower

RANK <- "Class"

### Use Kruskal test to identify the classes (relative abundance sum of core OTUs belong to this class)

dat %>%
  group_by_(RANK) %>%
  mutate(mean=mean(Abundance)) %>%
  do(tidy(kruskal.test(.$Abundance~.$Well))) %>%
  ungroup() %>%
  mutate(p.adjust=p.adjust(p.value), method="BH") -> results_quantitative

# Table of signficant taxa
data.frame(results_quantitative) %>%
  subset(p.adjust<0.05) -> sigtab

sigtab
write.table(sigtab, "sigtab_kruskal_BH_adjusted_class_well.txt", sep="\t", col.names = NA)


#### Statistical test of the selected OTUs, respectively
class <- "Nitrospira" 
#class <- "Thermodesulfovibrionia"
#class <- "Brocadiae" 
#class <- "Alphaproteobacteria" 

PT3 <- dunnTest(Abundance~Well,subset(dat, Class==class), method="bh")
PT3 <- PT3$res

cldList(comparison = PT3$Comparison,
        p.value    = PT3$P.adj, # or PT$P.unadj
        threshold  = 0.05)

###############################################################################################
################## Correlations between core taxa and Environmental variables #################
###############################################################################################

# import processed environmental data
df1 <- read.csv(file.path(dataPath, "env_data_NA_omit.csv"), sep=";", dec = ",", row.names=1, header=TRUE)
df2 <- df1[,16:32]
sample_data(my_data_core) <- df2

RANK <- "Class" # we can substitute the taxonomic levels and re-run the analysis
datt <- my_data_core %>%
  tax_glom(taxrank=RANK) %>%
  transform_sample_counts(function(x) 100*x/sum(x)) %>%
  psmelt() # the lower the taxonomic level, the slower

head(datt)
variable.names(datt)

## Select the variables for correlation analysis
datt1 <- datt[,c(1,3:20)] 

#########################################################################
values <- colnames(datt1)[3:19]
datt1$OTU <- as.factor(datt1$OTU)
corVals <- matrix(0, nrow = length(levels(datt1$OTU)), ncol = length(values))
corVals.p <- matrix(0, nrow = length(levels(datt1$OTU)), ncol = length(values))
rownames(corVals) <- levels(datt1$OTU)
colnames(corVals) <- values
rownames(corVals.p) <- levels(datt1$OTU)
colnames(corVals.p) <- values
for (i in values) {
  for (j in levels(datt1$OTU)) {
    corVals.p[j, i] = cor.test(datt1[datt1$OTU == j, "Abundance"], datt1[datt1$OTU == j, i], method = "spearman")$p.value
    corVals[j, i] = cor(datt1[datt1$OTU == j, "Abundance"], datt1[datt1$OTU == j, i], method = "spearman")
  }
}

corVals.p[] <- p.adjust(corVals.p, method="BH")

coe <- melt(corVals)
pv <- melt(corVals.p)

names(coe) <- c("OTU", "Parameter", "Spearman")
names(pv) <- c("OTU", "Parameter", "adjusted.p")

coe$adjust.p <- pv$adjusted.p
coe1 <- subset(coe, adjust.p<0.05)

round(coe1, 4)
write.table(coe1, "spearman_class_env.csv", sep=";", dec=",", col.names=NA)

p_cor <- ggplot(coe1, aes(Parameter, OTU)) +
  geom_tile(aes(fill = Spearman)) + 
  scale_fill_gradient(low = "purple", high = "yellow") +
  theme_classic() +
  scale_x_discrete(position = "top") +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) 
p_cor

ggsave("26heatmap_class_env_at_otu_level.pdf", width=8, height=5)


