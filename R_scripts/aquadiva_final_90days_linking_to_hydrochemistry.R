dataPath <- "/home/tech/Dokumente/AquaDiva/manuscript/data"
setwd(dataPath)

# Load libraries
library(phyloseq)
library(reshape2)
library(readr)
library(gridExtra)
library(ggplot2)
library(magrittr)
library(digest)

set_colors_cluster <- c("#8CFF8C", "#8FC5FF", "#FBC58C", "#DFE07F", "#D0A5fE")
set_colors_bcluster <- c("#8CFF8C", "cyan", "pink", "#DFE07F", "#D0A5fE")
set_colors_mcluster <- c("#8CFF8C", "cyan", "#FBC58C", "#DFE07F", "#D0A5fE", "#A2D9F7")
set_colors_well <- c(
 "H13"= "green1", "H14" ="green4" , # Well H1 
  # "#A2A980", "wheat", # Well H2
"H31"=  "cyan1","H32"= "cyan4", # well H3
"H41"=  "#A2D9F7", "H42"= "lightgoldenrod1","H43" = "lightgoldenrod4",# well H4
"H51"=  "#FBC58C","H52"= "mediumpurple3","H53"= "plum" # Well H5
) 

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

######################################################################
# Load data
load(file.path(dataPath, "data_rarefy.phyloseq"))

###########################################################
####   Data preparation: hydrochemical parameters    ######
###########################################################

dat_all <- data.frame(sample_data(data_rarefy))
u <- dat_all # for all env samples

## Remove variables that contain too many missing or 0 values: NO2.1, FE, MN, FE.1, MN.1
variable.names(u)

u1 <- u[ , -which(names(u) %in% c("NO2.1","FE", "MN", "FE.1", "MN.1", "P", "PO4"))]

# Remove the unnecessary variables
variable.names(u1)

u2 <- u1[ , -which(names(u1) %in% c("EC25","air_temperature__minimum", "air_temperature__maximum", "sunshine_duration",
                                   "soil_temperature__minimum","relative_humidity__minimum","relative_humidity__mean","relative_humidity__maximum",
                                   "wind_speed__mean","gust_speed__maximum", "Observed","Chao1","se.chao1", "ACE","se.ACE", "Simpson", "InvSimpson",
                                   "Fisher","Cluster" ))]

# select the necessary variables
variable.names(u2)
df <-subset(u2, select=EC:Shannon)

#### Change the variable names
variable.names(df)
colnames(df)[which(names(df) == "PH")] <- "pH"
colnames(df)[which(names(df) == "FE2")] <- "Iron2"
colnames(df)[which(names(df) == "NO2")] <- "Nitrite"
colnames(df)[which(names(df) == "NH4")] <- "Ammonium"
colnames(df)[which(names(df) == "NO3")] <- "Nitrate"
colnames(df)[which(names(df) == "PO4")] <- "Orthophosphate"
colnames(df)[which(names(df) == "SO4")] <- "Sulphate"
colnames(df)[which(names(df) == "Cl")] <- "Chloride"
colnames(df)[which(names(df) == "CA")] <- "Calcium"
colnames(df)[which(names(df) == "K")] <- "Potassium"
colnames(df)[which(names(df) == "MG")] <- "Magnesium"
colnames(df)[which(names(df) == "NA.")] <- "Sodium"
colnames(df)[which(names(df) == "P")] <- "Phosphorus"
colnames(df)[which(names(df) == "S")] <- "Sulfur"
colnames(df)[which(names(df) == "WL_mamsl")] <- "Water_level"
colnames(df)[which(names(df) == "precipitation__sum")] <- "precipitation"
colnames(df)[which(names(df) == "Temp")] <- "Water_temperature"

write.table(df, file.path(dataPath, "df_organized_data_data90_rarefied.csv"), sep=";", dec=",", col.names=NA)

#######################################################
####  correlation plot (pearson correlation) ##########
#######################################################

library(corrplot)
p.mat <- cor.mtest(df)$p
corr_mat=round(cor(df, use = "pairwise.complete.obs"),2)

corrplot(corr_mat, type = "upper", order = "hclust", mar = c(0,0,0,0), method = "square", outline = T, 
         addgrid.col = "darkgray", addrect = 8, rect.col = "black", rect.lwd = 5, tl.col = "black", 
         tl.cex = 0.75, cl.cex = 0.75,
         p.mat = p.mat, sig.level = 0.05, insig = "blank")

dev.copy(pdf, file.path(dataPath, "env_pearson.eps"), width = 6, height = 6)
dev.off()


###############################################################################
###                correlation between environemntal parameters             ###
###############################################################################

# constrained ordination: to see how environmental variables are associated with these changes in community composition. 
# 1) We constrain the ordination axes to linear combinations of environmental variables. 
# 2) We then plot the environmental scores onto the ordination

########################### CAP ##########################################
# see: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

# Remove data points with missing metadata
# Add the categorical variables
# Remove water temp and water level (because large data missing for H42 and H43)
# Remove "air_temperature__mean", because it has nothing to do with the samples
# remove Water_temperature, since it didn't affect the community compostion but got many missing values
# remove KS4.3,KB8.2 
# remove shannon, sodium

u3 <- cbind(u2[,c(1:13,39,40)],df[,-c(8, 9, 19, 21:22,24,25)])

df1 <-na.omit(u3) # 101 reduced to 84 samples (with water temp), 97 (without water temp)
write.table(df1, file.path(dataPath, "env_data_NA_omit.csv"), sep=";", dec=",", col.names = NA)
my_data_prop_sqrt_not_na <- data_rarefy

sample_data(my_data_prop_sqrt_not_na)<-df1

my_data_prop_sqrt_not_na <- prune_taxa(taxa_sums(my_data_prop_sqrt_not_na) > 0, my_data_prop_sqrt_not_na)

bray_not_na <- phyloseq::distance(physeq = my_data_prop_sqrt_not_na, method = "bray")
# 117 samples after removing NA (original 142)

############### 1) Environmental variables #########################
##### CAP for quantitative data #####
# Run DCA (Detrended correspondence analysis) to determine whether to use RDA or CCA
dca <- ordinate(physeq = my_data_prop_sqrt_not_na, 
                method = "DCA")
dca # or summary(dca) to check the "Axis lengths" of the DCA1 (Important) # 5.5844 or 4.0760 (after removing NA)
# If this value < 3, it is better to use RDA
# If this value > 4, it is better to use CCA
# If this value is between 3 and 4, either use CCA or RDA

# CAP ordinate (NOTE: CCA is really slow...)
# install.packages("vegan")

####################################################################################

##  Envfit + CCA/RDA: much more powerful when more complex factors are tested  ##

####################################################################################

library(vegan)
otutable <- vegan_otu(my_data_prop_sqrt_not_na) # import otu table
variable.names(df1)
dat <- data.frame(scale(sample_data(my_data_prop_sqrt_not_na)[,-c(1:15)])) # remove categorical data and Shannon

################# variance inflation factors ################### 
# Linear dependencies between constraints can be investigated via the variance inflation factor or VIF
# VIF is a measure of how much the variance of $\hat{\beta}_j$ is inflated by presence of other covariates
# Lots of rules of thumb
# VIF >= 20 indicates strong collinearity in constraints
# VIF >= 10 potnetially of concern & should be looked at
# they will be completely removed from the estimation, and no biplot scores or centroids are calculated for these aliased constraints. 

ord_all <- capscale(bray_not_na~., dat)
temp <- vif.cca(ord_all)
temp
select_para <- names(temp[temp < 10])
select_para

##### 2. RDA #####
dat1 <- dat[,select_para]
sample_data(my_data_prop_sqrt_not_na)<- dat1

bray_not_na <- phyloseq::distance(physeq = my_data_prop_sqrt_not_na, method = "bray")

cap_ord <- ordinate(
  physeq = my_data_prop_sqrt_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~.)

# Significance test
set.seed(123)
anova(cap_ord, by="margin", perm=999) # marginal effects of the terms 
drop1(cap_ord, test="perm")

#### Final model
select_para1 <- select_para[c(2, 3,4, 6, 7, 8, 9)]

dat2 <- dat[,select_para1]
sample_data(my_data_prop_sqrt_not_na)<- dat2

bray_not_na <- phyloseq::distance(physeq = my_data_prop_sqrt_not_na, method = "bray")

cap_ord <- ordinate(
  physeq = my_data_prop_sqrt_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ .)

# Significance test
anova(cap_ord, perm=999)
set.seed(123)
anova(cap_ord, by="margin", perm=999) # marginal effects of the terms 
anova(cap_ord, by="terms") # sequential
drop1(cap_ord, test="perm")
ordistep(cap_ord, perm.max = 200)
summary(cap_ord)
RsquareAdj(cap_ord) # ## RsquareAdj gives the same result as component [a] of varpart

################### Biplot #####################
# CAP plot
sample_data(my_data_prop_sqrt_not_na) <- u # add categorical data for plotting

cap_plot_env <- plot_ordination(
  physeq = my_data_prop_sqrt_not_na, 
  ordination = cap_ord, 
  color="Well",
  axes = c(1,2))+  theme_bw() +
  geom_point(aes(colour = Well), alpha = 0.8, size = 4) +
  scale_color_manual(values = set_colors_well) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14))
   
cap_plot_env
cap_plot_env$layers<-cap_plot_env$layer [-1]

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1*1, 
                 yend = CAP2*1, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.1 * CAP1, 
                 y = 1.1 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
p <- cap_plot_env + 
  geom_segment(
    mapping = arrow_map, 
    size = 0.75, 
    data = arrowdf, 
    color = "blue", 
    arrow = arrowhead
  ) + 
  geom_text(mapping = label_map, size = 5, color="blue", data = arrowdf, show.legend = FALSE) +
  ggtitle("CAP") +
  scale_y_continuous(breaks = c(-1.5, -1, -0.5,0,  0.5, 1, 1.5))+
  theme(legend.key = element_blank(),  # removes the box around each legend item
        legend.position = "right") +   # legend at the bottom
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=14), legend.text=element_text(size=14)) 
p

ggsave(file.path(dataPath, "cap_linking_env.pdf"), width = 7, height = 5)


#################################################################################################
########################  PCA: hydrochemical difference between samples  ########################
#################################################################################################

set.seed(123)
metadata_tb <- dat
variable.names(metadata_tb)
metadata_tb.pca <- prcomp(dat, center = TRUE, scale. = TRUE)
# center = TRUE: a logical value indicating whether the variables should be shifted to be zero centered
# center = TRUE, scale. = TRUE: the variables is adviced to be scaled to have unit variance before the analysis takes place

# To check the results
# The print method returns the standard deviation of each of the four PCs, and their rotation (or loadings), 
# which are the coefficients of the linear combinations of the continuous variables.
print(metadata_tb.pca)

# The summary method describe the importance of the PCs. 
# We can see there that the first two PCs accounts for more than 51% of the variance of the data.
summary(metadata_tb.pca)

metadata_tb.pca$rotation # the matrix of variable loadings 

# plot method: The plot method returns a plot of the variances (y-axis) associated with the PCs (x-axis). 
# The Figure below is useful to decide how many PCs to retain for further analysis. 
# we can see that the first two PCs explain most of the variability in the data.
plot(metadata_tb.pca, type = "l")
plot(metadata_tb.pca)

#### simple PCA plot
biplot(metadata_tb.pca, choices = 1:2, scale = 1)

## Better PCA plot
metadata_tb$Cluster <- sample_data(my_data_prop_sqrt_not_na)$bCluster
row.names(metadata_tb)
row.names(sample_data(my_data_prop_sqrt_not_na))
library(ggbiplot)
ggbiplot(metadata_tb.pca, 
         var.axes = TRUE, # Draw arrows for the variables?
         obs.scale = 1, var.scale = 1, 
  color=metadata_tb$Cluster,
  groups=metadata_tb$Cluster) +
  scale_colour_manual(name="Cluster", values= set_colors_cluster) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(dataPath, "PCA.png"), width = 8, height = 6)


##### calculate mean and standard deviation of samples
dx <- u2[,c(5, 14:37)]
dx <- na.omit(dx)
data_mean <- ddply(dx, .(Well), colwise(mean))
data_sd <- ddply(dx, .(Well), colwise(sd))
data_mean
data_sd
write.table(data_mean, file.path(dataPath, "mean_env_data.csv"), sep=";", dec=",", col.names=NA)
write.table(data_sd, file.path(dataPath, "sd_env_data.csv"), sep=";", dec=",", col.names=NA)

###############################
### Variance partitioning #####
###############################

#http://www.hiercourse.com/docs/variationPartitioning.html

datx <- df1
well.pcnm <- pcnm(dist(datx[,c("LON", "LAT")]))
otutable <- vegan_otu(my_data_prop_sqrt_not_na) # import otu table
mod <- varpart(vegdist(otutable), ~ .,  well.pcnm$vectors, data=dat[, select_para[1:9]])
mod
showvarparts(2, bg = c("hotpink", "skyblue"))
plot(mod, bg = c("hotpink", "skyblue"))
dev.off()

## Add sampling campaign (used)
mod1 <- varpart(vegdist(otutable), ~ ., datx$Samp_Cam, well.pcnm$vectors, data=dat[, select_para[1:9]])
mod1

showvarparts(3, bg = 2:4)
plot(mod1, bg = 2:4)
dev.copy(pdf, file.path(dataPath, "variance_partition_3factors.pdf"), width = 6, height = 6)
dev.off()

###########################################################################################################################################
############## Optional: Subset the data to the hydrochemically connected wells in HTL/HTU, and rerun the code ############################
### HTL ###
HTL_abun <- subset_samples(data_rarefy, Well=="H13"|Well=="H31"|Well=="H41"|Well=="H51")
HTL_abun <- filter_taxa(HTL_abun, function(x) sum(x) > 0, TRUE)
my_data_prop_sqrt_not_na <- HTL_abun

### HTU ###
HTU_abun <- subset_samples(data_rarefy, Well=="H42"| Well=="H43"|Well=="H53" )
HTU_abun <- filter_taxa(HTU_abun, function(x) sum(x) > 0, TRUE)
my_data_prop_sqrt_not_na <- HTU_abun

