library(vegan)   # vegdist
library(foreach)

##### Turnover of OTUs #### slightly modified (distance was changed to vegdist)

main_raup_crick_abundance = function(spXsite, null_bray_curtis, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE) {
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model

  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }

  ## count number of sites
  n_sites<-nrow(spXsite)

  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
    
  ##looping over each pairwise community combination:
  
  rcindices <- which(lower.tri(results))
  for (j in seq_along(rcindices)) {
  #for(null.one in 1:(nrow(spXsite)-1)){
  #  for(null.two in (null.one+1):nrow(spXsite)){
      i <- rcindices[j]
      null.two <- ifelse(i %% n_sites == 0, n_sites, i %% n_sites)
      null.one <- ifelse(null.two == n_sites, (i %/% n_sites), (i %/% n_sites) + 1)

      ## empirically observed bray curtis
      obs.bray = vegdist(spXsite[c(null.one,null.two),],'bray');
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis[,j] == obs.bray);
      
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis[,j] < obs.bray);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
    #}; ## end null.two loop
    
  }; ## end null.one loop
  
  if (as.distance.matrix) { ## return as distance matrix if so desired
    results <- as.dist(results)
  }	
  
  return(results)
}

####
# Start of processing
####
#rc_bray <- pre_raup_crick_abundance(otutable, plot_names_in_col1=FALSE, reps=1000)

# Set up the paths to the files
basePath <- "path/to/assembly/data"
rawPath <- file.path(basePath, "pnk75")
calcPath <- file.path(rawPath, "RC")
csvFile <- file.path(rawPath, "pnk75_otutable.csv")

otutable <- read.csv(csvFile, sep=";", dec = ",", row.names=1, header=TRUE)
rcfiles <- dir(calcPath, pattern="idiv_23277.*", full.names=T)

pre_rc_bray <- foreach(r = rcfiles, .combine='rbind', .multicombine=TRUE) %do% {
   readRDS(r)
}

dim(otutable)
dim(pre_rc_bray)

rc_bray <- main_raup_crick_abundance(otutable, pre_rc_bray, plot_names_in_col1=FALSE, reps=nrow(pre_rc_bray))
write.csv(as.matrix(rc_bray), file.path(calcPath, "rc_bray.csv"), quote=F)

