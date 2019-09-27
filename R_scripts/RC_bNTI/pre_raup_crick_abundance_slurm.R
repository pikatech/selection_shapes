path <- "path/to/assembly/data" 
setwd(path)
source("pre_raup_crick_abundance.R")
f <- commandArgs(T)
if (length(f) == 2) {
   csvFile <- f[1]
   reps <- as.integer(f[2])
   otutable <- read.csv(csvFile, sep=";", dec = ",", row.names=1, header=TRUE)
   rc_bray <- pre_raup_crick_abundance(otutable, plot_names_in_col1=FALSE, reps=reps)
   saveRDS(rc_bray, file=stdout(), ascii = TRUE)
} else {
   stop("Please specify a csv file and the number of repetitions!")
}
