#!/bin/bash
#SBATCH --mem-per-cpu=1G
srun -o idiv_%j_%t Rscript pre_raup_crick_abundance_slurm.R pnk75/pnk75_otutable.csv 20
