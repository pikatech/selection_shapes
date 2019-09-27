#!/bin/bash
#SBATCH --mem-per-cpu=26G
srun -o idiv_%j_%t Rscript bNTI_randomized.R 10xHT/HTU/HTU_otutable.csv 10xHT/HTU/HTU.tre 100
