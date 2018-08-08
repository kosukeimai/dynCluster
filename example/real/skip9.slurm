#!/bin/bash
#SBATCH --nodes=1  # node count  
#SBATCH --ntasks-per-node=8
#SBATCH -t 71:55:00     
#SBATCH -J skip9
#SBATCH -o log.%j
#SBATCH --mail-type=begin 
#SBATCH --mail-type=end 
#SBATCH --mail-user=hubertj@princeton.edu 

module load boost/1.53.0
./../src/ZTM ./config.txt 1

