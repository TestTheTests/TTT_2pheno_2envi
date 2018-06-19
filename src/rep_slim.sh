#!/bin/bash
set -e
set -u
set -o pipefail

######################################################################################
#
# File	  : rep_slim_script.sh
# History : 20180515  Created by Brett Ford (BF)
#
#######################################################################################
#
# This script runs simulations for all unique combinations of parameters
# provided
#
#######################################################################################

# set path to where simulation file is
mypath="/Users/brettford/Desktop/Northeastern/slim/TTT_2pheno_2envi/simfiles"
cd $mypath

#separate values in arrays by spaces
declare -a lista=(0.25 0.5 1 1.25 1.5)
declare -a listb=(-0.75 -0.5 -0.25 0.25 0.5 0.75)

#rep the for loop for as many constants that you intend on defining 
#at the command line
for i in "${lista[@]}" 
do
  for j in "${listb[@]}"
  do
	 slim -d "sigma_K=${i}" -d "QTL_cov=${j}" LocalAdapt2trait_1mut_2env_20180514_d0.009_c0.027_var.slim &
         sleep 10s
  done
done
