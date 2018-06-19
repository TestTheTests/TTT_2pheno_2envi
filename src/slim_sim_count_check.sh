/#!/bin/bash
set -e
set -u
set -o pipefail

######################################################################################
#
# File	  : slim_sim_count_check.sh
# History : 20180515  Created by Brett Ford (BF)
#
#######################################################################################
#
# This script prints the number of parameters from the lists provided
#
#######################################################################################

declare -a lista=(0.25 0.5 1 1.25 1.5)
declare -a listb=(-0.75 -0.5 -0.25 0.25 0.5 0.75)

count=0
for i in "${lista[@]}" 
do
  for j in "${listb[@]}"
  do
     echo $i $j
     count=$(( $count + 1 ))
  done
done

#print total number of simulations
echo "Total number of simulations:" $count
