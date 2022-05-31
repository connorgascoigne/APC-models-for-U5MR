#!/bin/bash    

for interaction in {1..4}; 
do 
  
  echo "$interaction"
  Rscript --vanilla post_fit_plots.R $interaction
  
done
