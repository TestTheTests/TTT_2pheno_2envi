######################################
#
# effects_proc_script.R
# Brett Ford
# Created 20180525
#
# This script takes the outputm2effects text files from slim
# and plots a scatterplot of effect sizes for phenotype 0 against
# effect sizes for phenotype1 with marginal histograms 
# for the effect sizes of each phenotype
#
# Usage: Rscript effects_proc_script.R
#
######################################

# Set paths to directories
# TIP: If you git clone this repository, you should just have to change the root path
root_path <- ("/Users/brettford/Desktop/Northeastern/slim/TTT_2pheno_2envi")
results_path <- paste0(root_path, "/results/")
figures_path <- paste0(root_path, "/figures/")

# Set working directory
setwd(results_path)

##Will have to specify directories relative to TTT_2pheno_2envi when I complete the script

# Load libraries:
library(Hmisc) # Package required to calculate Pearson correlation with p-val
library(SDMTools)
##Scatterplot-Histogram Function:
scatterhist <- function(x, y, z, mat, xlab = "", ylab = "", plottitle="",
                        xsize=1){
  #Zones creates a matrix that describes where to plot the figures below, in sequence,
  #for example, the code below plots the first figure in the entire top row (because there
  #are 3 columns)
  effects_table <- data.frame(e0=x, e1=y, f=z)
  fixed_loci <- effects_table[effects_table$f==1,]
  fixed_corr <- rcorr(fixed_loci$e0, fixed_loci$e1, type="pearson")
  fixed_corr <- round(fixed_corr$r[2],2)
  fixed_corr <- paste("r(fixed)=", fixed_corr)

  zones <- matrix(c(1,1,1,
                    8,5,7,
                    2,6,4,
                    0,3,0), ncol = 3, byrow = TRUE)
  layout(zones, widths=c(1,4,1), heights = c(1,3,10,1))

  # tuning to plot histograms nicely
  all_effects <- c(x,y)
  all_effects_min <- floor(range(all_effects)[1])
  all_effects_max <- ceiling(range(all_effects)[2])
  xhist = hist(x, breaks=c(seq(floor(range(x)[1]),
                               ceiling(range(x)[2]), 0.25)), plot=FALSE)

  yhist = hist(y, breaks=c(seq(floor(range(y)[1]),
                               ceiling(range(y)[2]), 0.25)),
               plot=FALSE)
  top <- max(c(xhist$counts, yhist$counts))

  # for all three titles:
  #   drop the axis titles and omit boxes, set up margins
  par(xaxt="n", yaxt="n", bty="n", mar = c(.3,.2,.3,0) +.05)
  # fig 1 from the layout
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0,paste(plottitle), cex=1.5)
  par(xaxt="n", yaxt="n", bty="n", mar = c(.3,.5,.3,0) +.05)
  # fig 2
  plot(x=1,y=1,type="n", ylim=c(-1, 1), xlim=c(- 1, 1))
  text(0.5,0,paste(ylab), cex=1.5, srt=90) #srt rotates the text 90 degrees
  # fig 3
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0,paste(xlab), cex=1.5)

  # fig 4, the first histogram, needs different margins
  # no margin on the left
  par(mar = c(2,0,1,1))
  #par(mar = c(2,2,2,2))
  barplot(yhist$counts, axes = FALSE, xlim = c(0, top),
          space = 0, horiz = TRUE)
  # fig 5, other histogram needs no margin on the bottom
  par(mar = c(0,2,1,1))
  barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0)
  # fig 6, finally, the scatterplot-- needs regular axes, different margins
  par(mar = c(2,2,.5,.5), xaxt="s", yaxt="s", bty="o") #n means to supress, otherwise s means plot
  # this color allows traparency & overplotting-- useful if a lot of points
  xsize=1
  x_floor <- floor(range(x)[1])
  x_ceiling <- ceiling(range(x)[2])
  y_floor <- floor(range(y)[1])
  y_ceiling <- ceiling(range(y)[2])
  x_tick <- (x_ceiling-x_floor)/0.5
  y_tick <- (y_ceiling-y_floor)/0.5
  effect_corr <- rcorr(x, y, type="pearson")
  effect_corr <- round(effect_corr$r[2],2)
  effect_corr <- paste("r=", effect_corr)
  rbPal <- colorRampPalette(c('red','blue'))
  pt_colors <- rbPal(10)[as.numeric(cut(z,breaks = 10))]
  plot(x, y , xlim= c(x_floor, x_ceiling), ylim=c(y_floor, y_ceiling), xaxp= c(x_floor, x_ceiling, x_tick), yaxp= c(y_floor, y_ceiling, y_tick), pch=19, col=pt_colors, cex=1.5)
  text((x_floor),(y_ceiling-0.25),paste0(effect_corr, "\n", fixed_corr), cex=1.25, pos=4)
                                                                                                                      

  par(mar = c(0.5,0.5,0.5,0.5))
  plot(x=1,y=1,type="n",ylim=c(0,1), xlim=c(0,1), axes = FALSE)
  pnts <- cbind(x =c(0.4,0.6,0.6,0.4), y =c(0.8,0.8,0.2,0.2))
  legend.gradient(pnts, cols = rbPal(100), limits = c(0, 1),
                  title = "Fixed?", cex=1.5)

  par(mar = c(0.5,0.5,0.5,0.5))
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1), axes = FALSE)
  text(0,0,mat, cex=1.25)
}


# Get all filenames for effect size files in results path
effect_size_filenames <- list.files(path=results_path, pattern=".+_outputm2effects.+")

#Plot for every effects file in results_path
for (i in effect_size_filenames) {
  #Read in effects table for simulation
  effects_table <- read.table(i, header=T, sep="\t")

  #Simulation details
  simulation_filename_contents <- unlist(strsplit(i, split="[_]"))
  simulation <- simulation_filename_contents[1]

  #Set var-cov matrix
  var <- grep('QTL_var',readLines(paste0(simulation, "_simulation_parameters.txt")), value=TRUE)
  var <- unlist(strsplit(var, split="\t"))
  covar <- grep('QTL_cov',readLines(paste0(simulation, "_simulation_parameters.txt")), value=TRUE)
  covar <- unlist(strsplit(covar, split="\t"))
  var <- var[2]
  covar <- covar[2]

  qtl_matrix <- matrix(c(var,covar,covar,
         var), nrow= 2, ncol = 2)
  qtl_matrix <- paste0(qtl_matrix[1], " ", qtl_matrix[2], "\n", qtl_matrix[3], " ", qtl_matrix[4])
  #Set image details before plotting
  png(paste0(figures_path,simulation, "_dist_effect_sizes.png"), width=6, height = 6, units="in", res=500)
  scatterhist(effects_table$e0, effects_table$e1, effects_table$f, qtl_matrix, xlab="Phenotype0 Effect Size", ylab="Phenotype1 Effect Size", 
              plottitle=paste(simulation, "Distribution of Effect Sizes"))

  # write png file
  dev.off()
}

