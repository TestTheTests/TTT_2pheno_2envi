######################################
#
# chrom_effects_viz.R
# Brett Ford
# Created 20180618
#
# This script takes the outputMuts2 files in the results folder
# and plots the effect sizes for each phenotype to visualize
# where mutations are arising, at what effect size, and 
# at what position relative to other mutations (i.e., modular?)
#
#
# Usage: Rscript chrom_effects_viz.R
#
######################################

# Set paths to directories
# TIP: If you git clone this repository, you should just have to change the root path
root_path <- ("/Users/brettford/Desktop/Northeastern/slim/TTT_2pheno_2envi")
results_path <- paste0(root_path, "/results/")
figures_path <- paste0(root_path, "/figures/")

#Obtain effect files
effects_filenames <- list.files(path=results_path, pattern=".+outputm2effects_")

  


for(s in 1:length(effects_filenames)) {
  #Get file with positions and plot
  m2_pos <- read.table(effects_filenames[s], header = TRUE, sep = "\t")
  
  #Some early files didn't have position in them, this is to make sure those files aren't analyzed:
  if (is.null(m2_pos$position)==TRUE){
    print("No position in file (early version of m2effects file)")
  }
  else{
  #Get simulation number to write with filename
  file_contents <- unlist(strsplit(effects_filenames[s], split="[_]"))
  sim <- file_contents[1]
  
  #Make new columns with colors to represent positive or negative effect
  m2_pos$e0_color <- "NA"
  m2_pos$e1_color <- "NA"
  
  #Assign colors based on whether value is negative or positive
  for(i in 1:length(m2_pos$e0)) {
    if (m2_pos$e0[i] < 0) {
      m2_pos$e0_color[i]<- "red"
    }
    else {
      m2_pos$e0_color[i]<- "blue"
    }
  }
  for(i in 1:length(m2_pos$e1)) {
    if (m2_pos$e1[i] < 0) {
      m2_pos$e1_color[i] <- "red"
    }
    else {
      m2_pos$e1_color[i] <- "blue"
    }
  }
  #Look at first few columns to make sure color assignment worked
  head(m2_pos)
  
  #Now take the absolute value of phenotypes; leave phenotype1 positive and convert phenotype2 to negative
  m2_pos$e0_abs <- "NA"
  m2_pos$e1_abs <- "NA"
  for(i in 1:length(m2_pos$e0)) {
    m2_pos$e0_abs[i]<- abs(m2_pos$e0[i])
  }
  for(i in 1:length(m2_pos$e1)) {
    m2_pos$e1_abs[i]<- -1*abs(m2_pos$e1[i])
  }
  
  # Make levels for colors
  levels(m2_pos$e0_color) <- c("red", "blue")
  levels(m2_pos$e1_color) <- c("red", "blue")
  
  
  png(paste0(figures_path, sim, "_mut_effect_allchrom.png"), width=8, height = 6, units="in", res=500)
  plot(m2_pos$position, m2_pos$e0_abs, yaxt="n", type="n", main="Whole Genome Molecular Architecture", xlab= "Chromosome Position", ylab= "Effect Size", xlim = c(0,350000), ylim = c(-3,3))
  axis(2, at=-3:3, labels=c(3, 2, 1, 0, 1, 2, 3))
  lines(m2_pos$position, m2_pos$e0_abs, type="h", col=c(m2_pos$e0_color))
  lines(m2_pos$position, m2_pos$e1_abs, type="h", col=c(m2_pos$e1_color))
  abline(0, 0, lwd=1)
  text(300000, 2.5, "Phenotype0")
  text(300000, -2.5, "Phenotype1")
  legend("topleft", legend=c("negative", "positive"), bty="n", col="gray32", pch=22, cex=1, pt.bg=c("red", "blue"), y.intersp=1.0)
  dev.off()
  
  #Visualize mutational effects for every other chromosome, where m2 mutations arise
  chr1 <- m2_pos[m2_pos$position<50000,]
  chr3 <- m2_pos[m2_pos$position<150000,]
  chr3 <- m2_pos[m2_pos$position>50000,]
  chr5 <- m2_pos[m2_pos$position<250000,]
  chr5 <- m2_pos[m2_pos$position>200000,]
  chr7 <- m2_pos[m2_pos$position<350000,]
  chr7 <- m2_pos[m2_pos$position>300000,]
  
  #Chromosome1
  png(paste0(figures_path, sim, "_mut_effect_chrom1.png"), width=8, height = 6, units="in", res=500)
  plot(chr1$position, chr1$e0_abs, yaxt="n", type="n", main="Chromosome 1 Molecular Architecture", xlab= "Chromosome Position", ylab= "Effect Size", xlim = c(0,50000), ylim = c(-3,3))
  axis(2, at=-3:3, labels=c(3, 2, 1, 0, 1, 2, 3))
  lines(chr1$position, chr1$e0_abs, type="h", col=c(chr1$e0_color), lwd=4)
  lines(chr1$position, chr1$e1_abs, type="h", col=c(chr1$e1_color), lwd=4)
  abline(0, 0, lwd=4)
  text(40000, 2.5, "Phenotype0")
  text(40000, -2.5, "Phenotype1")
  legend("topleft", legend=c("negative", "positive"), bty="n", col="gray32", pch=22, cex=1, pt.bg=c("red", "blue"), y.intersp=1.0)
  dev.off()
  
  #Chromosome3
  png(paste0(figures_path, sim, "_mut_effect_chrom3.png"), width=8, height = 6, units="in", res=500)
  plot(chr3$position, chr3$e0_abs, yaxt="n", type="n", main="Chromosome 3 Molecular Architecture", xlab= "Chromosome Position", ylab= "Effect Size", xlim = c(100000,150000), ylim = c(-3,3))
  axis(2, at=-3:3, labels=c(3, 2, 1, 0, 1, 2, 3))
  lines(chr3$position, chr3$e0_abs, type="h", col=c(chr3$e0_color), lwd=4)
  lines(chr3$position, chr3$e1_abs, type="h", col=c(chr3$e1_color), lwd=4)
  abline(0, 0, lwd=4)
  text(140000, 2.5, "Phenotype0")
  text(140000, -2.5, "Phenotype1")
  legend("topleft", legend=c("negative", "positive"), bty="n", col="gray32", pch=22, cex=1, pt.bg=c("red", "blue"), y.intersp=1.0)
  dev.off()
  
  #Chromosome5
  png(paste0(figures_path, sim, "_mut_effect_chrom5.png"), width=8, height = 6, units="in", res=500)
  plot(chr5$position, chr5$e0_abs, yaxt="n", type="n", main="Chromosome 5 Molecular Architecture", xlab= "Chromosome Position", ylab= "Effect Size", xlim = c(200000,250000), ylim = c(-3,3))
  axis(2, at=-3:3, labels=c(3, 2, 1, 0, 1, 2, 3))
  lines(chr5$position, chr5$e0_abs, type="h", col=c(chr5$e0_color), lwd=4)
  lines(chr5$position, chr5$e1_abs, type="h", col=c(chr5$e1_color), lwd=4)
  abline(0, 0, lwd=4)
  text(240000, 2.5, "Phenotype0")
  text(240000, -2.5, "Phenotype1")
  legend("topleft", legend=c("negative", "positive"), bty="n", col="gray32", pch=22, cex=1, pt.bg=c("red", "blue"), y.intersp=1.0)
  dev.off()
  
  #Chromosome7
  png(paste0(figures_path, sim, "_mut_effect_chrom7.png"), width=8, height = 6, units="in", res=500)
  plot(chr7$position, chr7$e0_abs, yaxt="n", type="n", main="Chromosome 7 Molecular Architecture", xlab= "Chromosome Position", ylab= "Effect Size", xlim = c(300000,350000), ylim = c(-3,3))
  axis(2, at=-3:3, labels=c(3, 2, 1, 0, 1, 2, 3))
  lines(chr7$position, chr7$e0_abs, type="h", col=c(chr7$e0_color), lwd=4)
  lines(chr7$position, chr7$e1_abs, type="h", col=c(chr7$e1_color), lwd=4)
  abline(0, 0, lwd=4)
  text(340000, 2.5, "Phenotype0")
  text(340000, -2.5, "Phenotype1")
  legend("topleft", legend=c("negative", "positive"), bty="n", col="gray32", pch=22, cex=1, pt.bg=c("red", "blue"), y.intersp=1.0)
  dev.off()
  }
}
