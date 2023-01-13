library(ggplot2)
library(reshape)
rm(list = ls())
setwd("/Users/fangshu/Desktop/research/design/simulation/result/plot/")
opendir <- '/Users/fangshu/Desktop/research/design/simulation/result/table/added_effect/'
trt <- c("NE","Trt1","Trt2","Trt3")
ylim <- c(0,1)
angle1 <- c(45,45,135) 
angle2 <- c(45,135,135) 
density <- c(0,13,10)
col <- 1 # rainbow(7)
ylab <- c("Type I Error for Treatment Effect","Power for Treatment Effect","Power for Treatment Effect","Power for Treatment Effect")
for(i in 1:4){
  dat <- read.csv(paste0(opendir,paste0(trt[i],"_sig_table.csv")),header = T)
  pp <- dat[,c(2:4)]
  pp <- as.matrix(pp)
  pp <- t(pp)
  colnames(pp) <- paste0("W ",c(1:5))
  
  png(paste0("sim_added_effect_",i,".png"))
  op <- par(mar=c(3,5,1,1))
  #barplot(pp, beside=TRUE, ylim=ylim, col=col, angle=angle1, density=density,axis.lty=1,font.axis = 2.2,ylab=ylab[i])
  barplot(pp, beside=TRUE, ylim=ylim, col=col, angle=angle1, density=density,axis.lty=1,font.axis = 2.2,ylab=ylab[i],
          cex.lab=1.6,cex.axis=1.6,cex.names = 1.6)
  barplot(pp, add=TRUE, beside=TRUE, ylim=ylim, col=col, angle=angle2, density=density,xaxt="none",ylab=ylab[i],
          cex.lab=1.6,cex.axis=1.6)
  if(i==2){
    legend("topright", legend=c('M1: no random effect','M2: only pen effect','M3: pen and cohort effect'), ncol=1, fill=TRUE, col=col, 
           angle=angle1, density=density,cex = 1.5)
    par(bg="transparent")
    legend("topright", legend=c('M1: no random effect','M2: only pen effect','M3: pen and cohort effect'), ncol=1, fill=TRUE, col=col, 
           angle=angle2, density=density,cex = 1.5)
  }
  par(op)
  dev.off()
}

opendir <- '/Users/fangshu/Desktop/research/design/simulation/result/table/binary/'
trt <- c("NE","TMS","ENR","TIL")
for(i in 1:4){
  dat <- read.csv(paste0(opendir,paste0(trt[i],"_sig_table.csv")),header = T)
  pp <- dat[,c(2:4)]
  pp <- as.matrix(pp)
  pp <- t(pp)
  colnames(pp) <- paste0("W ",c(1:5))
  
  png(paste0("sim_binary_",i,".png"))
  op <- par(mar=c(3,5,1,1))
  barplot(pp, beside=TRUE, ylim=ylim, col=col, angle=angle1, density=density,axis.lty=1,font.axis = 2.2,ylab=ylab[i],
          cex.lab=1.6,cex.axis=1.6,cex.names = 1.6)
  barplot(pp, add=TRUE, beside=TRUE, ylim=ylim, col=col, angle=angle2, density=density,xaxt="none",ylab=ylab[i],
          cex.lab=1.6,cex.axis=1.6)
  if(i==2){
    legend("topright", legend=c('M1: no random effect','M2: only pen effect','M3: pen and cohort effect'), ncol=1, fill=TRUE, col=col, 
           angle=angle1, density=density,cex = 1.5)
    par(bg="transparent")
    legend("topright", legend=c('M1: no random effect','M2: only pen effect','M3: pen and cohort effect'), ncol=1, fill=TRUE, col=col, 
           angle=angle2, density=density,cex = 1.5)
  }
  par(op)
  dev.off()
}

library(ggplot2)
library(reshape)
rm(list = ls())
setwd("/Users/fangshu/Desktop/research/design/simulation/result/plot/additional/")
opendir <- '/Users/fangshu/Desktop/research/design/simulation/result/table/added_equal_rf/'
trt <- c("NE","Trt1","Trt2","Trt3")
ylim <- c(0,1)
angle1 <- c(45,45,135) 
angle2 <- c(45,135,135) 
density <- c(0,13,10)
col <- 1 # rainbow(7)
ylab <- c("Type I Error for Treatment Effect","Power for Treatment Effect","Power for Treatment Effect","Power for Treatment Effect")
for(i in 1:4){
  dat <- read.csv(paste0(opendir,paste0(trt[i],"_sig_table.csv")),header = T)
  pp <- dat[,c(2:4)]
  pp <- as.matrix(pp)
  pp <- t(pp)
  colnames(pp) <- paste0("W ",c(1:5))
  
  png(paste0("sim_added_equal_",i,".png"))
  op <- par(mar=c(3,5,1,1))
  barplot(pp, beside=TRUE, ylim=ylim, col=col, angle=angle1, density=density,axis.lty=1,font.axis = 2.2,ylab=ylab[i],
          cex.lab=1.6,cex.axis=1.6,cex.names = 1.6)
  barplot(pp, add=TRUE, beside=TRUE, ylim=ylim, col=col, angle=angle2, density=density,xaxt="none",ylab=ylab[i],
          cex.lab=1.6,cex.axis=1.6)
  if(i==2){
    legend("topright", legend=c('M1: no random effect','M2: only pen effect','M3: pen and cohort effect'), ncol=1, fill=TRUE, col=col, 
           angle=angle1, density=density,cex = 1.5)
    par(bg="transparent")
    legend("topright", legend=c('M1: no random effect','M2: only pen effect','M3: pen and cohort effect'), ncol=1, fill=TRUE, col=col, 
           angle=angle2, density=density,cex = 1.5)
  }
  par(op)
  dev.off()
}


