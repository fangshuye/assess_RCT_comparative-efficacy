library(ggplot2)
library(ggtextures)
library(reshape)
rm(list = ls())
setwd("/Users/fangshu/Desktop/research/design/simulation/result/plot/")
images = c(
  m1 = "white.png",
  m2 = "cross.png",
  m3 = "dia.png")

opendir <- '/Users/fangshu/Desktop/research/design/simulation/result/table/added_effect/'
filenames <- list.files(opendir,pattern = "_sig_table.csv")
dat_all <- c()
for(i in 1:4){
  dat <- read.csv(paste0(opendir,filenames[i]),header = T)
  fullname <- unlist(strsplit(filenames[i], "_"))
  dat$treatment <- rep(fullname[1],nrow(dat))
  dat_all <- rbind(dat_all,dat)
}

dat_p <- melt(dat_all, id.vars=c('X',"treatment"))
dat_p$treatment <- factor(dat_p$treatment,levels = c("NE","Trt1","Trt2","Trt3"))
hum.names <- as_labeller(c('NE' ="No effect (mean difference = 0)",
                           'Trt1' = "Trt1 vs Trt2 (mean difference = 2.050)",
                           'Trt2' = "Trt1 vs Trt3 (mean difference = 7.778)",
                           'Trt3' = "Trt1 vs Trt4 (mean difference = 22.507)"))
p <- ggplot(dat_p,aes(x = X, y= value, image = variable))+
  facet_wrap(.~treatment,nrow=2,labeller=hum.names)+
  geom_textured_bar(stat="identity", width=.5,position = "dodge")+
  ylab("Type I Error or Power of treatment effect simulation in continuous situation")+
  xlab("")+
  scale_image_manual(values=images,
                     name = "Model Type",
                     labels = c("model with no random effect",
                                "model with only pen effect",
                                "model with all random effects"))+theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5))+
  ggtitle('')
ggsave("/Users/fangshu/Desktop/research/design/simulation/result/plot/sim_added_effect.png",p,dpi = 300,
       width =7, height = 7)



opendir <- '/Users/fangshu/Desktop/research/design/simulation/result/table/binary/'
filenames <- list.files(opendir,pattern = "_sig_table.csv")
dat_all <- c()
for(i in 1:4){
  dat <- read.csv(paste0(opendir,filenames[i]),header = T)
  fullname <- unlist(strsplit(filenames[i], "_"))
  dat$treatment <- rep(fullname[1],nrow(dat))
  dat_all <- rbind(dat_all,dat)
}
dat_p <- melt(dat_all, id.vars=c('X',"treatment"))
dat_p$treatment <- factor(dat_p$treatment,levels = c("NE","TMS","ENR","TIL"))
hum.names <- as_labeller(c('NE' ="No effect (odd ratio = 1)",
                           'TMS' = "PLC vs TMS (odd ratio = 0.7875)",
                           'ENR' = "PLC vs ENR (odd ratio = 0.5425)",
                           'TIL' = "PLC vs TIL (odd ratio = 0.3099)"))


p <- ggplot(dat_p,aes(x = X, y= value, image = variable))+
  facet_wrap(.~treatment,nrow=2,labeller=hum.names)+
  geom_textured_bar(stat="identity", width=.5,position = "dodge")+
  ylab("Type I Error or Power of treatment effect simulation")+
  xlab("")+
  scale_image_manual(values=images,
                     name = "Model Type",
                     labels = c("model with no random effect",
                                "model with only pen effect",
                                "model with all random effects"))+theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0.5))+
  ggtitle('')

ggsave("/Users/fangshu/Desktop/research/design/simulation/result/plot/sim_binary.png",p,dpi = 300,
       width =7, height = 7)
