#####################SIMULATION: 10000 FOR PLC vs TIL#########################
rm(list = ls())
library(lme4)
library(dplyr)
library(lmerTest)
library(doParallel)
set.seed(20200609, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)

### data size
# need to decide
pen_num <- 20
cattle_per_pen <- 40
# fix 
cattle_per_cohort <- 40 
total_cattle <- pen_num*cattle_per_pen
cohort_num <- total_cattle/cattle_per_cohort

# get fixed effect and random effect
var_cohort <- 57.73
var_pen <- 101.76
var_subject <- 5901.47

# PLC
PLC <- 903.274
Trt1 <- 903.274-2.050
Trt2 <- 903.274-7.778
Trt3 <- 903.274-22.507

# NE: no effect
NE <- PLC

# pull together
treatments <- c("Trt1","Trt2","Trt3","NE")
treatment_values <- c(Trt1,Trt2,Trt3,NE)
# choose trt
treatment <- treatments[3]
treatment_value <- treatment_values[3]

nrep <- 10000
dat_sim <- data.frame(subject=c(1:total_cattle),cohort=rep(1:cohort_num,each=cattle_per_cohort))

table_from_sim <- foreach (rep=1:nrep, .combine = rbind)%dopar%{
  
  ################### part 1: simulate the design data #######################
  # design 1
  sim_1 <- dat_sim
  for(c in 1:cohort_num){
    n1 <- (c-1)*cattle_per_cohort
    n2 <- c*cattle_per_cohort
    sim_1[(n1+1):n2,"trt"] <- sample(rep(1:2,each=cattle_per_cohort/2),replace = F)
    sim_1[(n1+1):n2,"pen"] <- sample(rep(1:pen_num,each=cattle_per_cohort/pen_num),replace = F)
  }
  
  # design 2
  sim_2 <- dat_sim
  for(c in 1:cohort_num){
    n1 <- (c-1)*cattle_per_cohort
    n2 <- c*cattle_per_cohort
    sim_2[(n1+1):n2,"trt"] <- sample(rep(1:2,each=cattle_per_cohort/2),replace = F)
  }
  sim_2$pen <- rep(1:pen_num,each=cattle_per_pen) # animals are systematically allocated to pens-sequentially
  
  # design 3
  sim_3 <- dat_sim
  sim_3$pen <- rep(1:pen_num,each=cattle_per_pen) # equal to cohort in this design
  pen_trt <- sample(rep(1:2,each=pen_num/2),replace = F)
  sim_3$trt <- rep(pen_trt,each=cattle_per_pen)

  #design 4
  sim_4 <- dat_sim
  pen_trt <- sample(rep(1:2,each=pen_num/2),replace = F)
  
  for(c in 1:cohort_num){
    n1 <- (c-1)*cattle_per_cohort
    n2 <- c*cattle_per_cohort
    sim_4[(n1+1):n2,"pen"] <- sample(rep(1:pen_num,each=cattle_per_cohort/pen_num),replace = F)
  }
  sim_4 <- sim_4[order(sim_4$pen),]
  sim_4$trt <- rep(pen_trt,each=cattle_per_pen)
  # design 5
  sim_5 <- dat_sim
  pen_trt_trt <- c()
  for(pb in 1:(pen_num/2)){
    pen_trt_trt <- c(pen_trt_trt,sample(1:2))
  }
  pen_trt <- data.frame(pen=c(1:pen_num),trt=pen_trt_trt)  
  for(c in 1:cohort_num){
    n1 <- (c-1)*cattle_per_cohort
    n2 <- c*cattle_per_cohort
    sim_5[(n1+1):n2,"trt"] <- sample(rep(1:2,each=cattle_per_cohort/2),replace = F)
  }
  pen_trt1 <- pen_trt[which(pen_trt$trt==1),"pen"]
  pen_trt2 <- pen_trt[which(pen_trt$trt==2),"pen"]
  sim_5 <- sim_5[order(sim_5$trt),]
  sim_5[1:(total_cattle/2),"pen"] <- rep(pen_trt1,each=cattle_per_pen)
  sim_5[(total_cattle/2+1):total_cattle,"pen"] <- rep(pen_trt2,each=cattle_per_pen)
  ################### part 2: simulate the outcome #######################
  
  sim_list <- list(sim_1,sim_2,sim_3,sim_4,sim_5)
  pvalue <- c()
  coef_m1 <- c()
  coef_m2 <- c()
  coef_m3 <- c()
  for(i in 1:length(sim_list)){
    PenEffect <- data.frame(pen=c(1:pen_num),pen_effect=rnorm(pen_num,0,sqrt(var_pen))) 
    CohortEffect <- data.frame(cohort=c(1:cohort_num),cohort_effect=rnorm(cohort_num,0,sqrt(var_cohort))) 
    SubjectEffect <- data.frame(subject=c(1:total_cattle),subject_effect=rnorm(total_cattle,0,sqrt(var_subject))) 
    
    sim_i <- sim_list[[i]]
    
    sim_i <- merge(sim_i,PenEffect,by="pen")
    sim_i <- merge(sim_i,CohortEffect,by="cohort")
    sim_i <- merge(sim_i,SubjectEffect,by="subject")
    sim_i <- sim_i %>%
      mutate(trt =  ifelse(trt==1,"PLC",treatment)) %>%      # change the treatment factor name
      mutate(trt_effect = ifelse(trt=="PLC",PLC,treatment_value)) %>%  # added fixed effect column
      mutate(y = trt_effect+cohort_effect+pen_effect+subject_effect) # added all effect
    
    sim_i$trt <- factor(sim_i$trt,levels = c("PLC",treatment))
    ################### part 3: analyze the simulation data #######################
    m1 <- lm(y ~ trt, data=sim_i)
    
    coef_1 <- summary(m1)$coef
    resid_1 <- var(resid(m1))
    #trt_2,se,var_resid
    coef_m1 <- c(coef_m1,coef_1[2,1],coef_1[2,2],resid_1) 
    
    m2 <- lmer(y ~ trt+(1|pen), data=sim_i)
    coef_2 <- summary(m2)$coefficients
    ranef_2 <- as.data.frame(summary(m2)$varcor)
    ranef_2 <- ranef_2$vcov # pen,Residual
    # (trt_2,se,pen,Residual)
    coef_m2 <- c(coef_m2,coef_2[2,1],coef_2[2,2],ranef_2)
    
    if(i==2|i==3){
      coef_m3 <- c(coef_m3,rep(NA,5))
      coef_3[2,5] <- NA
      # (trt_2,se,ranef_cohort,ranef_pen,ranef_subject)
    }else{
      m3 <- lmer(y ~ trt+(1|cohort)+(1|pen), data=sim_i)
      coef_3 <- summary(m3)$coefficients
      ranef_3 <- as.data.frame(summary(m3)$varcor)
      ranef_3 <- ranef_3$vcov
      coef_m3 <- c(coef_m3,coef_3[2,1],coef_3[2,2],ranef_3)
    }
    
    pvalue <- c(pvalue,coef_1[2,4],coef_2[2,5],coef_3[2,5])
     
  }
  
  pvalue <- ifelse(pvalue<=0.05,1,0)
  return(c(pvalue,coef_m1,coef_m2,coef_m3))
}

######################### output: sig table ###########################
sig_table <- table_from_sim[,1:15]

sig_table <- sig_table %>%
  apply(2,sum)/nrep %>%
  round(4)

sig_table <- matrix(sig_table,nrow = 5,byrow = T)
colnames(sig_table) = c("m1","m2","m3")
rownames(sig_table) = c("design 1",'design 2','design 3','design 4','design 5')

write.csv(sig_table,paste0(treatment,"_sig_table.csv"))

######################### output: coef of m1 ###########################
res_m1 <- table_from_sim[,16:30]
res_m1 <- res_m1 %>%
  apply(2,mean) %>%
  round(4) %>%
  matrix(nrow = 5,byrow = T)%>%
  as.data.frame(.)

true_row <- c(treatment_value-PLC,NA,var_cohort+var_pen+var_subject)
res_m1 <- rbind(res_m1,true_row)
colnames(res_m1) = c("mean_diff","s.e.","var_resid")
rownames(res_m1) = c("design 1",'design 2','design 3','design 4','design 5','true value')
res_m1$var_resid <- round(res_m1$var_resid,2)
write.csv(res_m1, paste0(treatment,'_coef_table_m1.csv'))
######################### output: coef of m2 ###########################
res_m2 <- table_from_sim[,31:50]

res_m2 <- res_m2 %>%
  apply(2,mean) %>%
  round(4) %>%
  matrix(nrow = 5,byrow = T) %>%
  as.data.frame(.)

true_row <- c(treatment_value-PLC,NA,var_pen,var_subject+var_cohort)
res_m2 <- rbind(res_m2,true_row)

colnames(res_m2) = c("mean_diff","s.e.","var_pen","var_residual")

res_m2 <- res_m2 %>%
  mutate(var_pen=round(var_pen,2)) %>%
  mutate(var_residual=round(var_residual,2))


rownames(res_m2) = c("design 1",'design 2','design 3','design 4','design 5','true value')  
write.csv(res_m2,paste0(treatment,'_coef_table_m2.csv'))


######################### output: coef of m3 ###########################
res_m3 <- table_from_sim[,51:75]

res_m3 <- res_m3 %>%
  apply(2,mean) %>%
  round(4) %>%
  matrix(nrow = 5,byrow = T) %>%
  as.data.frame(.)
true_row <- c(treatment_value-PLC,NA,var_cohort,var_pen,var_subject)
res_m3 <- rbind(res_m3,true_row)

colnames(res_m3) = c("mean_diff","s.e.","var_cohort","var_pen","var_subject")

res_m3 <- res_m3 %>%
  mutate(var_cohort=round(var_cohort,2)) %>%
  mutate(var_pen=round(var_pen,2)) %>%
  mutate(var_subject=round(var_subject,2))


rownames(res_m3) = c("design 1",'design 2','design 3','design 4','design 5','true value')  
write.csv(res_m3,paste0(treatment,'_coef_table_m3.csv'))

# res_m2 <- res_m2 %>%
#   mutate(var_cohort=round(var_cohort,2)) %>%
#   mutate(var_pen=round(var_pen,2)) %>%
#   mutate(var_subject=round(var_subject,2)) %>%
#   mutate(total_var=var_cohort+var_pen+var_subject) %>%
#   mutate(percent_cohort=var_cohort/total_var) %>%
#   mutate(percent_cohort=paste0(round(percent_cohort*100,2),"%")) %>%
#   mutate(percent_pen = var_pen/total_var) %>%
#   mutate(percent_pen=paste0(round(percent_pen*100,2),"%")) %>%
#   mutate(percent_subject= var_subject/total_var) %>%
#   mutate(percent_subject=paste0(round(percent_subject*100,2),"%")) %>%
#   mutate(mean=coef_Intercept+coef_trt) %>%
#   mutate(mean_PLC=coef_Intercept) %>%
#   mutate(total_var=NULL) %>%
#   mutate(coef_Intercept=NULL) %>%
#   mutate(coef_trt=NULL)
#   
# colnames(res_m2)[which(names(res_m2) == "mean")] <- paste0("mean_",treatment)