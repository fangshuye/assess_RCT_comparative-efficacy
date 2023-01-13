#####################SIMULATION: Binary outcome#########################
#####################Pen level analysis#########################
#####################M1 and M2#########################

rm(list = ls())
library(locfit)
library(lme4)
library(dplyr)
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
var_cohort <- 0.56
var_pen <- 0.83

OR <- exp(-0.3759)
risk_Baseline <- OR/(1+OR)
risk_PLC <- risk_Baseline
risk_TIL <- risk_Baseline/2.32
risk_ENR <- risk_Baseline/1.50
risk_TMS <- risk_Baseline/1.16
# PLC
Odd_PLC <- risk_PLC/(1-risk_PLC)
PLC <- log(Odd_PLC)
# TIL
Odd_TIL <- risk_TIL/(1-risk_TIL)
TIL <- round(log(Odd_TIL),4) 
# ENR
Odd_ENR <- risk_ENR/(1-risk_ENR)
ENR <- round(log(Odd_ENR),4) 
# TMS
Odd_TMS <- risk_TMS/(1-risk_TMS)
TMS <- round(log(Odd_TMS),4)  
# NE: no effect
Odd_NE <- Odd_PLC
NE <- PLC 

# pull together
treatments <- c("TIL","ENR","TMS","NE")
treatment_values <- c(TIL,ENR,TMS,NE)
Odd_trts <- c(Odd_TIL,Odd_ENR,Odd_TMS,Odd_NE)
# choose trt
treatment <- treatments[4]
treatment_value <- treatment_values[4]
Odd_trt <- Odd_trts[4]

##############################  simulation begin ###############################
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
  
  new <- c()
  for(i in 1:length(sim_list)){
    PenEffect <- data.frame(pen=c(1:pen_num),pen_effect=rnorm(pen_num,0,sqrt(var_pen))) 
    CohortEffect <- data.frame(cohort=c(1:cohort_num),cohort_effect=rnorm(cohort_num,0,sqrt(var_cohort))) 
    
    sim_i <- sim_list[[i]]
    
    sim_i <- merge(sim_i,PenEffect,by="pen")
    sim_i <- merge(sim_i,CohortEffect,by="cohort")
    sim_i <- sim_i %>%
      mutate(trt =  ifelse(trt==1,"PLC",treatment)) %>%      # change the treatment factor name
      mutate(trt_effect = ifelse(trt=="PLC",PLC,treatment_value)) %>%  # added fixed effect column
      mutate(added_effect = trt_effect+cohort_effect+pen_effect) %>% # added all effect
      mutate(risk = expit(added_effect))  # get the inverse of the logistic link function
    
    sim_i$y <- rbinom(n=total_cattle, size = 1, prob = sim_i$risk)
    sim_i$trt <- factor(sim_i$trt,levels = c("PLC",treatment))
    
    ##### create pen level dat #####
    sim_pen_i <- sim_i %>%
      select(pen,trt,y)
    sim_pen_i <- sim_pen_i %>%
      group_by(pen,trt) %>%
      summarise(event = sum(y),
                total = n())
    
    sim_i <- sim_i %>%
      select(cohort, pen, trt, y)
    ################### part 3: analyze the simulation data #######################
    # m_indi <- glm(y ~ trt, data=sim_i, family = "binomial")
    # m_pen <- glm(cbind(event,total-event) ~ trt, data=sim_pen_i, family = "binomial")
    # new <- c(new, c(summary(m_indi)$coef[2,1:2],summary(m_pen)$coef[2,1:2]))
    
    # (odd_ratio,se,ranef_pen)
    m_pen <-  glmer(cbind(event,total-event) ~ trt+(1|pen), data=sim_pen_i, family = "binomial")
    m_indi <-  glmer(y ~ trt+(1|pen), data=sim_i, family = "binomial")
    new <- c(new, c(summary(m_indi)$coef[2,1:2],summary(m_pen)$coef[2,1:2]))
    
    
  }
  
  pvalue <- ifelse(pvalue<=0.05,1,0)
  return(c(pvalue,coef_m1,coef_m2,coef_m3))
}

m1 <- new %>%
  round(4)
m1 <- matrix(m1,ncol = 2,byrow = T)
colnames(m1) = c("Effect size","Std.Error")
savedir <- "/Users/fangshu/Desktop/research/design/simulation/code/Pen_level_analysis/"
write.csv(m1,paste0(savedir,"M1.csv"))

