#####################SIMULATION: 10000 FOR PLC vs TIL#########################
library(locfit)
library(lme4)
library(dplyr)
library(doParallel)
set.seed(20200609, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)

# get fixed effect and random effect
sigma_pen <- sqrt(0.83)
sigma_cohort <- sqrt(0.56) 
sigma_subject <- sqrt(3.29)

OR <- exp(-0.3759)
risk_Baseline <- OR/(1+OR)
risk_PLC <- risk_Baseline
risk_TIL <- risk_Baseline/2.32
# PLC
PLC <- log(risk_PLC/(1-risk_PLC))
# TIL
TIL <- log(risk_TIL/(1-risk_TIL))


nrep <- 10000
dat_sim <- data.frame(subject=c(1:1200),cohort=rep(1:30,each=40))

table_from_sim <- foreach (rep=1:nrep, .combine = rbind)%dopar%{
  
  ################### part 1: simulate the design data #######################
  # design 1
  sim_1 <- dat_sim
  for(c in 1:30){
    n <- c*40
    sim_1[(n-39):n,"trt"] <- sample(rep(1:2,each=20),replace = F)
  }
  sim_1$random_row_num <- sample(1:1200,replace = FALSE)
  sim_1 <- sim_1[order(sim_1$random_row_num),]
  sim_1 <- sim_1 %>%
    mutate(pen=rep(1:30,each=40)) %>%
    mutate(random_row_num=NULL)
  
  # design 2
  sim_2 <- dat_sim
  sim_2$trt <- sample(rep(1:2,each=600),replace = F)
  sim_2$pen <- rep(1:30,each=40) # equal to cohort in this design
  
  # design 3
  sim_3 <- dat_sim
  sim_3$pen <- rep(1:30,each=40) # equal to cohort in this design
  pen_trt <- sample(rep(1:2,each=15),replace = F)
  sim_3$trt <- rep(pen_trt,each=40)

  # design 4
  # sim_4 <- dat_sim
  # for(c in 1:30){
  #   n <- c*40
  #   sim_4[(n-39):n,"pen"] <- sample(1:30,replace = T)
  # }
  
  # design 5 
  
  
  ################### part 2: simulate the outcome #######################
  
  sim_list <- list(sim_1,sim_2,sim_3)
    
  for(i in 1:length(sim_list)){
    PenEffect <- data.frame(pen=c(1:30),pen_effect=rnorm(30,0,sigma_pen)) 
    CohortEffect <- data.frame(cohort=c(1:30),cohort_effect=rnorm(30,0,sigma_cohort)) 
    SubjectEffect <- data.frame(subject=c(1:1200),subject_effect=rnorm(1200,0,sigma_subject)) 
    
    sim_i <- sim_list[[i]]
    
    sim_i <- merge(sim_i,PenEffect,by="pen")
    sim_i <- merge(sim_i,CohortEffect,by="cohort")
    sim_i <- merge(sim_i,SubjectEffect,by="subject")
    sim_i <- sim_i %>%
      mutate(trt =  ifelse(trt==1,"PLC","TIL")) %>%      # change the treatment factor name
      mutate(trt_effect = ifelse(trt=="PLC",PLC,TIL)) %>%  # added fixed effect column
      mutate(added_effect = trt_effect+cohort_effect+pen_effect+subject_effect) %>% # added all effect
      mutate(risk = expit(added_effect)) %>%  # get the inverse of the logistic link function
      select(subject,cohort,pen,trt,risk)
    sim_i$y <- rbinom(n=1200, size = 1, prob = sim_i$risk)
    
    sim_list[[i]] <- sim_i
  }
  
  
  ################### part 3: analyze the simulation data #######################
  
  # design 1
  m1 <- glm(y ~ trt, data=sim_list[[1]], family = "binomial")
  m2 <- glmer(y ~ trt+(1|cohort)+(1|pen), data=sim_list[[1]], family = "binomial")
  
  coef1_1 <- summary(m1)$coef
  ranef1_1 <- var(resid(m1))
  coef1_2 <- summary(m2)$coefficients
  ranef1_2 <- c(unlist(summary(m2)$varcor),var(resid(m2)))
  
  # design 2
  m1 <- glmer(y ~ trt+(1|cohort), data=sim_list[[2]], family = "binomial")
  m2 <- glmer(y ~ trt+(1|cohort)+(1|pen), data=sim_list[[2]], family = "binomial")
  
  coef2_1 <- summary(m1)$coefficients
  ranef2_1 <- c(unlist(summary(m1)$varcor),var(resid(m1)))
  
  coef2_2 <- summary(m2)$coefficients
  ranef2_2 <- c(unlist(summary(m2)$varcor),var(resid(m2)))
  
  # design 3
  
  m <- glmer(y ~ trt+(1|cohort)+(1|pen), data=sim_list[[3]], family = "binomial")
  coef3 <- summary(m)$coefficients
  ranef3 <- c(unlist(summary(m)$varcor),var(resid(m)))
  
  # design 4 and 5
  
  # m1 <- glm(y ~ trt+(1|pen), data=dat_sim, family = "binomial")
  # m2 <- glmer(y ~ trt+(1|cohort)+(1|pen), data=dat_sim, family = "binomial")
  
  pvalue <- c(coef1_1[2,4],coef1_2[2,4], # design 1 (m1,m2)
              coef2_1[2,4],coef2_2[2,4], # design 2 (m1,m2)
              coef3[2,4],coef3[2,4] # design 3 (m1,m2)
              )
  pvalue <- ifelse(pvalue<=0.05,1,0)
  
  coef_m1 <- c(coef1_1[,1],ranef1_1, # design 1 (Intercept,trtTIL,ranef_subject)
               coef2_1[,1],ranef2_1, # design 2 (Intercept,trtTIL,ranef_cohort,ranef_subject)
               coef3[,1],ranef3 # design 3 (Intercept,trtTIL,ranef_cohort,ranef_pen,ranef_subject)
               )
  coef_m2 <- c(coef1_2[,1],ranef1_2, # design 1 (Intercept,trtTIL,ranef_cohort,ranef_pen,ranef_subject)
               coef2_2[,1],ranef2_2, # design 2 (Intercept,trtTIL,ranef_cohort,ranef_pen,ranef_subject)
               coef3[,1],ranef3 # design 3 (Intercept,trtTIL,ranef_cohort,ranef_pen,ranef_subject)
               )
  
  return(c(pvalue,coef_m1,coef_m2))
}

######################### output: sig table ###########################
sig_table <- table_from_sim[,1:6]

sig_table <- sig_table %>%
  apply(2,sum)/nrep %>%
  round(4)

sig_table <- matrix(sig_table,nrow = 3,byrow = T)
colnames(sig_table) = c("m1","m2")
rownames(sig_table) = c("design 1",'design 2','design 3')

write.csv(sig_table,'sig_table.csv')

######################### output: coef of m1 ###########################
res_m1 <- table_from_sim[,7:18]
res_m1 <- apply(res_m1,2, mean)
coef_table_m1 <- matrix(NA,nrow = 3,ncol = 5)
colnames(coef_table_m1) = c("coef_Intercept","coef_trtTIL","var_cohort","var_pen","var_resid")
rownames(coef_table_m1) = c("design 1",'design 2','design 3')
coef_table_m1[1,c(1,2,5)] <- res_m1[1:3]
coef_table_m1[2,c(1,2,3,5)] <- res_m1[4:7]
coef_table_m1[3,] <- res_m1[8:12]

write.csv(coef_table_m1,'coef_table_m1.csv')
######################### output: coef of m2 ###########################
res_m2 <- table_from_sim[,19:33]

res_m2 <- res_m2 %>%
  apply(2,mean) %>%
  round(2) %>%
  matrix(nrow = 3,byrow = T)

colnames(res_m2) = c("coef_Intercept","coef_trtTIL","var_cohort","var_pen","var_resid")
rownames(res_m2) = c("design 1",'design 2','design 3')

write.csv(res_m2,'coef_table_m2.csv')
