#####################SIMULATION: 10000 FOR PLC vs TIL#########################
rm(list = ls())
library(lme4)
library(dplyr)
library(lmerTest)
library(doParallel)
library(ggplot2)
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
Trt3 <- 903.274-22.507

treatment <- "Trt3"
treatment_value <- Trt3

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
  m1_fix <- c()
  m1_se <- c()
  m1_T <- c()
  m1_df <- c()
  m1_res <- c()
  
  m2_fix <- c()
  m2_se <- c()
  m2_T <- c()
  m2_df <- c()
  m2_pen <- c()
  m2_res <- c()
  
  m3_fix <- c()
  m3_se <- c()
  m3_T <- c()
  m3_df <- c()
  m3_cohort <- c()
  m3_pen <- c()
  m3_res <- c()

  
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
    # m1
    m1 <- lm(y ~ trt, data=sim_i)
    coef_1 <- summary(m1)$coef
    # m2
    m2 <- lmer(y ~ trt+(1|pen), data=sim_i)
    coef_2 <- summary(m2)$coefficients
    ranef_2 <- as.data.frame(summary(m2)$varcor)
    ranef_2 <- ranef_2$vcov # pen,Residual
    # m3
    if(i==2|i==3){
      m3_fix <- c(m3_fix,NA)
      m3_se <- c(m3_se,NA)
      m3_T <- c(m3_T,NA)
      m3_df <- c(m3_df,NA)
      m3_cohort <- c(m3_cohort,NA)
      m3_pen <- c(m3_pen,NA)
      m3_res <- c(m3_res,NA)
      coef_3[2,5] <- NA
    }else{
      m3 <- lmer(y ~ trt+(1|cohort)+(1|pen), data=sim_i)
      coef_3 <- summary(m3)$coefficients
      ranef_3 <- as.data.frame(summary(m3)$varcor)
      ranef_3 <- ranef_3$vcov
      m3_fix <- c(m3_fix,coef_3[2,1])
      m3_se <- c(m3_se,coef_3[2,2])
      m3_T <- c(m3_T,coef_3[2,4])
      m3_df <- c(m3_df,coef_3[2,3])
      m3_cohort <- c(m3_cohort,ranef_3[1])
      m3_pen <- c(m3_pen,ranef_3[2])
      m3_res <- c(m3_res,ranef_3[3])
    }
    
    pvalue <- c(pvalue,coef_1[2,4],coef_2[2,5],coef_3[2,5])
    m1_fix <- c(m1_fix,coef_1[2,1])
    m1_se <- c(m1_se,coef_1[2,2])
    m1_T <- c(m1_T,coef_1[2,3])
    m1_df <- c(m1_df,m1$df.residual)
    m1_res <- c(m1_res,var(resid(m1)))
    
    m2_fix <- c(m2_fix,coef_2[2,1])
    m2_se <- c(m2_se,coef_2[2,2])
    m2_T <- c(m2_T,coef_2[2,4])
    m2_df <- c(m2_df,coef_2[2,3])
    m2_pen <- c(m2_pen,ranef_2[1])
    m2_res <- c(m2_res,ranef_2[2])
    
  }
  pvalue <- ifelse(pvalue<=0.05,1,0)
  return(c(pvalue, m1_fix, m1_se, m1_T,m1_df, m1_res, 
           m2_fix, m2_se, m2_T, m2_df, m2_pen, m2_res,
           m3_fix,  m3_se, m3_T, m3_df, m3_cohort, m3_pen, m3_res))
}
write.csv(table_from_sim,"check_sim_continuous_all_data.csv")

######################### output: sig table ###########################
sig_table <- table_from_sim[,1:15]

sig_table <- sig_table %>%
  apply(2,sum)/nrep %>%
  round(4)

sig_table <- matrix(sig_table,nrow = 5,byrow = T)
colnames(sig_table) = c("M1","M2","M3")
rownames(sig_table) = c("W1",'W2','W3','W4','W5')

write.csv(sig_table,paste0(treatment,"_sig_table.csv"))

#################### para table #############################

para_table <- table_from_sim[,16:105]

para_table_mean <- para_table %>%
  apply(2,mean) %>%
  matrix(ncol = 5,byrow = T)%>%
  as.data.frame(.)
  
colnames(para_table_mean) <- paste0("W",1:5)
rownames(para_table_mean) <- paste0(c("m1_fix", "m1_se", "m1_T","m1_df", "m1_res", 
                                      'm2_fix', "m2_se", "m2_T", "m2_df", "m2_pen", "m2_res",
                                      "m3_fix",  "m3_se", "m3_T", "m3_df", "m3_cohort", "m3_pen", "m3_res"),"_mean")  
  
# s.e. part
para_table_var <- para_table %>%
  apply(2,var) %>%
  matrix(ncol = 5,byrow = T)%>%
  as.data.frame(.)

colnames(para_table_var) <- paste0("W",1:5)
rownames(para_table_var) <- paste0(c("m1_fix", "m1_se", "m1_T","m1_df", "m1_res", 
                                      'm2_fix', "m2_se", "m2_T", "m2_df", "m2_pen", "m2_res",
                                      "m3_fix",  "m3_se", "m3_T", "m3_df", "m3_cohort", "m3_pen", "m3_res"),"_var")  

para_table_output <- rbind(para_table_mean,para_table_var[c(1:3,6:8,12:14),])
true_col <- c(treatment_value-PLC,rep(NA,3),var_cohort+var_pen+var_subject,
              treatment_value-PLC,rep(NA,3),var_pen,var_subject+var_cohort,
              treatment_value-PLC,rep(NA,3),var_cohort,var_pen,var_subject,
              rep(NA,9))

para_table_output[,6] <- true_col

# re-order the rows
para_table_output <- para_table_output[c(1,19,2,20,3,21,4,5,
                                         6,22,7,23,8,24,9:11,
                                         12,25,13,26,14,27,15:18),]
write.csv(para_table_output,"check_sim_continuous_para.csv")

# plot the t value
temp_table <- table_from_sim[,-c(1:15)]
T_value_m1_table <- temp_table[,11:15] %>%
  as.data.frame(.)
T_value_m2_table <- temp_table[,36:40] %>%
  as.data.frame(.)
T_value_m3_table <- temp_table[,66:70] %>%
  as.data.frame(.)
colnames(T_value_m1_table) <- paste0("W",1:5)
colnames(T_value_m2_table) <- paste0("W",1:5)
colnames(T_value_m3_table) <- paste0("W",1:5)

# reshape the data
library(tidyr)
m1_plot <- gather(T_value_m1_table)
m2_plot <- gather(T_value_m2_table)
m3_plot <- gather(T_value_m3_table)

p1 <- ggplot(m1_plot,aes(x=value,color=key))+geom_density()+labs(x="T value",color="Workflow")+
  geom_vline(xintercept = qnorm(0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  geom_vline(xintercept = qnorm(1-0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  scale_x_continuous(breaks = c(qt(0.025,798),0,qt(1-0.025,798)),labels = c("-1.96","0","1.96"))+
  scale_colour_manual(values = c("mediumpurple3","deepskyblue","orangered3","goldenrod1","forestgreen"))

ggsave("continuous_t_value_m1.png",
       p1)
