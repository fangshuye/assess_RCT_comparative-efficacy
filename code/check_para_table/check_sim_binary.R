#####################SIMULATION: 10000 FOR PLC vs TIL#########################
rm(list = ls())
library(locfit)
library(lme4)
library(dplyr)
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
var_cohort <- 0.56
var_pen <- 0.83

OR <- exp(-0.3759)
risk_Baseline <- OR/(1+OR)
risk_PLC <- risk_Baseline
risk_TMS <- risk_Baseline/1.16
# PLC
Odd_PLC <- risk_PLC/(1-risk_PLC)
PLC <- log(Odd_PLC)
# TMS
Odd_TMS <- risk_TMS/(1-risk_TMS)
TMS <- round(log(Odd_TMS),4)  

# choose trt
treatment <- "TMS"
treatment_value <- TMS
Odd_trt <- Odd_TMS

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
    ################### part 3: analyze the simulation data #######################
    m1 <- glm(y ~ trt, data=sim_i, family = "binomial")
    
    coef_1 <- summary(m1)$coef
    #log_odd_ratio,se,z value, 
    coef_m1 <- c(coef_m1,coef_1[2,1:3]) 
    
    # (odd_ratio,se,z,ranef_pen)
    m2 <-  glmer(y ~ trt+(1|pen), data=sim_i, family = "binomial")
    coef_2 <- summary(m2)$coefficients
    ranef_2 <- unlist(summary(m2)$varcor)
    coef_m2 <- c(coef_m2,coef_2[2,1:3],ranef_2)
    
    
    if(i==2|i==3){
      coef_m3 <- c(coef_m3,rep(NA,5))
      coef_3[2,4] <- NA
    }else{
      m3 <-  glmer(y ~ trt+(1|cohort)+(1|pen), data=sim_i, family = "binomial")
      coef_3 <- summary(m3)$coefficients
      ranef_3 <- unlist(summary(m3)$varcor)
      
      # (odd_ratio,se,z,ranef_cohort,ranef_pen)
      coef_m3 <- c(coef_m3,coef_3[2,1:3],ranef_3)
    }
    pvalue <- c(pvalue,coef_1[2,4],coef_2[2,4],coef_3[2,4])
  }
  
  pvalue <- ifelse(pvalue<=0.05,1,0)
  return(c(pvalue,coef_m1,coef_m2,coef_m3))
}

write.csv(table_from_sim,"check_sim_binary_all_data.csv")

######################### output: sig table ###########################
sig_table <- table_from_sim[,1:15]

sig_table <- sig_table %>%
  apply(2,sum)/nrep %>%
  round(4)

sig_table <- matrix(sig_table,nrow = 5,byrow = T)
colnames(sig_table) = c("M1","M2","M3")
rownames(sig_table) = c("W1",'W2','W3','W4','W5')

write.csv(sig_table,paste0(treatment,"_sig_table.csv"))

#################### para #########################
para <- table_from_sim[,16:75]
mean_para <- para %>%
  apply(2,mean) %>%
  round(4)
var_para <- para %>%
  apply(2,var)

# mean table
mean_para_m1 <- mean_para[1:15] %>%
  matrix(ncol = 5,byrow = F)%>%
  as.data.frame(.)
rownames(mean_para_m1) <- c("m1_logoddratio_mean","m1_se_mean","m1_zvalue_mean")

mean_para_m2 <- mean_para[16:35] %>%
  matrix(ncol = 5,byrow = F)%>%
  as.data.frame(.)
rownames(mean_para_m2) <- c("m2_logoddratio_mean","m2_se_mean","m2_tzvalue_mean","m2_varpen_mean")

mean_para_m3 <- mean_para[36:60] %>%
  matrix(ncol = 5,byrow = F)%>%
  as.data.frame(.)
rownames(mean_para_m3) <- c("m3_logoddratio_mean","m3_se_mean","m3_tzvalue_mean","m3_varcohort_mean","m3_varpen_mean")

mean_table <- rbind(mean_para_m1,
                    mean_para_m2,
                    mean_para_m3)
# var table
var_para_m1 <- var_para[1:15] %>%
  matrix(ncol = 5,byrow = F)%>%
  as.data.frame(.)
rownames(var_para_m1) <- c("m1_logoddratio_var","m1_se_var","m1_zvalue_var")

var_para_m2 <- var_para[16:35] %>%
  matrix(ncol = 5,byrow = F)%>%
  as.data.frame(.)
rownames(var_para_m2) <- c("m2_logoddratio_var","m2_se_var","m2_tzvalue_var","m2_varpen_var")

var_para_m3 <- var_para[36:60] %>%
  matrix(ncol = 5,byrow = F)%>%
  as.data.frame(.)
rownames(var_para_m3) <- c("m3_logoddratio_var","m3_se_var","m3_tzvalue_var","m3_varcohort_var","m3_varpen_var")

var_table <- rbind(var_para_m1,
                    var_para_m2,
                    var_para_m3)

output <- rbind(mean_table,var_table)
output <- output[order(row.names(output)),]

colnames(output) <- paste0("W",c(1:5))
log_odd_ratio_round <- round(log(Odd_trt/Odd_PLC),4)
true_col <- c(log_odd_ratio_round,NA,
              NA,NA,
              NA,NA,
              log_odd_ratio_round,NA,
              NA,NA,
              NA,NA,
              var_pen,NA,
              log_odd_ratio_round,NA,
              NA,NA,
              NA,NA,
              var_cohort,NA,
              var_pen,NA)
output[,6] <- true_col
output <- output[-c(14,22,24),]
write.csv(output,"check_sim_binary_para.csv")

# plot the distribution of z
temp_table <- table_from_sim[,-c(1:15)]
z_value_m1_table <- temp_table[,seq(3,15,3)] %>%
  as.data.frame(.)
temp_table <- temp_table[,-c(1:15)]
z_value_m2_table <- temp_table[,seq(3,20,4)]%>%
  as.data.frame(.)
temp_table <- temp_table[,-c(1:20)]
z_value_m3_table <- temp_table[,seq(3,25,5)]%>%
  as.data.frame(.)

colnames(z_value_m1_table) <- paste0("W",1:5)
colnames(z_value_m2_table) <- paste0("W",1:5)
colnames(z_value_m3_table) <- paste0("W",1:5)

# reshape the data
library(tidyr)
m1_plot <- gather(z_value_m1_table)
m2_plot <- gather(z_value_m2_table)
m3_plot <- gather(z_value_m3_table)


p1 <- ggplot(m1_plot,aes(x=value,color=key))+geom_density()+labs(x="Z value",color="Workflow")+
  geom_vline(xintercept = qnorm(0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  geom_vline(xintercept = qnorm(1-0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  scale_x_continuous(breaks = c(qnorm(0.025),0,qnorm(1-0.025)),labels = c("-1.96","0","1.96"))+
  scale_colour_manual(values = c("mediumpurple3","deepskyblue","orangered3","goldenrod1","forestgreen"))

p2 <- ggplot(m2_plot,aes(x=value,color=key))+geom_density()+labs(x="Z value",color="Workflow")+
  geom_vline(xintercept = qnorm(0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  geom_vline(xintercept = qnorm(1-0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  scale_x_continuous(breaks = c(qnorm(0.025),0,qnorm(1-0.025)),labels = c("-1.96","0","1.96"))+
  scale_colour_manual(values = c("mediumpurple3","deepskyblue","orangered3","goldenrod1","forestgreen"))

p3 <- ggplot(m3_plot,aes(x=value,color=key))+geom_density()+labs(x="Z value",color="Workflow")+
  geom_vline(xintercept = qnorm(0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  geom_vline(xintercept = qnorm(1-0.025) , linetype="dashed", 
             color = "red", size=1.5)+
  scale_x_continuous(breaks = c(qnorm(0.025),0,qnorm(1-0.025)),labels = c("-1.96","0","1.96"))+
  scale_colour_manual(values = c("mediumpurple3","goldenrod1","forestgreen"))

ggsave("binary_z_value_m1.png",
       p1)
ggsave("binary_z_value_m2.png",
       p2)
ggsave("binary_z_value_m3.png",
       p3)





