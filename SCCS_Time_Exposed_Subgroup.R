
library('data.table')
library('bit64')
library('gnm')
library('lmtest')
library('ggplot2')
library('dplyr')
`%!in%` <- Negate(`%in%`)


sccs <-fread("subgroups/single_expo_30pre_3060post/full_cohort/full_cohort.csv")
sccs1 <- sccs[ENTRY_COHORT_IsExcluded==FALSE, ]


compare_patient_time <-function(cohort, cutoff){
  
  sccs_e <- sccs1[SELF_CONTROLLED_PERIOD_RISK_1==TRUE,]
  sccs_ints <- as.data.frame(sccs_e %>% group_by(PID) %>% transmute(Total=sum(INTERVAL_LENGTH)))
  short_PIDS <- subset(sccs_ints, Total<=cutoff)$PID
  long_PIDS <- subset(sccs_ints, Total>cutoff)$PID
  short_cohort <- sccs1[ (sccs1$PID %in% short_PIDS), ]
  long_cohort <- sccs1[ (sccs1$PID %!in% short_PIDS), ]
  print("the length of short_cohort is ")
  print(length(short_cohort$PID))
  print("the length of long_cohort is ")
  print(length(long_cohort$PID))
  
  short_m = gnm(OUTCOME_CASES_0 ~ SELF_CONTROLLED_PERIOD_RISK_1 +
                  SELF_CONTROLLED_PERIOD_PRERISK +
                  SELF_CONTROLLED_PERIOD_WASHOUT_1 + 
                  SELF_CONTROLLED_PERIOD_WASHOUT_2,
                eliminate=factor(PID), data=short_cohort, family=poisson, offset=I(log(INTERVAL_LENGTH)))
  
  long_m = gnm(OUTCOME_CASES_0 ~ SELF_CONTROLLED_PERIOD_RISK_1 +
                 SELF_CONTROLLED_PERIOD_PRERISK +
                 SELF_CONTROLLED_PERIOD_WASHOUT_1 + 
                 SELF_CONTROLLED_PERIOD_WASHOUT_2,
               eliminate=factor(PID), data=long_cohort, family=poisson, offset=I(log(INTERVAL_LENGTH)))
  
  return(list("short" = exp(coef(short_m)), "long" = exp(coef(long_m))))
}
  

compare_intervals <-function(cohort, cutoff){
  short_cohort <-  cohort[ (INTERVAL_LENGTH <= cutoff & SELF_CONTROLLED_PERIOD_RISK_1 == TRUE) 
                           | SELF_CONTROLLED_PERIOD_RISK_1 == FALSE,]
  long_cohort <- cohort[ (INTERVAL_LENGTH > cutoff & SELF_CONTROLLED_PERIOD_RISK_1 == TRUE) 
                         | SELF_CONTROLLED_PERIOD_RISK_1 == FALSE,]
  short_m = gnm(OUTCOME_CASES_0 ~ SELF_CONTROLLED_PERIOD_RISK_1 +
                  SELF_CONTROLLED_PERIOD_PRERISK +
                  SELF_CONTROLLED_PERIOD_WASHOUT_1 + 
                  SELF_CONTROLLED_PERIOD_WASHOUT_2,
           eliminate=factor(PID), data=short_cohort, family=poisson, offset=I(log(INTERVAL_LENGTH)))
  
  long_m = gnm(OUTCOME_CASES_0 ~ SELF_CONTROLLED_PERIOD_RISK_1 +
                  SELF_CONTROLLED_PERIOD_PRERISK +
                  SELF_CONTROLLED_PERIOD_WASHOUT_1 + 
                  SELF_CONTROLLED_PERIOD_WASHOUT_2,
           eliminate=factor(PID), data=long_cohort, family=poisson, offset=I(log(INTERVAL_LENGTH)))
  
  return(list("short" = exp(coef(short_m)), "long" = exp(coef(long_m))))
}


short_long_list <-function(comparer, cohort, cutoffs){
  sl_list =list()
  for (i in cutoffs){
    sl_list[[i]] <- comparer(cohort, i) 
  }
  return(sl_list)
}

sl_list_to_df <-function(sl, interval_to_graph){
  short_list <- list()
  long_list <- list()
  for(i in c(1:length(sl))){
    short_list[i] <- sl[[i]]["short"][[1]][interval_to_graph]
    long_list[i] <-sl[[i]]["long"][[1]][interval_to_graph]
  }
  df <- data.frame(index = c(1:length(short_list)),short = unlist(short_list), long =unlist(long_list))
  return(df)
}



graph_df <- function(df){
  print(df)
  plot(x = df$index, y =df$short, type ="l", col ="blue", ylim = c(0, 3),
       main= "Relative Incidence for Subgroup of Cohort Above and Below Cutoff Point",
       ylab="Relative Incidence", xlab="Cutoff Point in Days (blue <= x < red)" )
  lines(x =df$index, y=df$long,type="l", col = "red")
}


c <- c(1:1000)



sl <- short_long_list(comparer = compare_patient_time, sccs1, c )
sl

sl_df = sl_list_to_df(sl, "SELF_CONTROLLED_PERIOD_RISK_1TRUE")

graph_df(sl_df)

