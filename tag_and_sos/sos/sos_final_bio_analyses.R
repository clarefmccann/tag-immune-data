################################################# 
####### SOS (TAG) bio data final analyses ####### 
################################################# 

library(utils)
library(tidyverse)
library(panelr)
library(readxl)
library(ggplot2)
library(viridis)
library(ggpubr)
library(rstatix)
library(lmerTest)
library(stringr)
library(Hmisc)
library(corrplot)

workdir <- 'S:/MNHS-Psych/ANDL-Lab-TAG-Study/SOS Study'
HR_Resting1 <- read.csv(file.path(workdir,"/physio/Mean_HR_Resting1.csv", fsep=""))
HR_Resting2 <- read.csv(file.path(workdir,"/physio/Mean_HR_Resting2.csv", fsep=""))
HR_Task <- read.csv(file.path(workdir,"/physio/Mean_HR_Task.csv", fsep=""))

HR_conditions_raw <- read_excel(file.path(workdir,"/physio/Heart_Rate_Data_Task_Sequence_All_Slices_Conditions.xlsx", fsep=""))

sos_crp <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_CRP_cleaned.csv", fsep="")) 
sos_cytokines <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_cytokines_cleaned.csv", fsep=""))
sos_hormones <- read.csv(file.path(workdir,"/saliva/Cleaned_datasets/SOS_hormones_cleaned.csv", fsep=""))

####### Heart rate ####### 

colnames(HR_Resting1)[1] <- "ID"
colnames(HR_Resting2)[1] <- "ID"
colnames(HR_Task)[1] <- "ID"


# get rid of the repeat 313
HR_Resting1 <- filter(HR_Resting1, ID!="TAG313p2_20191117_151413")
HR_Resting2 <- filter(HR_Resting2, ID!="TAG313p2_20191117_151413")
HR_Task <- filter(HR_Task, ID!="TAG313p2_20191117_151413")

# There are two subjects (077 and 188) who had very low task HRs (high 30s) 
# despite having normal range resting state HRs, so they need to be removed
# as this was likely artifact during the task.

HR_Resting1 <- filter(HR_Resting1, ID!="TAG077_20190823_145422")
HR_Resting2 <- filter(HR_Resting2, ID!="TAG077_20190823_145422")
HR_Task <- filter(HR_Task, ID!="TAG077_20190823_145422")

HR_Resting1 <- filter(HR_Resting1, ID!="TAG188_20190731_145115")
HR_Resting2 <- filter(HR_Resting2, ID!="TAG188_20190731_145115")
HR_Task <- filter(HR_Task, ID!="TAG188_20190731_145115")

HR_Resting1$sequence <- rep(c("Resting 1"), each=50)
HR_Resting2$sequence <- rep(c("Resting 2"), each=49)
HR_Task$sequence <- rep(c("Task"), each=50)

HR_all_long <- bind_rows(HR_Resting1, HR_Resting2, HR_Task)

# average both resting states (for ID 080, use resting 1 as average because
# they don't have resting 2)

HR_all_wide <- spread(HR_all_long, sequence, Mean_HR)
HR_all_wide$RestingAvg <- ifelse(is.na(HR_all_wide$'Resting 2'), 
                                      HR_all_wide$'Resting 1', 
                                      rowMeans(HR_all_wide[,c('Resting 1', 'Resting 2')]))

# Add the individual task conditions:
HR_conditions <- pivot_longer(HR_conditions_raw, SOS068_20190303_150846:TAG223_20191116_145342, names_to = "ID", values_to = "HR")

HR_conditions$Task_condition <- as.factor(HR_conditions$Task_condition)

HR_conditions$ID <- as.factor(HR_conditions$ID)

HR_conditions_avg <- HR_conditions %>% group_by(ID, Task_condition) %>% summarise(mean_HR = mean(HR))

HR_conditions_clean <- HR_conditions_avg %>% filter(Task_condition!="NA")
HR_conditions_clean$Task_condition <- factor(HR_conditions_clean$Task_condition)
HR_conditions_clean_tbl <- as_tibble(HR_conditions_clean)
# Remove two outliers - their task HR data was impossibly low (probably artifact)
HR_conditions_clean_tbl <- filter(HR_conditions_clean_tbl, ID!="TAG077_20190823_145422")
HR_conditions_clean_tbl <- filter(HR_conditions_clean_tbl, ID!="TAG188_20190731_145115")
# and remove the repeat 313 again from this version of the long HR data
HR_conditions_clean_tbl <- filter(HR_conditions_clean_tbl, ID!="TAG313p2_20191117_151413")

# Averages by condition
meanhr_resting1 <- mean(HR_all_wide$'Resting 1') # 73.439
meanhr_resting2 <- mean(HR_all_wide$'Resting 2', na.rm = TRUE) # 72.371
meanhr_restingavg <- mean(HR_all_wide$RestingAvg) # 72.927
meanhr_task <- mean(HR_all_wide$Task) # 76.942
meanhr_R <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="R"), ]$mean_HR) # 76.962
meanhr_other <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="1"), ]$mean_HR) # 76.009
meanhr_selfdis <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="2"), ]$mean_HR) # 77.145
meanhr_selfcon <- mean(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR) # 77.519

# Check significant differences (Question: should this be paired t-tests or repeated measures ANOVAs?)

# paired t-tests
# Average all task vs. average resting state 
t.test(HR_all_wide$Task, HR_all_wide$RestingAvg, alternative = "greater", 
       paired = TRUE) # t = 7.4707, df = 49, p-value = 6.21e-10

# Self-connected vs. R condition (or should it be general resting state?) 
t.test(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR, 
       HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="R"), ]$mean_HR,
       alternative = "greater", paired = TRUE) # t = 2.8223, df = 49, p-value = 0.003435

# Self-connected vs. self-disconnected
t.test(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR, 
       HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="2"), ]$mean_HR,
       alternative = "greater", paired = TRUE) # t = 2.187, df = 49, p-value = 0.01677

# Self-connected vs. Other
t.test(HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="3"), ]$mean_HR, 
       HR_conditions_clean_tbl[which(HR_conditions_clean_tbl$Task_condition=="1"), ]$mean_HR,
       alternative = "greater", paired = TRUE) # t = 6.1557, df = 49, p-value = 6.731e-08

# Plotting Heart rate

# Task vs. resting state: HR_all_wide
HR_avg_long <- HR_all_wide[,c("ID","Task","RestingAvg")]
HR_avg_long <- pivot_longer(HR_avg_long, cols = c("Task","RestingAvg"), 
                            names_to = "condition", 
                            values_to = "mean_HR")
HR_avg_long$condition <- factor(HR_avg_long$condition,
                                                 levels = c("Task","RestingAvg"),
                                                 labels = c("SOS Task","Resting State Scan"))

compare_means(mean_HR ~ condition, data = HR_avg_long, method = "t.test", paired = TRUE)

my_comparisons1 <- list( c("SOS Task","Resting State Scan"))

get_box_stats <- function(y, upper_limit = max(HR_avg_long$mean_HR) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

avgHR_plot <- 
  ggplot(HR_avg_long, aes(x = condition, y = mean_HR)) +
  geom_violin(aes(fill = condition), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 9) +
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons1, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons order above
  ylab("Mean Heart Rate (bpm)") +
  xlab("Condition") +
  theme(legend.position = "none") 

# Within task: HR_conditions_clean_tbl
HR_conditions_clean_tbl$Task_condition <- factor(HR_conditions_clean_tbl$Task_condition,
                                                  levels = c("1","2","3","R"),
                                                  labels = c("Other", 
                                                             "Self-Disconnected",
                                                             "Self-Connected",
                                                             "Rest"))

compare_means(mean_HR ~ Task_condition, data = HR_conditions_clean_tbl, method = "t.test", paired = TRUE)

my_comparisons <- list( c("Self-Connected","Rest"),
                        c("Self-Connected","Self-Disconnected"),
                        c("Self-Connected","Other"),
                        c("Self-Disconnected","Other"),
                        c("Self-Disconnected","Rest"),
                        c("Other","Rest"))

get_box_stats <- function(y, upper_limit = max(HR_conditions_clean_tbl$mean_HR) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

taskHR_plot <- 
ggplot(HR_conditions_clean_tbl, aes(x = Task_condition, y = mean_HR)) +
  geom_violin(aes(fill = Task_condition), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 7.2) +
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons order above
  ylab("Mean Heart Rate (bpm)") +
  xlab("Condition") +
  theme(legend.position = "none") 


# Get ready to create difference scores and then merge with biomarkers later
colnames(HR_all_wide)[5] <- "Resting State"
HR_task_wide <- spread(HR_conditions_clean_tbl, Task_condition, mean_HR)
HR_final_wide <- merge(HR_all_wide, HR_task_wide, by = "ID", all = TRUE)
HR_final_wide_reduced <- HR_final_wide[,c("ID","Task","Resting State","Other","Self-Disconnected","Self-Connected","Rest")]
HR_final_wide_reduced$ID <- substr(HR_final_wide_reduced$ID, 1,6)
colnames(HR_final_wide_reduced)[1] <- "SampleID"
HR_final_wide_reduced$SampleID <- ifelse(HR_final_wide_reduced$SampleID == "SOS068",
                                         "TAG068", HR_final_wide_reduced$SampleID)
HR_final_wide_reduced$SampleID <- ifelse(HR_final_wide_reduced$SampleID == "SOS076",
                                         "TAG076", HR_final_wide_reduced$SampleID)

# Create the following HR difference scores:

# task heart rate reactivity variable: SOS task vs. resting state
HR_final_wide_reduced$task_rsHR_difference <- 
  HR_final_wide_reduced$Task - HR_final_wide_reduced$'Resting State' # SOS task overall - baseline resting state: if positive, it increased during the task

# “connected vs. disconnected heart rate reactivity”
HR_final_wide_reduced$conn_disHR_difference <- 
  HR_final_wide_reduced$'Self-Connected' - HR_final_wide_reduced$'Self-Disconnected'

# “connected vs. other heart rate reactivity”
HR_final_wide_reduced$conn_otherHR_difference <-
  HR_final_wide_reduced$'Self-Connected' - HR_final_wide_reduced$Other

# “connected vs. rest (fixation cross) heart rate reactivity”
HR_final_wide_reduced$conn_restHR_difference <-
  HR_final_wide_reduced$'Self-Connected' - HR_final_wide_reduced$Rest

HR_difference <- HR_final_wide_reduced[,c("SampleID","task_rsHR_difference",
                                          "conn_disHR_difference","conn_otherHR_difference",
                                          "conn_restHR_difference")]

summary(HR_difference)
sd(HR_difference$task_rsHR_difference, na.rm = TRUE)
sd(HR_difference$conn_disHR_difference, na.rm = TRUE)
sd(HR_difference$conn_otherHR_difference, na.rm = TRUE)
sd(HR_difference$conn_restHR_difference, na.rm = TRUE)


####### Salivary markers ####### 

# Report some basic plots taken from:
# S:\MNHS-Psych\ANDL-Lab-TAG-Study\SOS Study\saliva\Descriptives_plots
# and
# S:\MNHS-Psych\ANDL-Lab-TAG-Study\SOS Study\saliva\Biodata_scripts_output

# Also show plots for whole group with repeated measures (with standard errors)
# and growth curve modeling with quadratic curve (if it works)

sos_crp_panel <- panel_data(sos_crp[ , c("SampleID", "Time","CRP_log")], 
                            id = "SampleID", wave = "Time")
sos_cytokines_panel <- panel_data(sos_cytokines[ , c("SampleID", "Time","IL10_log","IL6_log","TNFalpha_log")], 
                                  id = "SampleID", wave = "Time")
sos_hormones_panel <- panel_data(sos_hormones[ , c("SampleID", "Time","DHEA_log","Testosterone_log","Cortisol_log")], 
                                 id = "SampleID", wave = "Time")

mlm_crp <- lmer(CRP_log ~ 1 + Time + (1 + Time | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_crp_panel)
summary(mlm_crp) # Time point does not linearly predict CRP (p = .173)
mlm2_crp <- lmer(CRP_log ~ 1 + poly(Time,2) + (1 + poly(Time,2) | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_crp_panel)
summary(mlm2_crp) # Time point does not quadratically predict CRP (p = .125)
AIC(mlm_crp,mlm2_crp)
BIC(mlm_crp,mlm2_crp)

mlm_IL10 <- lmer(IL10_log ~ 1 + Time + (1 + Time | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_cytokines_panel)
summary(mlm_IL10) # Time point *DOES* linearly predict IL10 (p = .027) but model might have singularity problem :|
mlm2_IL10 <- lmer(IL10_log ~ 1 + poly(Time,2) + (1 + poly(Time,2) | SampleID),
                 na.action = na.omit,
                 REML = TRUE,
                 data = sos_cytokines_panel)
summary(mlm2_IL10) # Time point *DOES* quadratically predict IL10 (p < .001)
AIC(mlm_IL10,mlm2_IL10)
BIC(mlm_IL10,mlm2_IL10)

mlm_IL6 <- lmer(IL6_log ~ 1 + Time + (1 + Time | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_cytokines_panel)
summary(mlm_IL6) # Time point does not linearly predict IL6 (p = .291)
mlm2_IL6 <- lmer(IL6_log ~ 1 + poly(Time,2) + (1 + poly(Time,2) | SampleID),
                 na.action = na.omit,
                 REML = TRUE,
                 data = sos_cytokines_panel)
summary(mlm2_IL6) # Time point *DOES* quadratically predict IL6 (p < .001)
AIC(mlm_IL6,mlm2_IL6)
BIC(mlm_IL6,mlm2_IL6)

mlm_TNFalpha <- lmer(TNFalpha_log ~ 1 + Time + (1 + Time | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_cytokines_panel)
summary(mlm_TNFalpha) # Time point does not linearly predict TNFa (p = .616)
mlm2_TNFalpha <- lmer(TNFalpha_log ~ 1 + poly(Time,2) + (1 + poly(Time,2) | SampleID),
                 na.action = na.omit,
                 REML = TRUE,
                 data = sos_cytokines_panel)
summary(mlm2_TNFalpha) # Time point *DOES* quadratically predict TNFa (p < .001)
AIC(mlm_TNFalpha,mlm2_TNFalpha)
BIC(mlm_TNFalpha,mlm2_TNFalpha)

mlm_dhea <- lmer(DHEA_log ~ 1 + Time + (1 + Time | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_hormones_panel)
summary(mlm_dhea) # Time point does not linearly predict DHEA (p = .070)
mlm2_dhea <- lmer(DHEA_log ~ 1 + poly(Time,2) + (1 + poly(Time,2) | SampleID),
                 na.action = na.omit,
                 REML = TRUE,
                 data = sos_hormones_panel)
summary(mlm2_dhea) # Time point *DOES* quadratically predict DHEA (p = .032)
AIC(mlm_dhea,mlm2_dhea)
BIC(mlm_dhea,mlm2_dhea)

mlm_test <- lmer(Testosterone_log ~ 1 + Time + (1 + Time | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_hormones_panel)
summary(mlm_test) # Time point does not linearly predict testosterone (p = .167)
mlm2_test <- lmer(Testosterone_log ~ 1 + poly(Time,2) + (1 + poly(Time,2) | SampleID),
                 na.action = na.omit,
                 REML = TRUE,
                 data = sos_hormones_panel)
summary(mlm2_test) # Time point *DOES* quadratically predict testosterone (p = .035)
AIC(mlm_test,mlm2_test)
BIC(mlm_test,mlm2_test)

mlm_cort <- lmer(Cortisol_log ~ 1 + Time + (1 + Time | SampleID),
                na.action = na.omit,
                REML = TRUE,
                data = sos_hormones_panel)
summary(mlm_cort) # Time point *DOES* linearly predict cortisol (p = .022)
mlm2_cort <- lmer(Cortisol_log ~ 1 + poly(Time,2) + (1 + poly(Time,2) | SampleID),
                 na.action = na.omit,
                 REML = TRUE,
                 data = sos_hormones_panel)
summary(mlm2_cort) # Time point *DOES* quadratically predict cortisol (p < .001)
AIC(mlm_cort,mlm2_cort)
BIC(mlm_cort,mlm2_cort)

# plot all the quad curves since they were da best

sos_cytokines_longest <- gather(sos_cytokines_panel, biomarker, 
                               log_conc, IL10_log:TNFalpha_log)
cytokines_quad_plot <- ggplot(sos_cytokines_longest, aes(Time, log_conc, color = biomarker)) +
#  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  scale_colour_viridis_d() +
  ylab("Log transformed concentration") +
  xlab("Time point") +
  ggtitle("Quadratic model for cytokines")
# Use this one above for the cytokines.

sos_hormones_longest <- gather(sos_hormones_panel, biomarker, 
                               log_conc, DHEA_log:Cortisol_log)
hormones_quad_plot <- ggplot(sos_hormones_longest, aes(Time, log_conc, color = biomarker)) +
#  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  scale_colour_viridis_d() +
  ylab("Log transformed concentration") +
  xlab("Time point") +
  ggtitle("Quadratic model for hormones")
# Use this one above for the hormones

ggarrange(cytokines_quad_plot, hormones_quad_plot, 
          labels = c("Cytokines","Hormones"),
          ncol = 2, nrow =1) # stretches weirdly

# The rest below are not as useful but just here in case:
ggplot(data = sos_crp_panel, aes(x = Time, y = CRP_log)) +
  #  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  scale_colour_viridis_d()
# Don't use the CRP one because it wasn't significant. It's just here in case.

immuneall <- merge(sos_crp_panel, sos_cytokines_panel, by = c("SampleID","Time"))
biomarkersall <- merge(immuneall, sos_hormones_panel, by = c("SampleID","Time"))

biomarkers_longest <- gather(biomarkersall, biomarker, 
                               log_conc, CRP_log:Cortisol_log)
ggplot(biomarkers_longest, aes(Time, log_conc, color = biomarker)) +
  #  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  scale_colour_viridis_d()

immune_longest <- gather(immuneall, biomarker, 
                             log_conc, CRP_log:TNFalpha_log)
ggplot(immune_longest, aes(Time, log_conc, color = biomarker)) +
  #  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  scale_colour_viridis_d()

## For each analyte, calculate each person's peak from T3 OR T4
# We do this because of the inconsistencies in the literature regarding when
# peak should happen - could be different by analyte and due to individual 
# differences. Doing this (picking T3 or T4) is basically a simpler version
# of that they do here (I think): https://pubmed.ncbi.nlm.nih.gov/24754834/
# Basically, what they did (they had heaps of time points) was manually select
# each person's peak, then they did a time-adjusted two-piece growth curve 
# model where the intercept was the person's peak. (They call this "landmark
# registration"*). Then the first part of the model was the
# activation slope (pre-peak slope) and the second part was the "regulation"
# slope (post-peak). I probably don't have enough time points for this to break
# into two parts. I think it's better to just select if T3 or T4 are that 
# person's peak. T2 was not designed to be the peak, so maybe for that we 
# just acknowledge some people had highest peak there due to MRI stress.
# Alternatively, could we include the T2 value as a covariate in further 
# analyses?
#
# * from the paper: "When alignment is based on a specific feature of the 
# curve, such as a response peak, the process is called landmark 
# registration (Kneip & Gasser, 1992)."

sos_crp_wide <- widen_panel(sos_crp_panel)
sos_crp_wide$CRPlog_t3t4peak <- ifelse(sos_crp_wide$CRP_log_3 > sos_crp_wide$CRP_log_4, 
                                       sos_crp_wide$CRP_log_3, sos_crp_wide$CRP_log_4)

sos_cytokines_wide <- widen_panel(sos_cytokines_panel)
sos_cytokines_wide$IL10log_t3t4peak <- ifelse(sos_cytokines_wide$IL10_log_3 > sos_cytokines_wide$IL10_log_4, 
                                              sos_cytokines_wide$IL10_log_3, sos_cytokines_wide$IL10_log_4)
sos_cytokines_wide$IL6log_t3t4peak <- ifelse(sos_cytokines_wide$IL6_log_3 > sos_cytokines_wide$IL6_log_4, 
                                              sos_cytokines_wide$IL6_log_3, sos_cytokines_wide$IL6_log_4)
sos_cytokines_wide$TNFalog_t3t4peak <- ifelse(sos_cytokines_wide$TNFalpha_log_3 > sos_cytokines_wide$TNFalpha_log_4, 
                                              sos_cytokines_wide$TNFalpha_log_3, sos_cytokines_wide$TNFalpha_log_4)

sos_hormones_wide <- widen_panel(sos_hormones_panel)
sos_hormones_wide$DHEAlog_t3t4peak <- ifelse(sos_hormones_wide$DHEA_log_3 > sos_hormones_wide$DHEA_log_4, 
                                              sos_hormones_wide$DHEA_log_3, sos_hormones_wide$DHEA_log_4)
sos_hormones_wide$Testlog_t3t4peak <- ifelse(sos_hormones_wide$Testosterone_log_3 > sos_hormones_wide$Testosterone_log_4, 
                                             sos_hormones_wide$Testosterone_log_3, sos_hormones_wide$Testosterone_log_4)
sos_hormones_wide$Cortlog_t3t4peak <- ifelse(sos_hormones_wide$Cortisol_log_3 > sos_hormones_wide$Cortisol_log_4, 
                                              sos_hormones_wide$Cortisol_log_3, sos_hormones_wide$Cortisol_log_4)

## Check if that "peak" is significantly higher than T1 levels
mean(sos_crp_wide$CRPlog_t3t4peak, na.rm = TRUE) # -2.052981
mean(sos_crp_wide$CRP_log_1, na.rm = TRUE) # -2.543955
t.test(sos_crp_wide$CRPlog_t3t4peak, sos_crp_wide$CRP_log_1, alternative = "greater", paired = TRUE)
# OR - should these be repeated measures ANOVAs to do all time points at once?

mean(sos_cytokines_wide$IL10log_t3t4peak, na.rm = TRUE) # 1.626901
mean(sos_cytokines_wide$IL10_log_1, na.rm = TRUE) # 1.089902
t.test(sos_cytokines_wide$IL10log_t3t4peak, sos_cytokines_wide$IL10_log_1, alternative = "greater", paired = TRUE)

mean(sos_cytokines_wide$IL6log_t3t4peak, na.rm = TRUE) # 1.01348
mean(sos_cytokines_wide$IL6_log_1, na.rm = TRUE) # 0.7275804
t.test(sos_cytokines_wide$IL6log_t3t4peak, sos_cytokines_wide$IL6_log_1, alternative = "greater", paired = TRUE)

mean(sos_cytokines_wide$TNFalog_t3t4peak, na.rm = TRUE) # 2.152535
mean(sos_cytokines_wide$TNFalpha_log_1, na.rm = TRUE) # 1.796177
t.test(sos_cytokines_wide$TNFalog_t3t4peak, sos_cytokines_wide$TNFalpha_log_1, alternative = "greater", paired = TRUE)

mean(sos_hormones_wide$DHEAlog_t3t4peak, na.rm = TRUE) # 4.621193
mean(sos_hormones_wide$DHEA_log_1, na.rm = TRUE) # 4.340595
t.test(sos_hormones_wide$DHEAlog_t3t4peak, sos_hormones_wide$DHEA_log_1, alternative = "greater", paired = TRUE)

mean(sos_hormones_wide$Testlog_t3t4peak, na.rm = TRUE) # 3.741303
mean(sos_hormones_wide$Testosterone_log_1, na.rm = TRUE) # 3.684813
t.test(sos_hormones_wide$Testlog_t3t4peak, sos_hormones_wide$Testosterone_log_1, alternative = "greater", paired = TRUE)

mean(sos_hormones_wide$Cortlog_t3t4peak, na.rm = TRUE) # -1.904597
mean(sos_hormones_wide$Cortisol_log_1, na.rm = TRUE) # -1.785707
t.test(sos_hormones_wide$Cortlog_t3t4peak, sos_hormones_wide$Cortisol_log_1, alternative = "greater", paired = TRUE)

# New long datasets with just T1 and the T3/T4 peak in it
sos_crp_long <- sos_crp_wide[,c("SampleID","CRP_log_1","CRPlog_t3t4peak")]
sos_crp_long <- pivot_longer(sos_crp_long, cols = c("CRP_log_1","CRPlog_t3t4peak"), 
                                            names_to = "time_point", 
                                            values_to = "CRP_log_conc")

sos_crp_long$time_point <- factor(sos_crp_long$time_point,
                                                 levels = c("CRP_log_1","CRPlog_t3t4peak"),
                                                 labels = c("T1","T3/T4 peak"))
sos_IL10_long <- sos_cytokines_wide[,c("SampleID","IL10_log_1","IL10log_t3t4peak")]
sos_IL10_long <- pivot_longer(sos_IL10_long, cols = c("IL10_log_1","IL10log_t3t4peak"), 
                             names_to = "time_point", 
                             values_to = "IL10_log_conc")

sos_IL10_long$time_point <- factor(sos_IL10_long$time_point,
                                  levels = c("IL10_log_1","IL10log_t3t4peak"),
                                  labels = c("T1","T3/T4 peak"))

sos_IL6_long <- sos_cytokines_wide[,c("SampleID","IL6_log_1","IL6log_t3t4peak")]
sos_IL6_long <- pivot_longer(sos_IL6_long, cols = c("IL6_log_1","IL6log_t3t4peak"), 
                             names_to = "time_point", 
                             values_to = "IL6_log_conc")

sos_IL6_long$time_point <- factor(sos_IL6_long$time_point,
                                  levels = c("IL6_log_1","IL6log_t3t4peak"),
                                  labels = c("T1","T3/T4 peak"))

sos_TNFa_long <- sos_cytokines_wide[,c("SampleID","TNFalpha_log_1","TNFalog_t3t4peak")]
sos_TNFa_long <- pivot_longer(sos_TNFa_long, cols = c("TNFalpha_log_1","TNFalog_t3t4peak"), 
                             names_to = "time_point", 
                             values_to = "TNFa_log_conc")

sos_TNFa_long$time_point <- factor(sos_TNFa_long$time_point,
                                  levels = c("TNFalpha_log_1","TNFalog_t3t4peak"),
                                  labels = c("T1","T3/T4 peak"))

sos_dhea_long <- sos_hormones_wide[,c("SampleID","DHEA_log_1","DHEAlog_t3t4peak")]
sos_dhea_long <- pivot_longer(sos_dhea_long, cols = c("DHEA_log_1","DHEAlog_t3t4peak"), 
                             names_to = "time_point", 
                             values_to = "DHEA_log_conc")

sos_dhea_long$time_point <- factor(sos_dhea_long$time_point,
                                  levels = c("DHEA_log_1","DHEAlog_t3t4peak"),
                                  labels = c("T1","T3/T4 peak"))

sos_test_long <- sos_hormones_wide[,c("SampleID","Testosterone_log_1","Testlog_t3t4peak")]
sos_test_long <- pivot_longer(sos_test_long, cols = c("Testosterone_log_1","Testlog_t3t4peak"), 
                              names_to = "time_point", 
                              values_to = "test_log_conc")

sos_test_long$time_point <- factor(sos_test_long$time_point,
                                   levels = c("Testosterone_log_1","Testlog_t3t4peak"),
                                   labels = c("T1","T3/T4 peak"))

sos_cort_long <- sos_hormones_wide[,c("SampleID","Cortisol_log_1","Cortlog_t3t4peak")]
sos_cort_long <- pivot_longer(sos_cort_long, cols = c("Cortisol_log_1","Cortlog_t3t4peak"), 
                              names_to = "time_point", 
                              values_to = "cort_log_conc")

sos_cort_long$time_point <- factor(sos_cort_long$time_point,
                                   levels = c("Cortisol_log_1","Cortlog_t3t4peak"),
                                   labels = c("T1","T3/T4 peak"))

# Set up comparisons to put means/medians on plots
my_comparisons_biomarkers <- list( c("T3/T4 peak","T1"))


get_box_stats_crp <- function(y, upper_limit = max(sos_crp_long$CRP_log_conc, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

get_box_stats_IL10 <- function(y, upper_limit = max(sos_IL10_long$IL10_log_conc, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

get_box_stats_IL6 <- function(y, upper_limit = max(sos_IL6_long$IL6_log_conc, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

get_box_stats_TNFa <- function(y, upper_limit = max(sos_TNFa_long$TNFa_log_conc, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

get_box_stats_dhea <- function(y, upper_limit = max(sos_dhea_long$DHEA_log_conc, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

get_box_stats_test <- function(y, upper_limit = max(sos_test_long$test_log_conc, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

get_box_stats_cort <- function(y, upper_limit = max(sos_cort_long$cort_log_conc, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n",
      "Median =", round(median(y), 2), "\n"
    )
  ))
}

CRP_plot <- 
  ggplot(sos_crp_long, aes(x = time_point, y = CRP_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_crp, geom = "text", hjust = 0.5, vjust = 4.6, size = 3.2, na.rm = TRUE) + # vjust = change to 8.3 if only using single plot and remove size argument
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_biomarkers, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_biomarkers order above
  ylab("Mean log(CRP)") +
  xlab("Time point") +
  theme(legend.position = "none") 

IL10_plot <- 
  ggplot(sos_IL10_long, aes(x = time_point, y = IL10_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_IL10, geom = "text", hjust = 0.5, vjust = 4.6, size = 3.2, na.rm = TRUE) + # vjust = change to 8.3 if only using single plot and remove size argument
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_biomarkers, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_biomarkers order above
  ylab("Mean log(IL-10)") +
  xlab("Time point") +
  theme(legend.position = "none") 

IL6_plot <- 
  ggplot(sos_IL6_long, aes(x = time_point, y = IL6_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_IL6, geom = "text", hjust = 0.5, vjust = 4.6, size = 3.2, na.rm = TRUE) + # vjust = change to 8.3 if only using single plot and remove size argument
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_biomarkers, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_biomarkers order above
  ylab("Mean log(IL-6)") +
  xlab("Time point") +
  theme(legend.position = "none") 

TNFa_plot <- 
  ggplot(sos_TNFa_long, aes(x = time_point, y = TNFa_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_TNFa, geom = "text", hjust = 0.5, vjust = 5, size = 3.2, na.rm = TRUE) + # vjust = change to 8.3 if only using single plot and remove size argument
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_biomarkers, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_biomarkers order above
  ylab("Mean log(TNF-alpha)") +
  xlab("Time point") +
  theme(legend.position = "none") 

dhea_plot <- 
  ggplot(sos_dhea_long, aes(x = time_point, y = DHEA_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_dhea, geom = "text", hjust = 0.5, vjust = 4.5, na.rm = TRUE) + 
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_biomarkers, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_biomarkers order above
  ylab("Mean log(TNF-alpha)") +
  xlab("Time point") +
  theme(legend.position = "none") 

test_plot <- 
  ggplot(sos_test_long, aes(x = time_point, y = test_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_test, geom = "text", hjust = 0.5, vjust = 4.5, na.rm = TRUE) + 
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_biomarkers, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_biomarkers order above
  ylab("Mean log(TNF-alpha)") +
  xlab("Time point") +
  theme(legend.position = "none") 

cort_plot <- 
  ggplot(sos_cort_long, aes(x = time_point, y = cort_log_conc)) +
  geom_violin(aes(fill = time_point), trim = FALSE) +
  geom_boxplot(width = 0.2) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_cort, geom = "text", hjust = 0.5, vjust = 4, na.rm = TRUE) + 
  stat_compare_means(method = "t.test", paired = TRUE, 
                     comparisons = my_comparisons_biomarkers, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_biomarkers order above
  ylab("Mean log(TNF-alpha)") +
  xlab("Time point") +
  theme(legend.position = "none") 

immune_plots <- ggarrange(CRP_plot, IL10_plot, IL6_plot, TNFa_plot, 
          labels = c("CRP","IL-10","IL-6","TNF-alpha"),
          ncol = 2, nrow =2) 

hormone_plots <- ggarrange(dhea_plot, test_plot, cort_plot, 
                          labels = c("DHEA","Testosterone","Cortisol"),
                          ncol = 3, nrow = 1) # these are stretched weirdly - export as height 300

# Correct for multiple comparisons using FDR within immune marker (4) and within hormones (3)
p_immune <- c(0.034, 0.000000013, 0.0013, 0.00066)
fdr_immune <- p.adjust(p_immune, method = "fdr", n = length(p_immune))

p_hormones <- c(0.016, 0.4, 0.86)
fdr_hormones <- p.adjust(p_hormones, method = "fdr", n = length(p_hormones))

## Acknowledge that x ppl (by analyte) had peak at T1, so they will have lower
# levels at T3/T4 than T1 - see if those people have lower subjective stress ratings
# This may be especially important for testosterone and cortisol

# Some people had peak at T2 which is due to scanner environment, but T3/T4 
# may still be higher than T1 for those people, so it's ok. But may want to
# check significant differences between T2 and T1 for everyone and report

# Then calculate each person's "reactivity" variable for use in further analyses
# This could be a simple difference score from T1 to T3/T4, or regress T1 on outcomes
# Could also calculate area under the curve, or slope from growth curve, but
# growth curve would have to have intercept set at the individual's peak 
# (landmark registration) and the slope would be the slope before that intercept
# See: https://pubmed.ncbi.nlm.nih.gov/24754834/

# Simple difference score between T1 and T3/T4 peak:
sos_crp_wide$crp_difference <- sos_crp_wide$CRPlog_t3t4peak - sos_crp_wide$CRP_log_1 # if they increased, this difference will be positive
crp_difference <- sos_crp_wide[,c("SampleID","crp_difference")]
sos_cytokines_wide$IL10_difference <- sos_cytokines_wide$IL10log_t3t4peak - sos_cytokines_wide$IL10_log_1
sos_cytokines_wide$IL6_difference <- sos_cytokines_wide$IL6log_t3t4peak - sos_cytokines_wide$IL6_log_1
sos_cytokines_wide$TNFa_difference <- sos_cytokines_wide$TNFalog_t3t4peak - sos_cytokines_wide$TNFalpha_log_1
cytokines_difference <- sos_cytokines_wide[,c("SampleID","IL10_difference","IL6_difference","TNFa_difference")]
sos_hormones_wide$dhea_difference <- sos_hormones_wide$DHEAlog_t3t4peak - sos_hormones_wide$DHEA_log_1
sos_hormones_wide$test_difference <- sos_hormones_wide$Testlog_t3t4peak - sos_hormones_wide$Testosterone_log_1
sos_hormones_wide$cort_difference <- sos_hormones_wide$Cortlog_t3t4peak - sos_hormones_wide$Cortisol_log_1
hormones_difference <- sos_hormones_wide[,c("SampleID","dhea_difference","test_difference","cort_difference")]
immune_difference <- merge(crp_difference,cytokines_difference, by = "SampleID")
biomarkers_difference <- merge(immune_difference, hormones_difference, by = "SampleID")
summary(biomarkers_difference)
sd(biomarkers_difference$crp_difference, na.rm = TRUE)
sd(biomarkers_difference$IL10_difference, na.rm = TRUE)
sd(biomarkers_difference$IL6_difference, na.rm = TRUE)
sd(biomarkers_difference$TNFa_difference, na.rm = TRUE)
sd(biomarkers_difference$dhea_difference, na.rm = TRUE)
sd(biomarkers_difference$test_difference, na.rm = TRUE)
sd(biomarkers_difference$cort_difference, na.rm = TRUE)


bio_all_difference <- merge(HR_difference, biomarkers_difference, by = "SampleID", all = TRUE)


##### Subjective ratings #####

sos_subj_stress <- read.csv(file.path(workdir,"/questionnaires/sos_subjective_stress_ratings.csv", fsep=""))
subj_stress_long <- pivot_longer(sos_subj_stress, 
                                 c("selfconnected_subj_stress","selfdisconnected_subj_stress","other_subj_stress"),
                                 names_to = "condition", values_to = "stress")

get_box_stats_subj <- function(y, upper_limit = max(sos_subj_stress$stress, na.rm = TRUE) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Mean =", round(mean(y), 2), "\n"
    )
  ))
}

subj_stress_long$condition <- factor(subj_stress_long$condition,
                                levels = c("selfconnected_subj_stress","selfdisconnected_subj_stress", "other_subj_stress"),
                                labels = c("Self-Connected","Self-Disconnected","Other"))


my_comparisons_subj <- list( c("Self-Connected","Self-Disconnected"),
                             c("Self-Connected","Other"),
                             c("Self-Disconnected","Other"))

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summary_subj <- summarySE(subj_stress_long, measurevar="stress", groupvars="condition")

subj_plot <- 
ggplot(summary_subj, aes(x=condition, y=stress, fill=condition)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=stress-se, ymax=stress+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  scale_fill_viridis_d() +
  stat_summary(fun.data = get_box_stats_subj, geom = "text", hjust = 0.5, vjust = .3, size = 3, na.rm = TRUE) + # vjust = change to 8.3 if only using single plot and remove size argument
  stat_compare_means(data = subj_stress_long, method = "wilcox.test", paired = TRUE, 
                     comparisons = my_comparisons_subj, 
                     method.args = list(alternative = "greater")) + # one tailed based on my_comparisons_subj order above
  ylab("Mean stress rating") +
  xlab("Condition") +
  theme(legend.position = "none") 

# Difference scores for subjective ratings
# SC vs SD
sos_subj_stress$subj_scsd_difference <- sos_subj_stress$selfconnected_subj_stress - sos_subj_stress$selfdisconnected_subj_stress # if they increased, this difference will be positive
# SC vs. Other
sos_subj_stress$subj_scother_difference <- sos_subj_stress$selfconnected_subj_stress - sos_subj_stress$other_subj_stress # if they increased, this difference will be positive

subjective_difference <- sos_subj_stress[,c("SampleID","subj_scsd_difference","subj_scother_difference")]

summary(subjective_difference)
sd(subjective_difference$subj_scsd_difference, na.rm = TRUE)
sd(subjective_difference$subj_scother_difference, na.rm = TRUE)

all_difference <- merge(subjective_difference, bio_all_difference, by = "SampleID", all = TRUE)


#### Correlations between reactivity (difference score) measures

difference_corrs = cor(as.matrix(all_difference[,2:14]), use = "pairwise.complete.obs")
difference_rcorrs = rcorr(as.matrix(all_difference[,2:14]))
difference_rcorrs$r   
difference_rcorrs$P

corrplot(difference_corrs, method = "number", tl.col="black", col=colorRampPalette(c("#FDE725FF","#29AF7FFF","#440154FF"))(100))


# SAVE WORKSPACE #

save.image('S:/MNHS-Psych/ANDL-Lab-TAG-Study/SOS Study/biopaper_analyses.RData')
         