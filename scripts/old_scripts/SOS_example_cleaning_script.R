library(readxl)
library(dplyr)
library(stringr)
library(broom.helpers)
library(ggplot2)
library(psych)
library(rcompanion)
library(DescTools)
library(gridExtra)

# Import tab one of dataset without password
# This dataset has been manually edited to change TAG180 duplicates (x4) that were actually TAG081
# One TAG309 3i value was edited to be TAG 1i, with the lowest of duplicates considered baseline
# A TAG164 3i duplicate was also removed, with the value matching supplied GRIDs considered the 'real' value
# TAG007 4i was assayed twice in duplicate due to NaN values. Second assay without NaN taken (plate 5, cell 2)

CRP_SOS <- read_excel("human_cvd_panel_crp_raw_data.xlsx", 1)

summary(CRP_SOS$CRP_Ave)

# CRP assayed using Millipore Human Cardiovascular Disease (CVD) Magnetic Bead Panel 3 (Acute Phase) 96-Well Plate Assay Cat. # HCVD3MAG-67K
# Undetected values were replaced manually prior to dataset import with half the lower limit of detection

# 62 undetected CRP (in duplicates), replaced with half the lower limit of detection (i.e., 0.001/2 = 0.0005)
# Change all variables to numeric for downstream analysis

CRP_SOS <- CRP_SOS %>% 
  mutate_at(vars(CRP_Ave), as.numeric)

str(CRP_SOS)

# 78 unique participant IDs
length(unique(CRP_SOS$SampleID))

# How many participants have CRP at each time point?
sum(CRP_SOS$CRP_Ave != "NA" & CRP_SOS$Time == "1")
sum(CRP_SOS$CRP_Ave != "NA" & CRP_SOS$Time == "2")
sum(CRP_SOS$CRP_Ave != "NA" & CRP_SOS$Time == "3")
sum(CRP_SOS$CRP_Ave != "NA" & CRP_SOS$Time == "4")

# Calculate and report normality statistics (skew and kurtosis)
# Kurtosis and Skew should be -2/+2 (West, et al. 1995)
# All time points are non-normal
CRP_summary <- describeBy(CRP_SOS$CRP_Ave, group = CRP_SOS$Time)
CRP_summary

# If skew or kurtosis is < -2 or > 2 then transform
# Results were unacceptable for squareroot transformed data but acceptable after log transformation
# All time points now have skewness and kurtosis < +/-2
CRP_SOS$CRP_log <- log(CRP_SOS$CRP_Ave)
CRP_log_summary <- describeBy(CRP_SOS$CRP_log, group = CRP_SOS$Time)
CRP_log_summary

# Check how many outliers remaining after log transformation
CRPlog_up_limit1 <- (CRP_log_summary[[1]][["mean"]] + 3*(CRP_log_summary[[1]][["sd"]]))
CRPlog_lo_limit1 <- (CRP_log_summary[[1]][["mean"]] - 3*(CRP_log_summary[[1]][["sd"]]))

CRPlog_up_limit2 <- (CRP_log_summary[[2]][["mean"]] + 3*(CRP_log_summary[[2]][["sd"]]))
CRPlog_lo_limit2 <- (CRP_log_summary[[2]][["mean"]] - 3*(CRP_log_summary[[2]][["sd"]]))

CRPlog_up_limit3 <- (CRP_log_summary[[3]][["mean"]] + 3*(CRP_log_summary[[3]][["sd"]]))
CRPlog_lo_limit3 <- (CRP_log_summary[[3]][["mean"]] - 3*(CRP_log_summary[[3]][["sd"]]))

CRPlog_up_limit4 <- (CRP_log_summary[[4]][["mean"]] + 3*(CRP_log_summary[[4]][["sd"]]))
CRPlog_lo_limit4 <- (CRP_log_summary[[4]][["mean"]] - 3*(CRP_log_summary[[4]][["sd"]]))

# Check the number of outliers prior to winsorizing (where relevant)
# Time 1 zero upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 1 & CRP_SOS$CRP_log > CRPlog_up_limit1,na.rm=TRUE)
sum(CRP_SOS$Time == 1 & CRP_SOS$CRP_log < CRPlog_lo_limit1,na.rm=TRUE)

# Time 2 zero upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 2 & CRP_SOS$CRP_log > CRPlog_up_limit2,na.rm=TRUE)
sum(CRP_SOS$Time == 2 & CRP_SOS$CRP_log < CRPlog_lo_limit2,na.rm=TRUE)

# Time 3 zero upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 3 & CRP_SOS$CRP_log > CRPlog_up_limit3,na.rm=TRUE)
sum(CRP_SOS$Time == 3 & CRP_SOS$CRP_log < CRPlog_lo_limit3,na.rm=TRUE)

# Time 4 zero upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 4 & CRP_SOS$CRP_log > CRPlog_up_limit4,na.rm=TRUE)
sum(CRP_SOS$Time == 4 & CRP_SOS$CRP_log < CRPlog_lo_limit4,na.rm=TRUE)

# No need to winsorize CRP after log transformation, no outliers

# Michelle to discuss downstream options with collaborators, therefore log transformed and "raw + winsorized" to be exported
# Winsorize any variables that are +/-3 SDs from the mean (within time points)
# First calculate the +/-3 SDs from the mean at each time point for each variable
CRP_up_limit1 <- (CRP_summary[[1]][["mean"]] + 3*(CRP_summary[[1]][["sd"]]))
CRP_lo_limit1 <- (CRP_summary[[1]][["mean"]] - 3*(CRP_summary[[1]][["sd"]]))

CRP_up_limit2 <- (CRP_summary[[2]][["mean"]] + 3*(CRP_summary[[2]][["sd"]]))
CRP_lo_limit2 <- (CRP_summary[[2]][["mean"]] - 3*(CRP_summary[[2]][["sd"]]))

CRP_up_limit3 <- (CRP_summary[[3]][["mean"]] + 3*(CRP_summary[[3]][["sd"]]))
CRP_lo_limit3 <- (CRP_summary[[3]][["mean"]] - 3*(CRP_summary[[3]][["sd"]]))

CRP_up_limit4 <- (CRP_summary[[4]][["mean"]] + 3*(CRP_summary[[4]][["sd"]]))
CRP_lo_limit4 <- (CRP_summary[[4]][["mean"]] - 3*(CRP_summary[[4]][["sd"]]))

# Check the number of outliers prior to winsorizing (where relevant)
# Time 1 two upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 1 & CRP_SOS$CRP_Ave > CRP_up_limit1,na.rm=TRUE)
sum(CRP_SOS$Time == 1 & CRP_SOS$CRP_Ave < CRP_lo_limit1,na.rm=TRUE)

# Time 2 two upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 2 & CRP_SOS$CRP_Ave > CRP_up_limit2,na.rm=TRUE)
sum(CRP_SOS$Time == 2 & CRP_SOS$CRP_Ave < CRP_lo_limit2,na.rm=TRUE)

# Time 3 two upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 3 & CRP_SOS$CRP_Ave > CRP_up_limit3,na.rm=TRUE)
sum(CRP_SOS$Time == 3 & CRP_SOS$CRP_Ave < CRP_lo_limit3,na.rm=TRUE)

# Time 4 three upper limit outliers, zero lower limit
sum(CRP_SOS$Time == 4 & CRP_SOS$CRP_Ave > CRP_up_limit4,na.rm=TRUE)
sum(CRP_SOS$Time == 4 & CRP_SOS$CRP_Ave < CRP_lo_limit4,na.rm=TRUE)

# Winsorize by time point due to outliers, no need to winsorize lower limit as no outliers
CRP_SOS$CRP_nontran_win <- CRP_SOS$CRP_Ave

CRP_SOS$CRP_nontran_win <- ifelse(CRP_SOS$Time==1, (Winsorize(CRP_SOS$CRP_nontran_win, maxval = CRP_up_limit1)),CRP_SOS$CRP_nontran_win)
CRP_SOS$CRP_nontran_win <- ifelse(CRP_SOS$Time==2, (Winsorize(CRP_SOS$CRP_nontran_win, maxval = CRP_up_limit2)),CRP_SOS$CRP_nontran_win)
CRP_SOS$CRP_nontran_win <- ifelse(CRP_SOS$Time==3, (Winsorize(CRP_SOS$CRP_nontran_win, maxval = CRP_up_limit3)),CRP_SOS$CRP_nontran_win)
CRP_SOS$CRP_nontran_win <- ifelse(CRP_SOS$Time==4, (Winsorize(CRP_SOS$CRP_nontran_win, maxval = CRP_up_limit4)),CRP_SOS$CRP_nontran_win)

# Removing outliers on raw data does improve non-normality, but not entirely
CRP_nontran_win_summary <- describeBy(CRP_SOS$CRP_nontran_win, group = CRP_SOS$Time)

# Save the cleaned dataset
write.csv(CRP_SOS, file="SOS_CRP_cleaned.csv")

# Boxplots for CRP
CRP_rawboxplot <- ggplot(CRP_SOS, aes(x = "", y = CRP_Ave)) +
  geom_boxplot() +
  ylab("CRP (ng/mL)") +
  geom_smooth(method = 'lm', color = "black") + theme(axis.title.x=element_blank())

CRP_logboxplot <- ggplot(CRP_SOS, aes(x = "", y = CRP_log)) +   
  geom_boxplot() +
  ylab("CRP (log ng/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

CRP_winboxplot <- ggplot(CRP_SOS, aes(x = "", y = CRP_nontran_win)) +   
  geom_boxplot() +
  ylab("CRP win (ng/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

# Histograms of non-transformed and log-transformed variable
CRP_hist <-ggplot(CRP_SOS, aes(x=CRP_Ave)) + geom_histogram(color="black", fill = "violetred4", binwidth = 1.0) + theme_bw()
CRP_loghist <-ggplot(CRP_SOS, aes(x=CRP_log)) + geom_histogram(color="black", fill = "violetred4", binwidth = 0.50) + theme_bw()
CRP_winhist <-ggplot(CRP_SOS, aes(x=CRP_nontran_win)) + geom_histogram(color="black", fill = "violetred4", binwidth = 0.50) + theme_bw()

grid.arrange(CRP_rawboxplot, CRP_winboxplot, CRP_logboxplot, CRP_hist, CRP_winhist, CRP_loghist, ncol=3)


################################################################################


# Import the first excel tab of TAG dataset without password
# This dataset has been manually edited to change TAG180 duplicates (x4) that were actually TAG081
# One TAG309 3i value was edited to be TAG 1i, with the lowest of duplicates considered baseline
# A TAG164 3i duplicate was also removed, with the value matching supplied GRIDs included as the 'real value'

# Undetectables below the assay lower limits were manually replaced prior to average calculation
# per Millipore Human Cytokine/Chemokine Magnetic Bead Panel 96 Well Plate Assay
# 63 undetected IL-10 (in duplicates), replaced with half the lower limit of detection (i.e., 1.1/2 = 0.55)
# 106 undetected IL-6 (in duplicates), replaced with half the lower limit of detection (i.e., 0.9/2 = 0.45)
# 15 undetected TNF-alpha (in duplicates), replaced with half the lower limit of detection (i.e., 0.7/2 = 0.35)

cytokines_SOS <- read_excel("Human_cytokine_rawdata.xlsx", 1)

# Change all variables to numeric for downstream analysis
cytokines_SOS <- cytokines_SOS %>% 
  mutate_at(vars(IL10_Ave, IL6_Ave, TNFalpha_Ave), as.numeric)

str(cytokines_SOS)

# Number of unique participant IDs
length(unique(cytokines_SOS$SampleID))

# How many participants have IL-10 at each time point?
sum(cytokines_SOS$IL10_Ave != "NA" & cytokines_SOS$Time == "1")
sum(cytokines_SOS$IL10_Ave != "NA" & cytokines_SOS$Time == "2")
sum(cytokines_SOS$IL10_Ave != "NA" & cytokines_SOS$Time == "3")
sum(cytokines_SOS$IL10_Ave != "NA" & cytokines_SOS$Time == "4")

# How many participants have IL-6 at each time point?
sum(cytokines_SOS$IL6_Ave != "NA" & cytokines_SOS$Time == "1")
sum(cytokines_SOS$IL6_Ave != "NA" & cytokines_SOS$Time == "2")
sum(cytokines_SOS$IL6_Ave != "NA" & cytokines_SOS$Time == "3")
sum(cytokines_SOS$IL6_Ave != "NA" & cytokines_SOS$Time == "4")

# How many participants have TNF-alpha at each time point?
sum(cytokines_SOS$TNFalpha_Ave != "NA" & cytokines_SOS$Time == "1")
sum(cytokines_SOS$TNFalpha_Ave != "NA" & cytokines_SOS$Time == "2")
sum(cytokines_SOS$TNFalpha_Ave != "NA" & cytokines_SOS$Time == "3")
sum(cytokines_SOS$TNFalpha_Ave != "NA" & cytokines_SOS$Time == "4")

# Calculate and report normality statistics (skew and kurtosis)
# Kurtosis and Skew should be -2/+2 (West, et al. 1995)
# Most time points are non-normal
IL10_summary <- describeBy(cytokines_SOS$IL10_Ave, group = cytokines_SOS$Time)
IL10_summary

IL6_summary <- describeBy(cytokines_SOS$IL6_Ave, group = cytokines_SOS$Time)
IL6_summary

TNFalpha_summary <- describeBy(cytokines_SOS$TNFalpha_Ave, group = cytokines_SOS$Time)
TNFalpha_summary

# Boxplots for the cytokines
ggplot(cytokines_SOS, aes(x = "", y = IL10_Ave)) +   
  geom_boxplot() +
  ylab("IL-10 (pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

ggplot(cytokines_SOS, aes(x = "", y = IL6_Ave)) +   
  geom_boxplot() +
  ylab("IL-6 (pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

ggplot(cytokines_SOS, aes(x = "", y = TNFalpha_Ave)) +   
  geom_boxplot() +
  ylab("TNF-alpha (pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

# Calculate the log of each variable
cytokines_SOS$IL10_log <- log(cytokines_SOS$IL10_Ave)
cytokines_SOS$IL6_log <- log(cytokines_SOS$IL6_Ave)
cytokines_SOS$TNFalpha_log <- log(cytokines_SOS$TNFalpha_Ave)

# Calculate and report normality statistics (skew and kurtosis)
# All cytokines are now acceptable +/- 2 skew and kurtosis

IL10_summary_log <- describeBy(cytokines_SOS$IL10_log, group = cytokines_SOS$Time)
IL10_summary_log

IL6_summary_log <- describeBy(cytokines_SOS$IL6_log, group = cytokines_SOS$Time)
IL6_summary_log

TNFalpha_summary_log <- describeBy(cytokines_SOS$TNFalpha_log, group = cytokines_SOS$Time)
TNFalpha_summary_log

# IL-10
IL10log_up_limit1 <- (IL10_summary_log[[1]][["mean"]] + 3*(IL10_summary_log[[1]][["sd"]]))
IL10log_lo_limit1 <- (IL10_summary_log[[1]][["mean"]] - 3*(IL10_summary_log[[1]][["sd"]]))

IL10log_up_limit2 <- (IL10_summary_log[[2]][["mean"]] + 3*(IL10_summary_log[[2]][["sd"]]))
IL10log_lo_limit2 <- (IL10_summary_log[[2]][["mean"]] - 3*(IL10_summary_log[[2]][["sd"]]))

IL10log_up_limit3 <- (IL10_summary_log[[3]][["mean"]] + 3*(IL10_summary_log[[3]][["sd"]]))
IL10log_lo_limit3 <- (IL10_summary_log[[3]][["mean"]] - 3*(IL10_summary_log[[3]][["sd"]]))

IL10log_up_limit4 <- (IL10_summary_log[[4]][["mean"]] + 3*(IL10_summary_log[[4]][["sd"]]))
IL10log_lo_limit4 <- (IL10_summary_log[[4]][["mean"]] - 3*(IL10_summary_log[[4]][["sd"]]))

# IL-6
IL6log_up_limit1 <- (IL6_summary_log[[1]][["mean"]] + 3*(IL6_summary_log[[1]][["sd"]]))
IL6log_lo_limit1 <- (IL6_summary_log[[1]][["mean"]] - 3*(IL6_summary_log[[1]][["sd"]]))

IL6log_up_limit2 <- (IL6_summary_log[[2]][["mean"]] + 3*(IL6_summary_log[[2]][["sd"]]))
IL6log_lo_limit2 <- (IL6_summary_log[[2]][["mean"]] - 3*(IL6_summary_log[[2]][["sd"]]))

IL6log_up_limit3 <- (IL6_summary_log[[3]][["mean"]] + 3*(IL6_summary_log[[3]][["sd"]]))
IL6log_lo_limit3 <- (IL6_summary_log[[3]][["mean"]] - 3*(IL6_summary_log[[3]][["sd"]]))

IL6log_up_limit4 <- (IL6_summary_log[[4]][["mean"]] + 3*(IL6_summary_log[[4]][["sd"]]))
IL6log_lo_limit4 <- (IL6_summary_log[[4]][["mean"]] - 3*(IL6_summary_log[[4]][["sd"]]))

# TNF_alpha
TNFalphalog_up_limit1 <- (TNFalpha_summary_log[[1]][["mean"]] + 3*(TNFalpha_summary_log[[1]][["sd"]]))
TNFalphalog_lo_limit1 <- (TNFalpha_summary_log[[1]][["mean"]] - 3*(TNFalpha_summary_log[[1]][["sd"]]))

TNFalphalog_up_limit2 <- (TNFalpha_summary_log[[2]][["mean"]] + 3*(TNFalpha_summary_log[[2]][["sd"]]))
TNFalphalog_lo_limit2 <- (TNFalpha_summary_log[[2]][["mean"]] - 3*(TNFalpha_summary_log[[2]][["sd"]]))

TNFalphalog_up_limit3 <- (TNFalpha_summary_log[[3]][["mean"]] + 3*(TNFalpha_summary_log[[3]][["sd"]]))
TNFalphalog_lo_limit3 <- (TNFalpha_summary_log[[3]][["mean"]] - 3*(TNFalpha_summary_log[[3]][["sd"]]))

TNFalphalog_up_limit4 <- (TNFalpha_summary_log[[4]][["mean"]] + 3*(TNFalpha_summary_log[[4]][["sd"]]))
TNFalphalog_lo_limit4 <- (TNFalpha_summary_log[[4]][["mean"]] - 3*(TNFalpha_summary_log[[4]][["sd"]]))

# IL-10
# Check the number of outliers prior to winsorizing (where relevant)
# Time 1, zero upper limit outliers, zero lower limit. 
# Winsorization of raw data therefore loses a lot of variation
sum(cytokines_SOS$Time == 1 & cytokines_SOS$IL10_log > IL10log_up_limit1,na.rm=TRUE)
sum(cytokines_SOS$Time == 1 & cytokines_SOS$IL10_log < IL10log_lo_limit1,na.rm=TRUE)

# Time 2, zero upper limit outliers, zero lower limit
sum(cytokines_SOS$Time == 2 & cytokines_SOS$IL10_log > IL10log_up_limit2,na.rm=TRUE)
sum(cytokines_SOS$Time == 2 & cytokines_SOS$IL10_log < IL10log_lo_limit2,na.rm=TRUE)

# Time 3, zero upper limit outliers, zero lower limit
sum(cytokines_SOS$Time == 3 & cytokines_SOS$IL10_log > IL10log_up_limit3,na.rm=TRUE)
sum(cytokines_SOS$Time == 3 & cytokines_SOS$IL10_log < IL10log_lo_limit3,na.rm=TRUE)

# Time 4, zero upper limit outliers, zero lower limit
sum(cytokines_SOS$Time == 4 & cytokines_SOS$IL10_log > IL10log_up_limit4,na.rm=TRUE)
sum(cytokines_SOS$Time == 4 & cytokines_SOS$IL10_log < IL10log_lo_limit4,na.rm=TRUE)

# Boxplots for the cytokines
ggplot(cytokines_SOS, aes(x = "", y = IL10_log)) +   
  geom_boxplot() +
  ylab("log IL-10 (pg/mL)") +
  ggtitle("SOS IL-10 log transformed") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

ggplot(cytokines_SOS, aes(x = "", y = IL6_log)) +   
  geom_boxplot() +
  ylab("log IL-6 (pg/mL)") +
  ggtitle("SOS IL-6 log transformed") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

ggplot(cytokines_SOS, aes(x = "", y = TNFalpha_log)) +   
  geom_boxplot() +
  ylab("log TNF-alpha (pg/mL)") +
  ggtitle("SOS TNF-alpha log transformed") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

# Histograms of transformed variables
ggplot(cytokines_SOS, aes(x=IL10_log)) + geom_histogram(color="black", fill = "violetred4", binwidth = 0.50) + theme_bw()
ggplot(cytokines_SOS, aes(x=IL6_log)) + geom_histogram(color="black", fill = "violetred4", binwidth = 0.25) + theme_bw()
ggplot(cytokines_SOS, aes(x=TNFalpha_log)) + geom_histogram(color="black", fill = "violetred4", binwidth = 0.20) + theme_bw()

# Winsorize any observations that are +/-3 SDs (IL-10, IL-6, TNF-alpha)
# First calculate the +/-3 SDs from the mean for each variable
IL10_up_limit <- (IL10_summary$mean + 3*(IL10_summary$sd))
IL6_up_limit <- (IL6_summary$mean + 3*(IL6_summary$sd))
TNFalpha_up_limit <- (TNFalpha_summary$mean + 3*(TNFalpha_summary$sd))

IL10_lo_limit <- (IL10_summary$mean - 3*(IL10_summary$sd))
IL6_lo_limit <- (IL6_summary$mean - 3*(IL6_summary$sd))
TNFalpha_lo_limit <- (TNFalpha_summary$mean - 3*(TNFalpha_summary$sd))

# Michelle to discuss downstream options with collaborators, therefore log transformed and "raw" winsorized to be exported
# Winsorize any variables that are +/-3 SDs from the mean (within time points)
# First calculate the +/-3 SDs from the mean for each variable
# IL-10
IL10_up_limit1 <- (IL10_summary[[1]][["mean"]] + 3*(IL10_summary[[1]][["sd"]]))
IL10_lo_limit1 <- (IL10_summary[[1]][["mean"]] - 3*(IL10_summary[[1]][["sd"]]))

IL10_up_limit2 <- (IL10_summary[[2]][["mean"]] + 3*(IL10_summary[[2]][["sd"]]))
IL10_lo_limit2 <- (IL10_summary[[2]][["mean"]] - 3*(IL10_summary[[2]][["sd"]]))

IL10_up_limit3 <- (IL10_summary[[3]][["mean"]] + 3*(IL10_summary[[3]][["sd"]]))
IL10_lo_limit3 <- (IL10_summary[[3]][["mean"]] - 3*(IL10_summary[[3]][["sd"]]))

IL10_up_limit4 <- (IL10_summary[[4]][["mean"]] + 3*(IL10_summary[[4]][["sd"]]))
IL10_lo_limit4 <- (IL10_summary[[4]][["mean"]] - 3*(IL10_summary[[4]][["sd"]]))

# IL-6
IL6_up_limit1 <- (IL6_summary[[1]][["mean"]] + 3*(IL6_summary[[1]][["sd"]]))
IL6_lo_limit1 <- (IL6_summary[[1]][["mean"]] - 3*(IL6_summary[[1]][["sd"]]))

IL6_up_limit2 <- (IL6_summary[[2]][["mean"]] + 3*(IL6_summary[[2]][["sd"]]))
IL6_lo_limit2 <- (IL6_summary[[2]][["mean"]] - 3*(IL6_summary[[2]][["sd"]]))

IL6_up_limit3 <- (IL6_summary[[3]][["mean"]] + 3*(IL6_summary[[3]][["sd"]]))
IL6_lo_limit3 <- (IL6_summary[[3]][["mean"]] - 3*(IL6_summary[[3]][["sd"]]))

IL6_up_limit4 <- (IL6_summary[[4]][["mean"]] + 3*(IL6_summary[[4]][["sd"]]))
IL6_lo_limit4 <- (IL6_summary[[4]][["mean"]] - 3*(IL6_summary[[4]][["sd"]]))

# TNF_alpha
TNFalpha_up_limit1 <- (TNFalpha_summary[[1]][["mean"]] + 3*(TNFalpha_summary[[1]][["sd"]]))
TNFalpha_lo_limit1 <- (TNFalpha_summary[[1]][["mean"]] - 3*(TNFalpha_summary[[1]][["sd"]]))

TNFalpha_up_limit2 <- (TNFalpha_summary[[2]][["mean"]] + 3*(TNFalpha_summary[[2]][["sd"]]))
TNFalpha_lo_limit2 <- (TNFalpha_summary[[2]][["mean"]] - 3*(TNFalpha_summary[[2]][["sd"]]))

TNFalpha_up_limit3 <- (TNFalpha_summary[[3]][["mean"]] + 3*(TNFalpha_summary[[3]][["sd"]]))
TNFalpha_lo_limit3 <- (TNFalpha_summary[[3]][["mean"]] - 3*(TNFalpha_summary[[3]][["sd"]]))

TNFalpha_up_limit4 <- (TNFalpha_summary[[4]][["mean"]] + 3*(TNFalpha_summary[[4]][["sd"]]))
TNFalpha_lo_limit4 <- (TNFalpha_summary[[4]][["mean"]] - 3*(TNFalpha_summary[[4]][["sd"]]))

# IL-10
# Check the number of outliers prior to winsorizing (where relevant)
# Time 1, 1 upper limit outliers, zero lower limit. 
# Winsorization of raw data therefore loses a lot of variation
sum(cytokines_SOS$Time == 1 & cytokines_SOS$IL10_Ave > IL10_up_limit1,na.rm=TRUE)
sum(cytokines_SOS$Time == 1 & cytokines_SOS$IL10_Ave < IL10_lo_limit1,na.rm=TRUE)

# Time 2 1 upper limit outliers, zero lower limit
sum(cytokines_SOS$Time == 2 & cytokines_SOS$IL10_Ave > IL10_up_limit2,na.rm=TRUE)
sum(cytokines_SOS$Time == 2 & cytokines_SOS$IL10_Ave < IL10_lo_limit2,na.rm=TRUE)

# Time 3 1 upper limit outliers, zero lower limit
sum(cytokines_SOS$Time == 3 & cytokines_SOS$IL10_Ave > IL10_up_limit3,na.rm=TRUE)
sum(cytokines_SOS$Time == 3 & cytokines_SOS$IL10_Ave < IL10_lo_limit3,na.rm=TRUE)

# Time 4 1 upper limit outliers, zero lower limit
sum(cytokines_SOS$Time == 4 & cytokines_SOS$IL10_Ave > IL10_up_limit4,na.rm=TRUE)
sum(cytokines_SOS$Time == 4 & cytokines_SOS$IL10_Ave < IL10_lo_limit4,na.rm=TRUE)

# Winsorize by time point due to outliers, no need to winsorize lower limit as no outliers
# IL-10
cytokines_SOS$IL10_nontran_win <- cytokines_SOS$IL10_Ave

cytokines_SOS$IL10_nontran_win  <- ifelse(cytokines_SOS$Time==1, (Winsorize(cytokines_SOS$IL10_nontran_win, maxval = IL10_up_limit1)),cytokines_SOS$IL10_nontran_win)
cytokines_SOS$IL10_nontran_win  <- ifelse(cytokines_SOS$Time==2, (Winsorize(cytokines_SOS$IL10_nontran_win, maxval = IL10_up_limit2)),cytokines_SOS$IL10_nontran_win)
cytokines_SOS$IL10_nontran_win  <- ifelse(cytokines_SOS$Time==3, (Winsorize(cytokines_SOS$IL10_nontran_win, maxval = IL10_up_limit3)),cytokines_SOS$IL10_nontran_win)
cytokines_SOS$IL10_nontran_win  <- ifelse(cytokines_SOS$Time==4, (Winsorize(cytokines_SOS$IL10_nontran_win, maxval = IL10_up_limit4)),cytokines_SOS$IL10_nontran_win)

# IL-6
cytokines_SOS$IL6_nontran_win <- cytokines_SOS$IL6_Ave

cytokines_SOS$IL6_nontran_win  <- ifelse(cytokines_SOS$Time==1, (Winsorize(cytokines_SOS$IL6_nontran_win, maxval = IL6_up_limit1)),cytokines_SOS$IL6_nontran_win)
cytokines_SOS$IL6_nontran_win  <- ifelse(cytokines_SOS$Time==2, (Winsorize(cytokines_SOS$IL6_nontran_win, maxval = IL6_up_limit2)),cytokines_SOS$IL6_nontran_win)
cytokines_SOS$IL6_nontran_win  <- ifelse(cytokines_SOS$Time==3, (Winsorize(cytokines_SOS$IL6_nontran_win, maxval = IL6_up_limit3)),cytokines_SOS$IL6_nontran_win)
cytokines_SOS$IL6_nontran_win  <- ifelse(cytokines_SOS$Time==4, (Winsorize(cytokines_SOS$IL6_nontran_win, maxval = IL6_up_limit4)),cytokines_SOS$IL6_nontran_win)

# TNF_alpha
cytokines_SOS$TNFalpha_nontran_win <- cytokines_SOS$TNFalpha_Ave

cytokines_SOS$TNFalpha_nontran_win  <- ifelse(cytokines_SOS$Time==1, (Winsorize(cytokines_SOS$TNFalpha_nontran_win, maxval = TNFalpha_up_limit1)),cytokines_SOS$TNFalpha_nontran_win)
cytokines_SOS$TNFalpha_nontran_win  <- ifelse(cytokines_SOS$Time==2, (Winsorize(cytokines_SOS$TNFalpha_nontran_win, maxval = TNFalpha_up_limit2)),cytokines_SOS$TNFalpha_nontran_win)
cytokines_SOS$TNFalpha_nontran_win  <- ifelse(cytokines_SOS$Time==3, (Winsorize(cytokines_SOS$TNFalpha_nontran_win, maxval = TNFalpha_up_limit3)),cytokines_SOS$TNFalpha_nontran_win)
cytokines_SOS$TNFalpha_nontran_win  <- ifelse(cytokines_SOS$Time==4, (Winsorize(cytokines_SOS$TNFalpha_nontran_win, maxval = TNFalpha_up_limit4)),cytokines_SOS$TNFalpha_nontran_win)


# Removing outliers on raw data does improve non-normality, but not entirely
IL6_nontran_win_summary <- describeBy(cytokines_SOS$IL6_nontran_win, group = cytokines_SOS$Time)
IL10_nontran_win_summary <- describeBy(cytokines_SOS$IL10_nontran_win, group = cytokines_SOS$Time)
TNFalpha_nontran_win_summary <- describeBy(cytokines_SOS$TNFalpha_nontran_win, group = cytokines_SOS$Time)

# Not necessary to winsorize as no outliers
# Save the cleaned dataset
write.csv(cytokines_SOS, file="SOS_cytokines_cleaned.csv")


# Boxplots for IL6
IL6_rawboxplot <- ggplot(cytokines_SOS, aes(x = "", y = IL6_Ave)) +
  geom_boxplot() +
  ylab("IL6 (pg/mL)") +
  geom_smooth(method = 'lm', color = "black") + theme(axis.title.x=element_blank())

IL6_logboxplot <- ggplot(cytokines_SOS, aes(x = "", y = IL6_log)) +   
  geom_boxplot() +
  ylab("IL6 (log pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

IL6_winboxplot <- ggplot(cytokines_SOS, aes(x = "", y = IL6_nontran_win)) +   
  geom_boxplot() +
  ylab("IL6 win (pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

# IL6 histograms of non-transformed and log-transformed variable
IL6_hist <-ggplot(cytokines_SOS, aes(x=IL6_Ave)) + geom_histogram(color="black", fill = "violetred4", binwidth = 5.0) + theme_bw()
IL6_loghist <-ggplot(cytokines_SOS, aes(x=IL6_log)) + geom_histogram(color="black", fill = "violetred4", binwidth = 0.50) + theme_bw()
IL6_winhist <-ggplot(cytokines_SOS, aes(x=IL6_nontran_win)) + geom_histogram(color="black", fill = "violetred4", binwidth = 2.5) + theme_bw()

grid.arrange(IL6_rawboxplot, IL6_winboxplot, IL6_logboxplot, IL6_hist, IL6_winhist, IL6_loghist, ncol=3)

# Boxplots for IL10
IL10_rawboxplot <- ggplot(cytokines_SOS, aes(x = "", y = IL10_Ave)) +
  geom_boxplot() +
  ylab("IL10 (pg/mL)") +
  geom_smooth(method = 'lm', color = "black") + theme(axis.title.x=element_blank())

IL10_logboxplot <- ggplot(cytokines_SOS, aes(x = "", y = IL10_log)) +   
  geom_boxplot() +
  ylab("IL10 (log pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

IL10_winboxplot <- ggplot(cytokines_SOS, aes(x = "", y = IL10_nontran_win)) +   
  geom_boxplot() +
  ylab("IL10 win (pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())


# IL10 histograms of non-transformed and log-transformed variable
IL10_hist <-ggplot(cytokines_SOS, aes(x=IL10_Ave)) + geom_histogram(color="black", fill = "violetred4", binwidth = 5.0) + theme_bw()
IL10_loghist <-ggplot(cytokines_SOS, aes(x=IL10_log)) + geom_histogram(color="black", fill = "violetred4", binwidth = 0.50) + theme_bw()
IL10_winhist <-ggplot(cytokines_SOS, aes(x=IL10_nontran_win)) + geom_histogram(color="black", fill = "violetred4", binwidth = 2.5) + theme_bw()

grid.arrange(IL10_rawboxplot, IL10_winboxplot, IL10_logboxplot, IL10_hist, IL10_winhist, IL10_loghist, ncol=3)


# Boxplots for TNFalpha
TNFalpha_rawboxplot <- ggplot(cytokines_SOS, aes(x = "", y = TNFalpha_Ave)) +
  geom_boxplot() +
  ylab("TNFalpha (pg/mL)") +
  geom_smooth(method = 'lm', color = "black") + theme(axis.title.x=element_blank())

TNFalpha_logboxplot <- ggplot(cytokines_SOS, aes(x = "", y = TNFalpha_log)) +   
  geom_boxplot() +
  ylab("TNFalpha (log pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())

TNFalpha_winboxplot <- ggplot(cytokines_SOS, aes(x = "", y = TNFalpha_nontran_win)) +   
  geom_boxplot() +
  ylab("TNFalpha win (pg/mL)") +
  geom_smooth(method='lm', color="black") + theme(axis.title.x=element_blank())


# TNFalpha histograms of non-transformed and log-transformed variable
TNFalpha_hist <-ggplot(cytokines_SOS, aes(x=TNFalpha_Ave)) + geom_histogram(color="black", fill = "violetred4", binwidth = 10) + theme_bw()
TNFalpha_loghist <-ggplot(cytokines_SOS, aes(x=TNFalpha_log)) + geom_histogram(color="black", fill = "violetred4", binwidth = .50) + theme_bw()
TNFalpha_winhist <-ggplot(cytokines_SOS, aes(x=TNFalpha_nontran_win)) + geom_histogram(color="black", fill = "violetred4", binwidth = 5.0) + theme_bw()

grid.arrange(TNFalpha_rawboxplot, TNFalpha_winboxplot, TNFalpha_logboxplot, TNFalpha_hist, TNFalpha_winhist, TNFalpha_loghist, ncol=3)