setwd("S:/MNHS-Psych/ANDL-Lab-TAG-Study/SOS Study/physio/")

library(readxl)
library(tidyverse)

# load data
hr_data <- read_excel("Heart_Rate_Data_Task_Sequence_All_Slices_Conditions.xlsx")

# reformat & average data
hr_data_longer <- pivot_longer(hr_data, SOS068_20190303_150846:TAG223_20191116_145342, names_to = "ID", values_to = "HR")

hr_data_longer$Task_condition <- as.factor(hr_data_longer$Task_condition)

hr_data_longer$ID <- as.factor(hr_data_longer$ID)

hr_data_longer_avg <- hr_data_longer %>% group_by(ID, Task_condition) %>% summarize(mean_HR = mean(HR))

hr_data_clean <- hr_data_longer_avg %>% filter(Task_condition!="NA")
hr_data_clean$Task_condition <- factor(hr_data_clean$Task_condition)
hr_data_clean_tbl <- as_tibble(hr_data_clean)


# plot data 

p <- hr_data_clean_tbl %>% group_by(Task_condition) %>% summarize(m=mean(mean_HR), s=sd(mean_HR)) %>% 
  ggplot(aes(x=Task_condition, y=m)) + geom_bar(stat="identity", width=1, color='black', fill='#B71C1C', alpha=0.9, lwd=0.2)+
  geom_point(data=hr_data_clean_tbl, aes(x=Task_condition,y=mean_HR), alpha=0.5, size=1, fill='black', color='black')+
  geom_errorbar(aes(ymin=m-s, ymax=m+s),width=0.08)+
  ylab("Heart Rate")+xlab("")+theme_classic(base_size = 16)

p

library(ggpubr)
library(rstatix)
stat.test <- hr_data_clean_tbl %>% t_test(mean_HR ~ Task_condition, paired = TRUE) %>% add_significance() %>% add_xy_position(x = "Task_condition")
pp <-  p + stat_pvalue_manual(stat.test, label = "p.adj", bracket.nudge.y = 5, size=4, step.increase = 0.1) 
ggsave('Plot_HR_per_condition_pvals.png', pp, width=150, height=120, units=c("mm"), dpi=300)
