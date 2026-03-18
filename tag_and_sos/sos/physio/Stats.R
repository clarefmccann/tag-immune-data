setwd('/home/clairea/wa26/sos/Heart_rate_analysis/Step3_Run_PhysIO_Toolbox/')

library(tidyverse)

# Import/organise/clean data

HR_Resting1 <- read_csv('Resting1/Mean_HR_Resting1.csv')
HR_Resting2 <- read_csv('Resting2/Mean_HR_Resting2.csv')
HR_Task <- read_csv('Task/Mean_HR_Task.csv')

colnames(HR_Resting1)[1] <- "ID"
colnames(HR_Resting2)[1] <- "ID"
colnames(HR_Task)[1] <- "ID"

# We couldn't get HR data for Resting2 for TAG 080. To make things easier, delete that subject across all sequences:
HR_Resting1 <- filter(HR_Resting1, ID!="TAG080_20190329_150631")
HR_Task <- filter(HR_Task, ID!="TAG080_20190329_150631")

# get rid of the repeat 313
HR_Resting1 <- filter(HR_Resting1, ID!="TAG313p2_20191117_151413")
HR_Resting2 <- filter(HR_Resting2, ID!="TAG313p2_20191117_151413")
HR_Task <- filter(HR_Task, ID!="TAG313p2_20191117_151413")

# get rid of the low outliers (task)
#HR_Resting1 <- filter(HR_Resting1, ID!="TAG077_20190823_145422")
#HR_Resting2 <- filter(HR_Resting2, ID!="TAG077_20190823_145422")
#HR_Task <- filter(HR_Task, ID!="TAG077_20190823_145422")

#HR_Resting1 <- filter(HR_Resting1, ID!="TAG188_20190731_145115")
#HR_Resting2 <- filter(HR_Resting2, ID!="TAG188_20190731_145115")
#HR_Task <- filter(HR_Task, ID!="TAG188_20190731_145115")

HR_Resting1$sequence <- rep(c("Resting 1"), each=51)
HR_Resting2$sequence <- rep(c("Resting state"), each=51)
HR_Task$sequence <- rep(c("Task"), each=51)

HR_all_sequences <- bind_rows(HR_Resting1, HR_Resting2, HR_Task)

# just resting2 and task
HR_rs2_task <- bind_rows(HR_Resting2, HR_Task)


# plotting 

library(ggplot2)

p <- HR_all_sequences %>% 
  group_by(sequence) %>% 
  summarize(m = mean(Mean_HR),
            s = sd(Mean_HR)) %>% 
  ggplot(aes(x = sequence, y = m)) + 
  geom_point(data=HR_all_sequences, aes(x=sequence,y=Mean_HR), alpha=0.5, size=2, fill='black', color='black')+
  geom_bar(stat = "identity", width=1, color='black', fill='#8A3324', alpha=0.5, lwd=0.2)+
  geom_errorbar(aes(ymin = m-s, ymax = m+s), width=0.2) +
  ylab("Mean Heart Rate")+xlab("")+
  theme_classic(base_size=16) 

ggsave('Plot_all_sequences.png', p, width=120, height=100, units=c("mm"), dpi=300)

 p <-   HR_rs2_task %>% 
  group_by(sequence) %>% 
  summarize(m = mean(Mean_HR),
            s = sd(Mean_HR)) %>% 
  ggplot(aes(x = sequence, y = m)) +
  geom_bar(stat = "identity", width=1, color='black', fill='#B71C1C', alpha=0.9, lwd=0.2)+
  geom_point(data=HR_rs2_task, aes(x=sequence,y=Mean_HR), alpha=0.5, size=1, fill='black', color='black')+
  geom_errorbar(aes(ymin = m-s, ymax = m+s), width=0.08) +
  ylab("Heart Rate")+xlab("")+
  theme_classic(base_size=16) 
 
# add p-val to plot 
library(ggpubr)
library(rstatix)
HR_rs2_task$sequence <- as.factor(HR_rs2_task$sequence)
stat.test <- HR_rs2_task %>% t_test(Mean_HR ~ sequence, paired = TRUE) %>% add_significance() %>% add_xy_position(x = "sequence") 
pp <- p + stat_pvalue_manual(stat.test, label = "p.signif", bracket.nudge.y = 5, size=5) 
ggsave('Plot_resting2_task.png', pp, width=120, height=100, units=c("mm"), dpi=300)





# paired t-tests

### resting1 vs resting2
HR_resting_only <- filter(HR_all_sequences, sequence!="Task")
t.test(Mean_HR ~ sequence, data = HR_resting_only, paired = TRUE)

### resting2 vs task
HR_resting2_and_task <- filter(HR_all_sequences, sequence!="Resting 1")
t.test(Mean_HR ~ sequence, data = HR_resting2_and_task, paired = TRUE)

### resting1 vs task
HR_resting1_and_task <- filter(HR_all_sequences, sequence!="Resting state")
t.test(Mean_HR ~ sequence, data = HR_resting1_and_task, paired = TRUE)




