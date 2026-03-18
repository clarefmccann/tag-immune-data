setwd('/Users/clairekelly/Documents/Work/Michelle/Physio_proper_run/')

load("all_pulse_data.RData")
load("all_resp_data.RData")


library(tidyverse)

# plot mean pulse and resp by sequence

library(ggplot2)

all_pulse_data %>% 
  group_by(sequence) %>% 
  summarize(m = mean(mean_pulse),
            s = sd(mean_pulse)) %>% 
  ggplot(aes(x = sequence, y = m)) + 
  geom_point(data=all_pulse_data, aes(x=sequence,y=mean_pulse), alpha=0.5, size=2, fill='black', color='black')+
  geom_bar(stat = "identity", width=1, color='black', fill='#8A3324', alpha=0.5, lwd=0.2)+
  geom_errorbar(aes(ymin = m-s, ymax = m+s), width=0.2) +
  ylab("Mean Pulse")+xlab("")+
  theme_classic(base_size=16) 

all_resp_data %>% 
  group_by(sequence) %>% 
  summarize(m = mean(mean_resp),
            s = sd(mean_resp)) %>% 
  ggplot(aes(x = sequence, y = m)) + 
  geom_point(data=all_resp_data, aes(x=sequence,y=mean_resp), alpha=0.5, size=2, fill='black', color='black')+
  geom_bar(stat = "identity", width=1, color='black', fill='#8A3324', alpha=0.5, lwd=0.2)+
  geom_errorbar(aes(ymin = m-s, ymax = m+s), width=0.2) +
  ylab("Mean Resp")+xlab("")+
  theme_classic(base_size=16) 


# run paired t-tests- to compare pulse and resp of each individual between sequences

### PULSE
### resting1 vs resting2
all_pulse_data_resting_only <- filter(all_pulse_data, sequence!="task")
t.test(mean_pulse ~ sequence, data = all_pulse_data_resting_only, paired = TRUE)
### resting2 vs task
all_pulse_data_resting2_and_task <- filter(all_pulse_data, sequence!="resting1")
t.test(mean_pulse ~ sequence, data = all_pulse_data_resting2_and_task, paired = TRUE)


### RESP
### resting1 vs resting2
all_resp_data_resting_only <- filter(all_resp_data, sequence!="task")
t.test(mean_resp ~ sequence, data = all_resp_data_resting_only, paired = TRUE)
### resting2 vs task
all_resp_data_resting2_and_task <- filter(all_resp_data, sequence!="resting1")
t.test(mean_resp ~ sequence, data = all_resp_data_resting2_and_task, paired = TRUE)
### resting1 vs task
all_resp_data_resting1_and_task <- filter(all_resp_data, sequence!="resting2")
t.test(mean_resp ~ sequence, data = all_resp_data_resting1_and_task, paired = TRUE)




