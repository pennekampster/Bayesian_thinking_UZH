library(here)

env <- read.csv(here("Group_projects", "SDM_Rhone_snails","data","gastero_environment.csv"), sep=";", header = T,row.names = 1)
dim(env)
fauna <- read.csv(here("Group_projects", "SDM_Rhone_snails","data","gastero_fauna.csv"), sep=";", header=T, row.names = 1)
dim(fauna)
samples <- read.csv(here("Group_projects", "SDM_Rhone_snails","data","gastero_samples.csv"), sep=";", header =T)
dim(samples)
restoration <- read.csv(here("Group_projects", "SDM_Rhone_snails", "data", "restoration_info.csv"), sep= ",", header = T)

all_data<-rbind(fauna,env)

all_data<-t(all_data)
all_data<-cbind(all_data,samples)
all_data$year <- all_data$year + 2000
all_data$site <- gsub("AV","DO",all_data$site)
all_data$site <- gsub("AM","UP",all_data$site)

#remove double downstream sites
all_data %>% filter(!grepl("N",all_data$site)) -> all_data

left_join(all_data, restoration, by=c("channel" = "Channel","site" = "Site")) %>% select(Year) -> all_data$restoration_year
left_join(all_data, restoration, by=c("channel" = "Channel","site" = "Site")) %>% select(Type,channel,site) -> a


restored <- all_data %>% filter(year<=all_data$restoration_year) %>% mutate(restoration = "pre")
restored_after <- all_data %>% filter(year>all_data$restoration_year) %>% mutate(restoration = "post")
not_restored <- all_data %>% filter(is.na(restoration_year)) %>% mutate(restoration = NA)

all_data_final <- rbind(restored, restored_after, not_restored)

all_equal(summary(all_data),summary(all_data_final[,-53]))

save(all_data_final,file=here("Group_projects","SDM_Rhone_snails","data","data.RData"))
