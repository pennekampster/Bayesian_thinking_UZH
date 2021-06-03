library(here)

env <- read.csv(here("Group_projects", "SDM_Rhone_snails","data","gastero_environment.csv"), sep=";", header = T,row.names = 1)
dim(env)
fauna <- read.csv(here("Group_projects", "SDM_Rhone_snails","data","gastero_fauna.csv"), sep=";", header=T, row.names = 1)
dim(fauna)
samples <- read.csv(here("Group_projects", "SDM_Rhone_snails","data","gastero_samples.csv"), sep=";", header =T)
dim(samples)


all_data<-rbind(fauna,env)

all_data<-t(all_data)
all_data<-cbind(all_data,samples)
save(all_data,file=here("Group_projects","SDM_Rhone_snails","data","data.RData"))
