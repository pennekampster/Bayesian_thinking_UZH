Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.0")
# 2 for loops:
# 1 for the plots and 1 for the inner plot

# The likelihood will be
# recruits ~ a * pow(adults, b)

# adults are latent variable, need to be continuos
# put a prior to adults as well, it has to be positive and with a mean
library(rjags)
library(BayesianTools)

modelCode = "model{

  for(i in 1:nplots){
    A[i] ~ dunif(0,50)
    
    for(j in 1:subplots){
      rec[4*(i-1)+j] ~ dpois(a*pow(A[i], b))
      obsA[4*(i-1)+j] ~ dpois(A[i])
    }
  }
  
  a ~ dunif(0,100)
  b ~ dunif(0,3)

}"

# Simulate data according to model

nplots = 100
subplots = 4
a = 4
b = .9
density = 1
data = list()
data$nplots = nplots
data$subplots = subplots
data$plot = rep(1:nplots, each = subplots)
data$subplot = rep(1:subplots, nplots)
data$trueA = rep(NA, nplots)
data$obsA = rep(NA, nplots * subplots)
data$rec = rep(NA, nplots * subplots)

for(i in 1:nplots){
  trueA = rlnorm(1, density)
  sel = (4*i-3):(4*i)
  data$trueA[i] = trueA
  data$obsA[sel] = rpois(subplots, trueA)
  data$rec[sel] = rpois(subplots, a * trueA^b)
}

hist(data$obsA)
plot(data$obsA, data$rec)
# run the Jags model with Data2 (NAs for observations -> prior predictive)
jagsModel <- jags.model(file= textConnection(modelCode), data=data, n.chains = 4)
para.names <- c("a","b", "A")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)

x = BayesianTools::getSample(Samples)
marginalPlot(x, which =which(grepl("a|b",colnames(x)))) # sampled prior parameters
par(pty = "s")
plot(apply(x[,1:100],2,median), data$trueA, asp = 1, xlab = "simA", ylab = "trueA")
abline(0,1)


priorPredDistr = x[,1:100] # prior predictive distribution
hist(priorPredDistr, breaks = 100)
hist(log(priorPredDistr), breaks = 100)


###################################################################################################
## The real thing
###################################################################################################

library(data.table)
library(tabulizer)
library(dplyr)
fia_data <- fread("Group_projects/Density_dependence_in_forest_communities/FIAData.txt")[,-10]
names(fia_data) <- c("Latitude", "Longitude", "CellID", "PlotID", "subplot nr.", "species nr.",
                     "adults", "cumBA", "seedlings")

# download species from here
# https://www.fia.fs.fed.us/program-features/urban/docs/Appendix%203a%20CORE%20Species%20List.pdf
species_table <- extract_tables("Group_projects/Density_dependence_in_forest_communities/Appendix 3a CORE Species List.pdf")

species_table2 <- data.table()
for (i_page in species_table) {
  species_table2 <- rbind(species_table2,data.table(i_page[-c(1,2),]))
}
names(species_table2) <-gsub("\\r","",paste0(i_page[c(1),], species_table[[1]][c(2),]))

species_table2[,species := paste(Genus, Species),]
species_table2 <- species_table2[,.("FIAcode" = as.integer(FIACodeCode), species)]

names(fia_data)[names(fia_data) == "species nr."] <- "FIAcode"
names(fia_data)[names(fia_data) == "subplot nr."] <- "subplot"
names(fia_data)[names(fia_data) == "PlotID"] <- "plot"
fia_data <- merge(fia_data,species_table2, by = "FIAcode")

fia_data_sp <- fia_data[species == "Abies balsamea"]

fia_data_sp[,i_plot := as.integer(factor(fia_data_sp$plot)),]
fia_data_sp <- fia_data_sp[order(i_plot)]

str(as.list(fia_data_sp))

data = as.list(fia_data_sp)
data$nsubplots = as.vector(table(data$i_plot))
data$idx = cumsum(data$nsubplots) + 1 -data$nsubplots[1]

data$nplots = uniqueN(data$plot)

data$obsA = fia_data_sp$adults
data$rec = fia_data_sp$seedlings

modelCode = "model{
  
  for(i in nplots){
    A[i] ~ dunif(0,50)
    for(j in nsubplots[i]){
      seedlings[idx[i]+j-1] ~ dpois(a*pow(A[i], b))
      adults[idx[i]+j-1] ~ dpois(A[i])
    }
  }
  a ~ dunif(0,10)
  b ~ dunif(0,2)

}"

# run the Jags model with Data2 (NAs for observations -> prior predictive)
jagsModel <- jags.model(file= textConnection(modelCode), data=data, n.chains = 3)
para.names <- c("a","b", "A")
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)
plot(Samples)

x = BayesianTools::getSample(Samples)
marginalPlot(x, which =which(grepl("a|b",colnames(x)))) # sampled prior parameters
par(pty = "s")
plot(apply(x[,1:100],2,median), data$trueA, asp = 1, xlab = "simA", ylab = "trueA")
abline(0,1)

unique(fia_data_sp$CellID)
par(mfrow = c(3,4))
for(i in 1:uniqueN(fia_data_sp$CellID)){
  {
    plot(seedlings~adults, data = fia_data_sp[CellID == unique(fia_data_sp$CellID)[i]])
    abline(lm(seedlings~adults,data = fia_data_sp[CellID == unique(fia_data_sp$CellID)[i]]))
  }
}

priorPredDistr = x[,1:100] # prior predictive distribution
hist(priorPredDistr, breaks = 100)
hist(log(priorPredDistr), breaks = 100)

