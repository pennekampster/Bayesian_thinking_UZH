#Adapted from work by Maeve McGillycuddy

library(glmmTMB)
library(bbmle)

library(mvabund); data(spider)
citation("mvabund")


## organize data into long format

sppTot <- sort(colSums(spider$abund), decreasing = TRUE)

tmp <- cbind(spider$abund, spider$x)

tmp$id <- 1:nrow(tmp)

spiderDat <- reshape(tmp,
           idvar = "id",
           timevar = "Species",
           times =  colnames(spider$abund),
           varying = list(colnames(spider$abund)),
           v.names = "abund",
           direction = "long")

#longdat <- spiderDat
#save(longdat, file="spider_longdat.Rdata")                       
#Alopacce <- subset(spiderDat, Species=="Alopacce")
#save(Alopacce, file="spider_ Alopacce.Rdata")

## fit rank-reduced models with varying dimension

fit_list <- lapply(2:10,
  function(d) {
    fit.rr <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = d),
    data = spiderDat, family=poisson()) #really this should be nbinom2, but it's slower
})

## compare fits via AIC
aic_vec <- sapply(fit_list, AIC)
aic_vec - min(aic_vec, na.rm = TRUE)

#estimated correlations between species across sites
round( attr((VarCorr(fit_list[[4]])$cond$id), "correlation"), 2) 

# Try other models (covariates and distributions)

fit_nb <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 4),
	data = spiderDat, family=nbinom2())

fit_nb_leaves <- glmmTMB(abund ~ Species* fallen.leaves + rr(Species + 0|id, d = 4),
	data = spiderDat, family=nbinom2())

fit_pois <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 4),
	data = spiderDat, family= poisson())

fit_pois_leaves <- glmmTMB(abund ~ Species* fallen.leaves + rr(Species + 0|id, d = 4),
	data = spiderDat, family=poisson())
	
AICtab(fit_nb, fit_nb_leaves, fit_pois, fit_pois_leaves)
	
summary(fit_nb_leaves)

round( attr((VarCorr(fit_nb_leaves)$cond$id), "correlation"), 2) #correlations between species across sites

# Just because that was the best d for one model, 
# doesn't mean it holds for new families and covariates.
# We could fit rank-reduced models with varying dimension
# for our new covariate and family, but that would be
# bleeding-edge research, so we don't need to do that here. 

fit_list <- lapply(2:10,
  function(d) {
  fit.rr <- glmmTMB(abund ~ Species*fallen.leaves + rr(Species + 0|id, d = d),
    data = spiderDat, family=nbinom2()) 
})

## compare fits via AIC
aic_vec <- sapply(fit_list, AIC)
aic_vec - min(aic_vec, na.rm = TRUE)

