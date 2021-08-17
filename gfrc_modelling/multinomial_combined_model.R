# This file fits the binomial any animal model, then a multinomial species model, 
# It then combines output in to give combined probabilities and all output images
# base directory to save all output
base_dir_on_comp <- "C:/Users/christina/"
base_dir_sim <- base_dir_on_comp + "gfrc_models/"
# subdirectory to output and generated images
img_save_dir <- paste0(base_dir_sim, "image_files/")
# this directory should be directory where input data file downloaded from GFRC data repository is saved
file_dir <- paste0(base_dir_sim, "data_files/")
# code currently in multinomial branch of MRSea so need to download package from github to base dir on comp 
# then switch to multinomial branch
devtools::load_all(paste0(base_dir_on_comp, "MRsea"))

library(knitr)
library(tidyverse)
library(kableExtra)
library(viridis)
library(hexbin)
library(car)
library(mvtnorm)

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output <- sub("1\\*", "", output) # Removes 1*10 but leaves anything else raised to 10 
  parse(text=output)
}

rez_for_var_list <- function(var_list_in, init_mod, distMats){
  salsa1dlist_in <- list(
    fitnessMeasure = 'BIC', 
    minKnots_1d=rep(1,length(var_list_in)),
    maxKnots_1d = rep(4,length(var_list_in)),
    startKnots_1d = rep(2,length(var_list_in)),
    degree=rep(2,length(var_list_in)),
    maxIterations = 100, 
    gaps=rep(0,length(var_list_in))
  )
  PA_1d_mod_in <- runSALSA1D(
    initialModel=init_mod, 
    salsa1dlist_in,
    varlist=var_list_in,
    datain=sectionz_ndvi, 
    removal=FALSE,
    suppress.printout=TRUE
  )
  salsa2dlist_in <- list(
    fitnessMeasure = 'BIC', 
    knotgrid = knotgrid, 
    startKnots=6, 
    minKnots=2,
    maxKnots=12, 
    gap=0
  )
  PA_2d_mod_in <- runSALSA2D(
    PA_1d_mod_in$bestModel, 
    salsa2dlist_in, 
    d2k=distMats$dataDist, 
    k2k=distMats$knotDist,
    suppress.printout = TRUE
  )
  an_res_pa_in <- Anova(PA_2d_mod_in$bestModel, test='F')
  rownames(an_res_pa_in) <- c(PA_2d_mod_in$bestModel$varshortnames, "2dsmooth", "Residuals")
  rez_out <- list(mod_2d=PA_2d_mod_in, anova_out=an_res_pa_in)
  return(rez_out)
}

cbbPalette <- c("#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")

#BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
#rmarkdown::render("Data_exploration_update.Rmd")

file_name <- "final_image_data.csv"
sectionz_ndvi<- read.csv(paste0(file_dir, file_name))

### prescence abcense model

sectionz_ndvi$response <- sectionz_ndvi$AnimalsPA
PA_mod <- glm(
  response~1, 
  family='binomial', 
  data=sectionz_ndvi
)

var_list48 <- c("dist_dam_clip", "dist_minor_clip")

set.seed(3004)
knotgrid <- getKnotgrid(coordData = cbind(sectionz_ndvi$x.pos, sectionz_ndvi$y.pos), numKnots = 300, plot=FALSE)

distMats <- makeDists(
  cbind(sectionz_ndvi$x.pos, sectionz_ndvi$y.pos),
  na.omit(knotgrid)
)

salsa1dlist_in <- list(
  fitnessMeasure = 'BIC', 
  minKnots_1d=rep(1,length(var_list48)),
  maxKnots_1d = rep(4,length(var_list48)),
  startKnots_1d = rep(2,length(var_list48)),
  degree=rep(2,length(var_list48)),
  maxIterations = 100, 
  gaps=rep(0,length(var_list48))
)

PA_1d_mod48 <- runSALSA1D(
  initialModel=PA_mod, 
  salsa1dlist_in,
  varlist=var_list48,
  datain=sectionz_ndvi, 
  removal=FALSE,
  suppress.printout=TRUE
)

salsa2dlist_in <- list(
  fitnessMeasure = 'BIC', 
  knotgrid = knotgrid, 
  startKnots=6, 
  minKnots=2,
  maxKnots=12, 
  gap=0
)

PA_2d_mod48 <- runSALSA2D(
  PA_1d_mod48$bestModel, 
  salsa2dlist_in, 
  d2k=distMats$dataDist, 
  k2k=distMats$knotDist,
  suppress.printout = TRUE
)

an_res_pa48 <- Anova(PA_2d_mod48$bestModel, test='F')
rownames(an_res_pa48) <- c(PA_2d_mod48$bestModel$varshortnames, "2dsmooth", "Residuals")

coefz48 <- PA_2d_mod48$bestModel$coefficients
covmat48 <- summary(PA_2d_mod48$bestModel)$cov.robust
rcoefs48 <- rmvnorm(1000, coefz48, sigma=covmat48)

pred_bootz48_all <- predict.gamMRSea(object=PA_2d_mod48$bestModel, coeff=rcoefs48[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=PA_2d_mod48$bestModel, coeff=rcoefs48[bt,])
  pred_bootz48_all <- cbind(pred_bootz48_all, pred_boot)
}
pred_ci48 <- t(apply(pred_bootz48_all, 1, quantile, probs=c(0.025, 0.975)))

pred_48 <- predict.gamMRSea(object=PA_2d_mod48$bestModel)
data_48 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_48, pred_ci48))
colnames(data_48) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


pred_bootz48_all_link <- predict.gamMRSea(object=PA_2d_mod48$bestModel, coeff=rcoefs48[1,], type="link")
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=PA_2d_mod48$bestModel, coeff=rcoefs48[bt,], type="link")
  pred_bootz48_all <- cbind(pred_bootz48_all_link, pred_boot)
}
pred_ci48_link <- t(apply(pred_bootz48_all_link, 1, quantile, probs=c(0.025, 0.975)))

pred_48_link <- predict.gamMRSea(object=PA_2d_mod48$bestModel, type="link")
data_48_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_48_link, pred_ci48_link))
colnames(data_48_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


### Prepare data for multinomial

sectionz_ndvi_mn <- sectionz_ndvi[sectionz_ndvi$AnimalsPA, ]
sectionz_ndvi_mn$SpringbokProb <- sectionz_ndvi_mn$Springbok / sectionz_ndvi_mn$Animals
sectionz_ndvi_mn$OryxProb <- sectionz_ndvi_mn$Oryx / sectionz_ndvi_mn$Animals
sectionz_ndvi_mn$ZebraProb <- sectionz_ndvi_mn$Zebra / sectionz_ndvi_mn$Animals
sectionz_ndvi_mn$KuduProb <- sectionz_ndvi_mn$Kudu / sectionz_ndvi_mn$Animals
sectionz_ndvi_mn$OstrichProb <- sectionz_ndvi_mn$Ostrich / sectionz_ndvi_mn$Animals
sectionz_ndvi_mn$OtherProb <- sectionz_ndvi_mn$Other / sectionz_ndvi_mn$Animals

library(VGAM)

mod_int_only <- vglm(
  cbind(Springbok, Oryx, Zebra, Kudu, Ostrich, Other)~1, 
  family=multinomial, 
  data=sectionz_ndvi_mn
)

knotgrid <- getKnotgrid(coordData = cbind(sectionz_ndvi_mn$x.pos, sectionz_ndvi_mn$y.pos), numKnots = 300)
dists_in <- makeDists(
  cbind(sectionz_ndvi_mn$x.pos, sectionz_ndvi_mn$y.pos),
  na.omit(knotgrid)
)

salsa2dlist <- list(
  fitnessMeasure = "BIC", 
  knotgrid = knotgrid, 
  startKnots=6, 
  minKnots=2,
  maxKnots=15, 
  gap=0,
  noradii=10
)

base_dir_on_comp = "C:/Users/kryzi/"
#base_dir_on_comp = "C:/Users/christina/"

devtools::load_all(paste0(base_dir_on_comp, "Documents/MRsea"))

mod_mn_sim <- runSALSA2Dmn(
  mod_int_only, 
  salsa2dlist,
  vglmdatain=sectionz_ndvi_mn,
  d2k=dists_in$dataDist, 
  k2k=dists_in$knotDist,
  suppress.printout=FALSE,
  chooserad=FALSE,
  hdetest=TRUE
)

sectionz_ndvi_dists <- makeDists(cbind(sectionz_ndvi$x.pos, sectionz_ndvi$y.pos), na.omit(knotgrid))

devtools::load_all(paste0(base_dir_on_comp, "Documents/MRsea"))

mn_function <- function(var_list_in){
  
  salsa1dlist_in <- list(
    fitnessMeasure = 'BIC', 
    minKnots_1d=rep(1,length(var_list_in)),
    maxKnots_1d = rep(4,length(var_list_in)),
    startKnots_1d = rep(2,length(var_list_in)),
    degree=rep(2,length(var_list_in)),
    maxIterations = 100, 
    gaps=rep(0,length(var_list_in))
  )
  
  MN_1d_mod <- runSALSA1Dmn(
    initialModel=mod_int_only, 
    salsa1dlist_in,
    varlist=var_list_in,
    datain=sectionz_ndvi_mn, 
    removal=FALSE,
    suppress.printout=TRUE
  )
  
  salsa2dlist <- list(
    fitnessMeasure = "BIC", 
    knotgrid = knotgrid, 
    startKnots=6, 
    minKnots=2,
    maxKnots=15, 
    gap=0,
    noradii=10
  )
  
  MN_2d_mod <- runSALSA2Dmn(
    MN_1d_mod$bestModel, 
    salsa2dlist, 
    vglmdatain=sectionz_ndvi_mn,
    d2k=dists_in$dataDist, 
    k2k=dists_in$knotDist,
    suppress.printout = FALSE
  )
  
  an_res_mn <- anova(MN_2d_mod$bestModel)
  rez_out <- list(mod_2d=MN_2d_mod, anova_out=an_res_mn)
  
  return(rez_out)
}



### Lowest BIC model

var_list51 <- c("Elev", "dist_river_clip", "NDVI250", "dist_fence_clip", "dist_minor_clip")
mod_2d_51 <- mn_function(var_list51)

pred_mn_bic <- predict.vglmMRSea(newdata=sectionz_ndvi[,c("x.pos","y.pos", var_list51)], object=mod_2d_51$mod_2d$bestModel, newdists=sectionz_ndvi_dists$dataDist, type="response", conf_int=T)

latlong <- sectionz_ndvi[,c("Lat", "Long")]

mn_bic_pred <- cbind(latlong, pred_mn_bic$predictions)
mn_bic_low <- cbind(latlong, pred_mn_bic$lower_limit)
mn_bic_high <- cbind(latlong, pred_mn_bic$higher_limit)

# link scale

pred_mn_bic_link <- predict.vglmMRSea(newdata=sectionz_ndvi[,c("x.pos","y.pos", var_list51)], object=mod_2d_51$mod_2d$bestModel, newdists=sectionz_ndvi_dists$dataDist, type="link", conf_int=T)

mn_bic_pred_link <- cbind(latlong, pred_mn_bic_link$predictions)
mn_bic_low_link <- cbind(latlong, pred_mn_bic_link$lower_limit)
mn_bic_high_link <- cbind(latlong, pred_mn_bic_link$higher_limit)

colnames(mn_bic_pred_link) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")


# 2d smooth only

# filter out any NAs
testcoefs_mn_bic <- mod_2d_51$mod_2d$bestModel@coefficients
cfmsk_mn_bc <- !is.na(testcoefs_mn_bic)
tstcfs_mn_bic <- testcoefs_mn_bic[cfmsk_mn_bc]

sp_col_mn_bc <- mod_2d_51$mod_2d$bestModel@splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_mn_bic <- str_sub(names(tstcfs_mn_bic),1,strlenclnam)
coefmask_mn_bic <- coefnamsplt_mn_bic == nam2d

# create radial gaussian bases
radii_mn_bc <- sp_col_mn_bc$radii
radiiInd_mn_bc <- sp_col_mn_bc$radiusIndices
aR_mn_bc <- sp_col_mn_bc$knotPos
lrf_mn_bc <- LRF.g(radiiInd_mn_bc, distMats$dataDist, radii_mn_bc, aR_mn_bc)

# combine coefmask and facts
coefz_mn_bc <- tstcfs_mn_bic[coefmask_mn_bic]
# get predicted values on link scale
dim(coefz_mn_bc) <- c(5,5)
coefz_mn_bc <- t(coefz_mn_bc)
predtm_mn_bc <- lrf_mn_bc %*% coefz_mn_bc
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

summ_mn_bc <- summary(mod_2d_51$mod_2d$bestModel)
covmat_mn_bc <- summ_mn_bc@cov.unscaled
rcoefs_mn_bc <- rmvnorm(1000, testcoefs_mn_bic, sigma=covmat_mn_bc)


bootcoefz_mn_bc <- rcoefs_mn_bc[1,coefmask_mn_bic]
dim(bootcoefz_mn_bc) <- c(5,5)
bootcoefz_mn_bc <- t(bootcoefz_mn_bc)
predboot_mn_bc_2d <- lrf_mn_bc %*% bootcoefz_mn_bc
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  print(bt)
  bootcf <- rcoefs_mn_bc[bt,coefmask_mn_bic]
  dim(bootcf) <- c(5,5)
  bootcf <- t(bootcf)
  predbt <- lrf_mn_bc %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_mn_bc_2d <- cbind(predboot_mn_bc_2d, predbt)
}

predboot_mn_bc_2d_1 <- predboot_mn_bc_2d[,seq(1,5000,5)]
predboot_mn_bc_2d_2 <- predboot_mn_bc_2d[,seq(2,5000,5)]
predboot_mn_bc_2d_3 <- predboot_mn_bc_2d[,seq(3,5000,5)]
predboot_mn_bc_2d_4 <- predboot_mn_bc_2d[,seq(4,5000,5)]
predboot_mn_bc_2d_5 <- predboot_mn_bc_2d[,seq(5,5000,5)]

pred_2d_ci_mn_bc_1 <- t(apply(predboot_mn_bc_2d_1, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_bc_2 <- t(apply(predboot_mn_bc_2d_2, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_bc_3 <- t(apply(predboot_mn_bc_2d_3, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_bc_4 <- t(apply(predboot_mn_bc_2d_4, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_bc_5 <- t(apply(predboot_mn_bc_2d_5, 1, quantile, probs=c(0.025, 0.975)))

data_2d_mn_bc_1 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_bc[,1], pred_2d_ci_mn_bc_1))
colnames(data_2d_mn_bc_1) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_bc_2 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_bc[,2], pred_2d_ci_mn_bc_2))
colnames(data_2d_mn_bc_2) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_bc_3 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_bc[,3], pred_2d_ci_mn_bc_3))
colnames(data_2d_mn_bc_3) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_bc_4 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_bc[,4], pred_2d_ci_mn_bc_4))
colnames(data_2d_mn_bc_4) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_bc_5 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_bc[,5], pred_2d_ci_mn_bc_5))
colnames(data_2d_mn_bc_5) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

# calc values for 1d variables
data_1d_mn_pred <- mn_bic_pred_link[,3:7] - cbind(data_2d_mn_bc_1$Prediction, data_2d_mn_bc_2$Prediction, 
                                                  data_2d_mn_bc_3$Prediction, data_2d_mn_bc_4$Prediction, 
                                                  data_2d_mn_bc_5$Prediction)
data_1d_mn_pred <- cbind(mn_bic_pred_link[,1:2], data_1d_mn_pred)
colnames(data_1d_mn_pred) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")


data_1d_mn_low <- mn_bic_low_link[,3:7] - cbind(data_2d_mn_bc_1$LowerCI, data_2d_mn_bc_2$LowerCI, 
                                                  data_2d_mn_bc_3$LowerCI, data_2d_mn_bc_4$LowerCI, 
                                                  data_2d_mn_bc_5$LowerCI)
data_1d_mn_low <- cbind(mn_bic_low_link[,1:2], data_1d_mn_low)
colnames(data_1d_mn_low) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")


data_1d_mn_high <- mn_bic_high_link[,3:7] - cbind(data_2d_mn_bc_1$UpperCI, data_2d_mn_bc_2$UpperCI, 
                                                  data_2d_mn_bc_3$UpperCI, data_2d_mn_bc_4$UpperCI, 
                                                  data_2d_mn_bc_5$UpperCI)
data_1d_mn_high <- cbind(mn_bic_high_link[,1:2], data_1d_mn_high)
colnames(data_1d_mn_high) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")


# combine models
combined_mn_bic_pred <- pred_mn_bic$predictions * cbind(data_48$Prediction, data_48$Prediction, data_48$Prediction, data_48$Prediction, data_48$Prediction, data_48$Prediction)
combined_mn_bic_low <- pred_mn_bic$lower_limit * cbind(data_48$LowerCI, data_48$LowerCI, data_48$LowerCI, data_48$LowerCI, data_48$LowerCI, data_48$LowerCI)
combined_mn_bic_high <- pred_mn_bic$higher_limit * cbind(data_48$UpperCI, data_48$UpperCI, data_48$UpperCI, data_48$UpperCI, data_48$UpperCI, data_48$UpperCI)

combined_mn_bic_pred <- cbind(latlong, combined_mn_bic_pred)
combined_mn_bic_low <- cbind(latlong, combined_mn_bic_low)
combined_mn_bic_high <- cbind(latlong, combined_mn_bic_high)

mn_bic_pred <- as.data.frame(cbind(mn_bic_pred, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(mn_bic_pred) <- c(colnames(mn_bic_pred[1:8]), "UTMX", "UTMY", "distX", "distY")
mn_bic_low <- as.data.frame(cbind(mn_bic_low, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(mn_bic_low) <- c(colnames(mn_bic_low[1:8]), "UTMX", "UTMY", "distX", "distY")
mn_bic_high <- as.data.frame(cbind(mn_bic_high, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(mn_bic_high) <- c(colnames(mn_bic_high[1:8]), "UTMX", "UTMY", "distX", "distY")

mn_bic_pred_link <- as.data.frame(cbind(mn_bic_pred_link, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(mn_bic_pred_link) <- c(colnames(mn_bic_pred_link[1:7]), "UTMX", "UTMY", "distX", "distY")
data_1d_mn_pred <- as.data.frame(cbind(data_1d_mn_pred, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_1d_mn_pred) <- c(colnames(data_1d_mn_pred[1:7]), "UTMX", "UTMY", "distX", "distY")
data_2d_mn_bc_1 <- as.data.frame(cbind(data_2d_mn_bc_1, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_2d_mn_bc_1) <- c(colnames(data_2d_mn_bc_1[1:5]), "UTMX", "UTMY", "distX", "distY")
data_2d_mn_bc_2 <- as.data.frame(cbind(data_2d_mn_bc_2, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_2d_mn_bc_2) <- c(colnames(data_2d_mn_bc_2[1:5]), "UTMX", "UTMY", "distX", "distY")
data_2d_mn_bc_3 <- as.data.frame(cbind(data_2d_mn_bc_3, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_2d_mn_bc_3) <- c(colnames(data_2d_mn_bc_3[1:5]), "UTMX", "UTMY", "distX", "distY")
data_2d_mn_bc_4 <- as.data.frame(cbind(data_2d_mn_bc_4, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_2d_mn_bc_4) <- c(colnames(data_2d_mn_bc_4[1:5]), "UTMX", "UTMY", "distX", "distY")
data_2d_mn_bc_5 <- as.data.frame(cbind(data_2d_mn_bc_5, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_2d_mn_bc_5) <- c(colnames(data_2d_mn_bc_5[1:5]), "UTMX", "UTMY", "distX", "distY")

combined_mn_bic_pred <- as.data.frame(cbind(combined_mn_bic_pred, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(combined_mn_bic_pred) <- c(colnames(combined_mn_bic_pred[1:8]), "UTMX", "UTMY", "distX", "distY")
combined_mn_bic_low <- as.data.frame(cbind(combined_mn_bic_low, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(combined_mn_bic_low) <- c(colnames(combined_mn_bic_low[1:8]), "UTMX", "UTMY", "distX", "distY")
combined_mn_bic_high <- as.data.frame(cbind(combined_mn_bic_high, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(combined_mn_bic_high) <- c(colnames(combined_mn_bic_high[1:8]), "UTMX", "UTMY", "distX", "distY")

library(scales)

# create graphs
bok_pred_mn_bic <- ggplot(mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_mn_bic <- ggplot(mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zeb_pred_mn_bic <- ggplot(mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_mn_bic <- ggplot(mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ost_pred_mn_bic <- ggplot(mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oth_pred_mn_bic <- ggplot(mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

# Confidence interval

bok_low_mn_bic <- ggplot(mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

bok_high_mn_bic <- ggplot(mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_low_mn_bic <- ggplot(mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_high_mn_bic <- ggplot(mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zeb_low_mn_bic <- ggplot(mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zeb_high_mn_bic <- ggplot(mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_low_mn_bic <- ggplot(mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_high_mn_bic <- ggplot(mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ost_low_mn_bic <- ggplot(mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ost_high_mn_bic <- ggplot(mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oth_low_mn_bic <- ggplot(mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oth_high_mn_bic <- ggplot(mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

## Whole model just smooth

pred_mn_bic1 <- ggplot(mn_bic_pred_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=First), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_bic2 <- ggplot(mn_bic_pred_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Second), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_bic3 <- ggplot(mn_bic_pred_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Third), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_bic4 <- ggplot(mn_bic_pred_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fourth), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_bic5 <- ggplot(mn_bic_pred_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fifth), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))


pred_mn_2d_bic1 <- ggplot(data_2d_mn_bc_1, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_2d_bic2 <- ggplot(data_2d_mn_bc_2, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_2d_bic3 <- ggplot(data_2d_mn_bc_3, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_2d_bic4 <- ggplot(data_2d_mn_bc_4, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_2d_bic5 <- ggplot(data_2d_mn_bc_5, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

###

pred_mn_1d_bic1 <- ggplot(data_1d_mn_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=First), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_1d_bic2 <- ggplot(data_1d_mn_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Second), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_1d_bic3 <- ggplot(data_1d_mn_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Third), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_1d_bic4 <- ggplot(data_1d_mn_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fourth), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

pred_mn_1d_bic5 <- ggplot(data_1d_mn_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) +  
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fifth), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

## 1d variables

#c("Elev", "dist_river_clip", "NDVI250", "dist_fence_clip", "dist_minor_clip")

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
ndvi_predvals <- seq(min(sectionz_ndvi$NDVI250), max(sectionz_ndvi$NDVI250), length.out=200)
fence_predvals <- seq(min(sectionz_ndvi$dist_fence_clip), max(sectionz_ndvi$dist_fence_clip), length.out=200)
minor_predvals <- seq(min(sectionz_ndvi$dist_minor_clip), max(sectionz_ndvi$dist_minor_clip), length.out=200)
river_predvals <- seq(min(sectionz_ndvi$dist_river_clip), max(sectionz_ndvi$dist_river_clip), length.out=200)


# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
elev_normal <- (sectionz_ndvi$Elev - mean(sectionz_ndvi$Elev)) / sd(sectionz_ndvi$Elev)
ndvi_normal <- (sectionz_ndvi$NDVI250 - mean(sectionz_ndvi$NDVI250)) / sd(sectionz_ndvi$NDVI250)
fence_normal <- (sectionz_ndvi$dist_fence_clip - mean(sectionz_ndvi$dist_fence_clip)) / sd(sectionz_ndvi$dist_fence_clip)
minor_normal <- (sectionz_ndvi$dist_minor_clip - mean(sectionz_ndvi$dist_minor_clip)) / sd(sectionz_ndvi$dist_minor_clip)
river_normal <- (sectionz_ndvi$dist_river_clip - mean(sectionz_ndvi$dist_river_clip)) / sd(sectionz_ndvi$dist_river_clip)


dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_elev <- (elev_normal - median(elev_normal))^2
dist_med_ndvi <- (ndvi_normal - median(ndvi_normal))^2
dist_med_fence <- (fence_normal - median(fence_normal))^2
dist_med_minor <- (minor_normal - median(minor_normal))^2
dist_med_river <- (river_normal - median(river_normal))^2



dist_med_mn_bc <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_ndvi + dist_med_fence + 
  dist_med_minor + dist_med_river
index_med_mn_bc <- which.min(dist_med_mn_bc)
median_vals_mn_bc <- sectionz_ndvi[index_med_mn_bc,]
distMats_pred_mn_bc <- makeDists(
  cbind(rep(median_vals_mn_bc$x.pos, 200), rep(median_vals_mn_bc$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_mn_bc <- as.data.frame(cbind(pred_mn_bic$predictions, sectionz_ndvi$Elev, sectionz_ndvi$NDVI250, 
                                        sectionz_ndvi$dist_fence_clip, sectionz_ndvi$dist_minor_clip, 
                                        sectionz_ndvi$dist_river_clip))
colnames(rug_vals_mn_bc) <- c("Springbok", "Oryx", "Zebra", "Kudu", "Ostrich", "Other", "Elev", "NDVI250", 
                              "dist_fence", "dist_minor", "dist_river")


predict_mn_bc <- as.data.frame(cbind(rep(median_vals_mn_bc$x.pos, 200), 
                                       rep(median_vals_mn_bc$y.pos, 200), 
                                       rep(median_vals_mn_bc$Elev, 200),
                                       rep(median_vals_mn_bc$NDVI250, 200),
                                       rep(median_vals_mn_bc$dist_fence_clip, 200),
                                       rep(median_vals_mn_bc$dist_minor_clip, 200),
                                       rep(median_vals_mn_bc$dist_river_clip, 200)))
colnames(predict_mn_bc) <- c("x.pos", "y.pos", "Elev", "NDVI250", "dist_fence_clip", 
                               "dist_minor_clip", "dist_river_clip")

predict_elev_mn_bc <- predict_mn_bc
predict_elev_mn_bc$Elev <- elev_predvals
pred_bootz_mn_bc_elev <- predict.vglmMRSea(newdata=predict_elev_mn_bc, object=mod_2d_51$mod_2d$bestModel, newdists=distMats_pred_mn_bc$dataDist, conf_int=T)


pred_df_elev_mn_bc_bok <- cbind(predict_elev_mn_bc, pred_bootz_mn_bc_elev$predictions[,1], pred_bootz_mn_bc_elev$lower_limit[,1], pred_bootz_mn_bc_elev$higher_limit[,1])
colnames(pred_df_elev_mn_bc_bok) <- c(colnames(predict_elev_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_bc_oryx <- cbind(predict_elev_mn_bc, pred_bootz_mn_bc_elev$predictions[,2], pred_bootz_mn_bc_elev$lower_limit[,2], pred_bootz_mn_bc_elev$higher_limit[,2])
colnames(pred_df_elev_mn_bc_oryx) <- c(colnames(predict_elev_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_bc_zeb <- cbind(predict_elev_mn_bc, pred_bootz_mn_bc_elev$predictions[,3], pred_bootz_mn_bc_elev$lower_limit[,3], pred_bootz_mn_bc_elev$higher_limit[,3])
colnames(pred_df_elev_mn_bc_zeb) <- c(colnames(predict_elev_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_bc_kudu <- cbind(predict_elev_mn_bc, pred_bootz_mn_bc_elev$predictions[,4], pred_bootz_mn_bc_elev$lower_limit[,4], pred_bootz_mn_bc_elev$higher_limit[,4])
colnames(pred_df_elev_mn_bc_kudu) <- c(colnames(predict_elev_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_bc_ost <- cbind(predict_elev_mn_bc, pred_bootz_mn_bc_elev$predictions[,5], pred_bootz_mn_bc_elev$lower_limit[,5], pred_bootz_mn_bc_elev$higher_limit[,5])
colnames(pred_df_elev_mn_bc_ost) <- c(colnames(predict_elev_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_bc_oth <- cbind(predict_elev_mn_bc, pred_bootz_mn_bc_elev$predictions[,6], pred_bootz_mn_bc_elev$lower_limit[,6], pred_bootz_mn_bc_elev$higher_limit[,6])
colnames(pred_df_elev_mn_bc_oth) <- c(colnames(predict_elev_mn_bc), "predictions", "LowerCI", "UpperCI")

plotout_elev_mn_bc_bok <- ggplot(pred_df_elev_mn_bc_bok, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_bc_oryx <- ggplot(pred_df_elev_mn_bc_oryx, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_bc_zeb <- ggplot(pred_df_elev_mn_bc_zeb, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_bc_kudu <- ggplot(pred_df_elev_mn_bc_kudu, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_bc_ost <- ggplot(pred_df_elev_mn_bc_ost, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_bc_oth <- ggplot(pred_df_elev_mn_bc_oth, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_elev_mn_bc_bok, rep("Springbok", nrow(pred_df_elev_mn_bc_bok)))
oryx_tmp <- cbind(pred_df_elev_mn_bc_oryx, rep("Oryx", nrow(pred_df_elev_mn_bc_oryx)))
zeb_tmp <- cbind(pred_df_elev_mn_bc_zeb, rep("Zebra", nrow(pred_df_elev_mn_bc_zeb)))
kudu_tmp <- cbind(pred_df_elev_mn_bc_kudu, rep("Kudu", nrow(pred_df_elev_mn_bc_kudu)))
ost_tmp <- cbind(pred_df_elev_mn_bc_ost, rep("Ostrich", nrow(pred_df_elev_mn_bc_ost)))
oth_tmp <- cbind(pred_df_elev_mn_bc_oth, rep("Other", nrow(pred_df_elev_mn_bc_oth)))
colnames(bok_tmp)[11] <- "Animal"
colnames(oryx_tmp)[11] <- "Animal"
colnames(zeb_tmp)[11] <- "Animal"
colnames(kudu_tmp)[11] <- "Animal"
colnames(ost_tmp)[11] <- "Animal"
colnames(oth_tmp)[11] <- "Animal"
pred_df_elev_mn_bc_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_elev_mn_bc <- ggplot(pred_df_elev_mn_bc_comb, aes(x=Elev, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_line() +
  scale_colour_manual(values=cbbPalette) + theme(text = element_text(size=20)) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))



predict_ndvi_mn_bc <- predict_mn_bc
predict_ndvi_mn_bc$NDVI250 <- ndvi_predvals
pred_bootz_mn_bc_ndvi <- predict.vglmMRSea(newdata=predict_ndvi_mn_bc, object=mod_2d_51$mod_2d$bestModel, newdists=distMats_pred_mn_bc$dataDist, conf_int=T)

pred_df_ndvi_mn_bc_bok <- cbind(predict_ndvi_mn_bc, pred_bootz_mn_bc_ndvi$predictions[,1], pred_bootz_mn_bc_ndvi$lower_limit[,1], pred_bootz_mn_bc_ndvi$higher_limit[,1])
colnames(pred_df_ndvi_mn_bc_bok) <- c(colnames(predict_ndvi_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_bc_oryx <- cbind(predict_ndvi_mn_bc, pred_bootz_mn_bc_ndvi$predictions[,2], pred_bootz_mn_bc_ndvi$lower_limit[,2], pred_bootz_mn_bc_ndvi$higher_limit[,2])
colnames(pred_df_ndvi_mn_bc_oryx) <- c(colnames(predict_ndvi_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_bc_zeb <- cbind(predict_ndvi_mn_bc, pred_bootz_mn_bc_ndvi$predictions[,3], pred_bootz_mn_bc_ndvi$lower_limit[,3], pred_bootz_mn_bc_ndvi$higher_limit[,3])
colnames(pred_df_ndvi_mn_bc_zeb) <- c(colnames(predict_ndvi_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_bc_kudu <- cbind(predict_ndvi_mn_bc, pred_bootz_mn_bc_ndvi$predictions[,4], pred_bootz_mn_bc_ndvi$lower_limit[,4], pred_bootz_mn_bc_ndvi$higher_limit[,4])
colnames(pred_df_ndvi_mn_bc_kudu) <- c(colnames(predict_ndvi_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_bc_ost <- cbind(predict_ndvi_mn_bc, pred_bootz_mn_bc_ndvi$predictions[,5], pred_bootz_mn_bc_ndvi$lower_limit[,5], pred_bootz_mn_bc_ndvi$higher_limit[,5])
colnames(pred_df_ndvi_mn_bc_ost) <- c(colnames(predict_ndvi_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_bc_oth <- cbind(predict_ndvi_mn_bc, pred_bootz_mn_bc_ndvi$predictions[,6], pred_bootz_mn_bc_ndvi$lower_limit[,6], pred_bootz_mn_bc_ndvi$higher_limit[,6])
colnames(pred_df_ndvi_mn_bc_oth) <- c(colnames(predict_ndvi_mn_bc), "predictions", "LowerCI", "UpperCI")

plotout_ndvi_mn_bc_bok <- ggplot(pred_df_ndvi_mn_bc_bok, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_bc_oryx <- ggplot(pred_df_ndvi_mn_bc_oryx, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_bc_zeb <- ggplot(pred_df_ndvi_mn_bc_zeb, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_bc_kudu <- ggplot(pred_df_ndvi_mn_bc_kudu, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_bc_ost <- ggplot(pred_df_ndvi_mn_bc_ost, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_bc_oth <- ggplot(pred_df_ndvi_mn_bc_oth, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_ndvi_mn_bc_bok, rep("Springbok", nrow(pred_df_ndvi_mn_bc_bok)))
oryx_tmp <- cbind(pred_df_ndvi_mn_bc_oryx, rep("Oryx", nrow(pred_df_ndvi_mn_bc_oryx)))
zeb_tmp <- cbind(pred_df_ndvi_mn_bc_zeb, rep("Zebra", nrow(pred_df_ndvi_mn_bc_zeb)))
kudu_tmp <- cbind(pred_df_ndvi_mn_bc_kudu, rep("Kudu", nrow(pred_df_ndvi_mn_bc_kudu)))
ost_tmp <- cbind(pred_df_ndvi_mn_bc_ost, rep("Ostrich", nrow(pred_df_ndvi_mn_bc_ost)))
oth_tmp <- cbind(pred_df_ndvi_mn_bc_oth, rep("Other", nrow(pred_df_ndvi_mn_bc_oth)))
colnames(bok_tmp)[11] <- "Animal"
colnames(oryx_tmp)[11] <- "Animal"
colnames(zeb_tmp)[11] <- "Animal"
colnames(kudu_tmp)[11] <- "Animal"
colnames(ost_tmp)[11] <- "Animal"
colnames(oth_tmp)[11] <- "Animal"
pred_df_ndvi_mn_bc_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_ndvi_mn_bc <- ggplot(pred_df_ndvi_mn_bc_comb, aes(x=NDVI250, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_line() +
  scale_colour_manual(values=cbbPalette) + theme(text = element_text(size=20)) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_fence_mn_bc <- predict_mn_bc
predict_fence_mn_bc$dist_fence_clip <- fence_predvals
pred_bootz_mn_bc_fence <- predict.vglmMRSea(newdata=predict_fence_mn_bc, object=mod_2d_51$mod_2d$bestModel, newdists=distMats_pred_mn_bc$dataDist, conf_int=T)

pred_df_fence_mn_bc_bok <- cbind(predict_fence_mn_bc, pred_bootz_mn_bc_fence$predictions[,1], pred_bootz_mn_bc_fence$lower_limit[,1], pred_bootz_mn_bc_fence$higher_limit[,1])
colnames(pred_df_fence_mn_bc_bok) <- c(colnames(predict_fence_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_bc_oryx <- cbind(predict_fence_mn_bc, pred_bootz_mn_bc_fence$predictions[,2], pred_bootz_mn_bc_fence$lower_limit[,2], pred_bootz_mn_bc_fence$higher_limit[,2])
colnames(pred_df_fence_mn_bc_oryx) <- c(colnames(predict_fence_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_bc_zeb <- cbind(predict_fence_mn_bc, pred_bootz_mn_bc_fence$predictions[,3], pred_bootz_mn_bc_fence$lower_limit[,3], pred_bootz_mn_bc_fence$higher_limit[,3])
colnames(pred_df_fence_mn_bc_zeb) <- c(colnames(predict_fence_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_bc_kudu <- cbind(predict_fence_mn_bc, pred_bootz_mn_bc_fence$predictions[,4], pred_bootz_mn_bc_fence$lower_limit[,4], pred_bootz_mn_bc_fence$higher_limit[,4])
colnames(pred_df_fence_mn_bc_kudu) <- c(colnames(predict_fence_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_bc_ost <- cbind(predict_fence_mn_bc, pred_bootz_mn_bc_fence$predictions[,5], pred_bootz_mn_bc_fence$lower_limit[,5], pred_bootz_mn_bc_fence$higher_limit[,5])
colnames(pred_df_fence_mn_bc_ost) <- c(colnames(predict_fence_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_bc_oth <- cbind(predict_fence_mn_bc, pred_bootz_mn_bc_fence$predictions[,6], pred_bootz_mn_bc_fence$lower_limit[,6], pred_bootz_mn_bc_fence$higher_limit[,6])
colnames(pred_df_fence_mn_bc_oth) <- c(colnames(predict_fence_mn_bc), "predictions", "LowerCI", "UpperCI")

plotout_fence_mn_bc_bok <- ggplot(pred_df_fence_mn_bc_bok, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_bc_oryx <- ggplot(pred_df_fence_mn_bc_oryx, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_bc_zeb <- ggplot(pred_df_fence_mn_bc_zeb, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_bc_kudu <- ggplot(pred_df_fence_mn_bc_kudu, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_bc_ost <- ggplot(pred_df_fence_mn_bc_ost, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_bc_oth <- ggplot(pred_df_fence_mn_bc_oth, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_fence_mn_bc_bok, rep("Springbok", nrow(pred_df_fence_mn_bc_bok)))
oryx_tmp <- cbind(pred_df_fence_mn_bc_oryx, rep("Oryx", nrow(pred_df_fence_mn_bc_oryx)))
zeb_tmp <- cbind(pred_df_fence_mn_bc_zeb, rep("Zebra", nrow(pred_df_fence_mn_bc_zeb)))
kudu_tmp <- cbind(pred_df_fence_mn_bc_kudu, rep("Kudu", nrow(pred_df_fence_mn_bc_kudu)))
ost_tmp <- cbind(pred_df_fence_mn_bc_ost, rep("Ostrich", nrow(pred_df_fence_mn_bc_ost)))
oth_tmp <- cbind(pred_df_fence_mn_bc_oth, rep("Other", nrow(pred_df_fence_mn_bc_oth)))
colnames(bok_tmp)[11] <- "Animal"
colnames(oryx_tmp)[11] <- "Animal"
colnames(zeb_tmp)[11] <- "Animal"
colnames(kudu_tmp)[11] <- "Animal"
colnames(ost_tmp)[11] <- "Animal"
colnames(oth_tmp)[11] <- "Animal"
pred_df_fence_mn_bc_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_fence_mn_bc <- ggplot(pred_df_fence_mn_bc_comb, aes(x=dist_fence_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to fence') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_minor_mn_bc <- predict_mn_bc
predict_minor_mn_bc$dist_minor_clip <- minor_predvals
pred_bootz_mn_bc_minor <- predict.vglmMRSea(newdata=predict_minor_mn_bc, object=mod_2d_51$mod_2d$bestModel, newdists=distMats_pred_mn_bc$dataDist, conf_int=T)

pred_df_minor_mn_bc_bok <- cbind(predict_minor_mn_bc, pred_bootz_mn_bc_minor$predictions[,1], pred_bootz_mn_bc_minor$lower_limit[,1], pred_bootz_mn_bc_minor$higher_limit[,1])
colnames(pred_df_minor_mn_bc_bok) <- c(colnames(predict_minor_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_bc_oryx <- cbind(predict_minor_mn_bc, pred_bootz_mn_bc_minor$predictions[,2], pred_bootz_mn_bc_minor$lower_limit[,2], pred_bootz_mn_bc_minor$higher_limit[,2])
colnames(pred_df_minor_mn_bc_oryx) <- c(colnames(predict_minor_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_bc_zeb <- cbind(predict_minor_mn_bc, pred_bootz_mn_bc_minor$predictions[,3], pred_bootz_mn_bc_minor$lower_limit[,3], pred_bootz_mn_bc_minor$higher_limit[,3])
colnames(pred_df_minor_mn_bc_zeb) <- c(colnames(predict_minor_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_bc_kudu <- cbind(predict_minor_mn_bc, pred_bootz_mn_bc_minor$predictions[,4], pred_bootz_mn_bc_minor$lower_limit[,4], pred_bootz_mn_bc_minor$higher_limit[,4])
colnames(pred_df_minor_mn_bc_kudu) <- c(colnames(predict_minor_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_bc_ost <- cbind(predict_minor_mn_bc, pred_bootz_mn_bc_minor$predictions[,5], pred_bootz_mn_bc_minor$lower_limit[,5], pred_bootz_mn_bc_minor$higher_limit[,5])
colnames(pred_df_minor_mn_bc_ost) <- c(colnames(predict_minor_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_bc_oth <- cbind(predict_minor_mn_bc, pred_bootz_mn_bc_minor$predictions[,6], pred_bootz_mn_bc_minor$lower_limit[,6], pred_bootz_mn_bc_minor$higher_limit[,6])
colnames(pred_df_minor_mn_bc_oth) <- c(colnames(predict_minor_mn_bc), "predictions", "LowerCI", "UpperCI")

plotout_minor_mn_bc_bok <- ggplot(pred_df_minor_mn_bc_bok, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_bc_oryx <- ggplot(pred_df_minor_mn_bc_oryx, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_bc_zeb <- ggplot(pred_df_minor_mn_bc_zeb, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_bc_kudu <- ggplot(pred_df_minor_mn_bc_kudu, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_bc_ost <- ggplot(pred_df_minor_mn_bc_ost, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_bc_oth <- ggplot(pred_df_minor_mn_bc_oth, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_minor_mn_bc_bok, rep("Springbok", nrow(pred_df_minor_mn_bc_bok)))
oryx_tmp <- cbind(pred_df_minor_mn_bc_oryx, rep("Oryx", nrow(pred_df_minor_mn_bc_oryx)))
zeb_tmp <- cbind(pred_df_minor_mn_bc_zeb, rep("Zebra", nrow(pred_df_minor_mn_bc_zeb)))
kudu_tmp <- cbind(pred_df_minor_mn_bc_kudu, rep("Kudu", nrow(pred_df_minor_mn_bc_kudu)))
ost_tmp <- cbind(pred_df_minor_mn_bc_ost, rep("Ostrich", nrow(pred_df_minor_mn_bc_ost)))
oth_tmp <- cbind(pred_df_minor_mn_bc_oth, rep("Other", nrow(pred_df_minor_mn_bc_oth)))
colnames(bok_tmp)[11] <- "Animal"
colnames(oryx_tmp)[11] <- "Animal"
colnames(zeb_tmp)[11] <- "Animal"
colnames(kudu_tmp)[11] <- "Animal"
colnames(ost_tmp)[11] <- "Animal"
colnames(oth_tmp)[11] <- "Animal"
pred_df_minor_mn_bc_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_minor_mn_bc <- ggplot(pred_df_minor_mn_bc_comb, aes(x=dist_minor_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to minor river') + geom_line() +
  scale_colour_manual(values=cbbPalette) + theme(text = element_text(size=20)) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_river_mn_bc <- predict_mn_bc
predict_river_mn_bc$dist_river_clip <- river_predvals
pred_bootz_mn_bc_river <- predict.vglmMRSea(newdata=predict_river_mn_bc, object=mod_2d_51$mod_2d$bestModel, newdists=distMats_pred_mn_bc$dataDist, conf_int=T)

pred_df_river_mn_bc_bok <- cbind(predict_river_mn_bc, pred_bootz_mn_bc_river$predictions[,1], pred_bootz_mn_bc_river$lower_limit[,1], pred_bootz_mn_bc_river$higher_limit[,1])
colnames(pred_df_river_mn_bc_bok) <- c(colnames(predict_river_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_bc_oryx <- cbind(predict_river_mn_bc, pred_bootz_mn_bc_river$predictions[,2], pred_bootz_mn_bc_river$lower_limit[,2], pred_bootz_mn_bc_river$higher_limit[,2])
colnames(pred_df_river_mn_bc_oryx) <- c(colnames(predict_river_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_bc_zeb <- cbind(predict_river_mn_bc, pred_bootz_mn_bc_river$predictions[,3], pred_bootz_mn_bc_river$lower_limit[,3], pred_bootz_mn_bc_river$higher_limit[,3])
colnames(pred_df_river_mn_bc_zeb) <- c(colnames(predict_river_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_bc_kudu <- cbind(predict_river_mn_bc, pred_bootz_mn_bc_river$predictions[,4], pred_bootz_mn_bc_river$lower_limit[,4], pred_bootz_mn_bc_river$higher_limit[,4])
colnames(pred_df_river_mn_bc_kudu) <- c(colnames(predict_river_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_bc_ost <- cbind(predict_river_mn_bc, pred_bootz_mn_bc_river$predictions[,5], pred_bootz_mn_bc_river$lower_limit[,5], pred_bootz_mn_bc_river$higher_limit[,5])
colnames(pred_df_river_mn_bc_ost) <- c(colnames(predict_river_mn_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_bc_oth <- cbind(predict_river_mn_bc, pred_bootz_mn_bc_river$predictions[,6], pred_bootz_mn_bc_river$lower_limit[,6], pred_bootz_mn_bc_river$higher_limit[,6])
colnames(pred_df_river_mn_bc_oth) <- c(colnames(predict_river_mn_bc), "predictions", "LowerCI", "UpperCI")

plotout_river_mn_bc_bok <- ggplot(pred_df_river_mn_bc_bok, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_bc_oryx <- ggplot(pred_df_river_mn_bc_oryx, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_bc_zeb <- ggplot(pred_df_river_mn_bc_zeb, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_bc_kudu <- ggplot(pred_df_river_mn_bc_kudu, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_bc_ost <- ggplot(pred_df_river_mn_bc_ost, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_bc_oth <- ggplot(pred_df_river_mn_bc_oth, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + theme(text = element_text(size=20)) # + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_river_mn_bc_bok, rep("Springbok", nrow(pred_df_river_mn_bc_bok)))
oryx_tmp <- cbind(pred_df_river_mn_bc_oryx, rep("Oryx", nrow(pred_df_river_mn_bc_oryx)))
zeb_tmp <- cbind(pred_df_river_mn_bc_zeb, rep("Zebra", nrow(pred_df_river_mn_bc_zeb)))
kudu_tmp <- cbind(pred_df_river_mn_bc_kudu, rep("Kudu", nrow(pred_df_river_mn_bc_kudu)))
ost_tmp <- cbind(pred_df_river_mn_bc_ost, rep("Ostrich", nrow(pred_df_river_mn_bc_ost)))
oth_tmp <- cbind(pred_df_river_mn_bc_oth, rep("Other", nrow(pred_df_river_mn_bc_oth)))
colnames(bok_tmp)[11] <- "Animal"
colnames(oryx_tmp)[11] <- "Animal"
colnames(zeb_tmp)[11] <- "Animal"
colnames(kudu_tmp)[11] <- "Animal"
colnames(ost_tmp)[11] <- "Animal"
colnames(oth_tmp)[11] <- "Animal"
pred_df_river_mn_bc_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_river_mn_bc <- ggplot(pred_df_river_mn_bc_comb, aes(x=dist_river_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_line() +
  scale_colour_manual(values=cbbPalette) + theme(text = element_text(size=20)) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))




# Combined models

bok_pred_comb_bic <- ggplot(combined_mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_comb_bic <- ggplot(combined_mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zeb_pred_comb_bic <- ggplot(combined_mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_comb_bic <- ggplot(combined_mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ost_pred_comb_bic <- ggplot(combined_mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oth_pred_comb_bic <- ggplot(combined_mn_bic_pred, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))


bok_low_comb_bic <-  ggplot(combined_mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_low_comb_bic <-  ggplot(combined_mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zeb_low_comb_bic <-  ggplot(combined_mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_low_comb_bic <-  ggplot(combined_mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ost_low_comb_bic <-  ggplot(combined_mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oth_low_comb_bic <-  ggplot(combined_mn_bic_low, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))


bok_high_comb_bic <-  ggplot(combined_mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_high_comb_bic <- ggplot(combined_mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zeb_high_comb_bic <- ggplot(combined_mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_high_comb_bic <- ggplot(combined_mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ost_high_comb_bic <- ggplot(combined_mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oth_high_comb_bic <- ggplot(combined_mn_bic_high, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

# write out graphs
png(filename=paste0(img_save_dir, "MnModelBokPred.png"))
bok_pred_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOryxPred.png"))
oryx_pred_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelZebPred.png"))
zeb_pred_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelKuduPred.png"))
kudu_pred_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOstPred.png"))
ost_pred_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOthPred.png"))
oth_pred_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelBokLow.png"))
bok_low_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOryxLow.png"))
oryx_low_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelZebLow.png"))
zeb_low_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelKuduLow.png"))
kudu_low_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOstLow.png"))
ost_low_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOthLow.png"))
oth_low_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelBokHigh.png"))
bok_high_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOryxHigh.png"))
oryx_high_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelZebHigh.png"))
zeb_high_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelKuduHigh.png"))
kudu_high_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOstHigh.png"))
ost_high_mn_bic
dev.off()
png(filename=paste0(img_save_dir, "MnModelOthHigh.png"))
oth_high_mn_bic
dev.off()

png(filename=paste0(img_save_dir, "MnModelLink1.png"))
pred_mn_bic1
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink2.png"))
pred_mn_bic2
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink3.png"))
pred_mn_bic3
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink4.png"))
pred_mn_bic4
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink5.png"))
pred_mn_bic5
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth1.png"))
pred_mn_2d_bic1
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth2.png"))
pred_mn_2d_bic2
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth3.png"))
pred_mn_2d_bic3
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth4.png"))
pred_mn_2d_bic4
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth5.png"))
pred_mn_2d_bic5
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth1.png"))
pred_mn_1d_bic1
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth2.png"))
pred_mn_1d_bic2
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth3.png"))
pred_mn_1d_bic3
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth4.png"))
pred_mn_1d_bic4
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth5.png"))
pred_mn_1d_bic5
dev.off()

png(filename=paste0(img_save_dir, "MnModelRiver1d.png"))
plotout_river_mn_bc
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1d.png"))
plotout_minor_mn_bc
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1d.png"))
plotout_fence_mn_bc
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1d.png"))
plotout_ndvi_mn_bc
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1d.png"))
plotout_elev_mn_bc
dev.off()

png(filename=paste0(img_save_dir, "MnModelRiver1dBok.png"))
plotout_river_mn_bc_bok
dev.off()
png(filename=paste0(img_save_dir, "MnModelRiver1dOryx.png"))
plotout_river_mn_bc_oryx
dev.off()
png(filename=paste0(img_save_dir, "MnModelRiver1dZeb.png"))
plotout_river_mn_bc_zeb
dev.off()
png(filename=paste0(img_save_dir, "MnModelRiver1dKudu.png"))
plotout_river_mn_bc_kudu
dev.off()
png(filename=paste0(img_save_dir, "MnModelRiver1dOst.png"))
plotout_river_mn_bc_ost
dev.off()
png(filename=paste0(img_save_dir, "MnModelRiver1dOth.png"))
plotout_river_mn_bc_oth
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1dBok.png"))
plotout_minor_mn_bc_bok
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1dOryx.png"))
plotout_minor_mn_bc_oryx
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1dZeb.png"))
plotout_minor_mn_bc_zeb
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1dKudu.png"))
plotout_minor_mn_bc_kudu
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1dOst.png"))
plotout_minor_mn_bc_ost
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1dOth.png"))
plotout_minor_mn_bc_oth
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1dBok.png"))
plotout_fence_mn_bc_bok
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1dOryx.png"))
plotout_fence_mn_bc_oryx
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1dZeb.png"))
plotout_fence_mn_bc_zeb
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1dKudu.png"))
plotout_fence_mn_bc_kudu
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1dOst.png"))
plotout_fence_mn_bc_ost
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1dOth.png"))
plotout_fence_mn_bc_oth
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1dBok.png"))
plotout_ndvi_mn_bc_bok
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1dOryx.png"))
plotout_ndvi_mn_bc_oryx
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1dZeb.png"))
plotout_ndvi_mn_bc_zeb
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1dKudu.png"))
plotout_ndvi_mn_bc_kudu
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1dOst.png"))
plotout_ndvi_mn_bc_ost
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1dOth.png"))
plotout_ndvi_mn_bc_oth
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1dBok.png"))
plotout_elev_mn_bc_bok
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1dOryx.png"))
plotout_elev_mn_bc_oryx
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1dZeb.png"))
plotout_elev_mn_bc_zeb
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1dKudu.png"))
plotout_elev_mn_bc_kudu
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1dOst.png"))
plotout_elev_mn_bc_ost
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1dOth.png"))
plotout_elev_mn_bc_oth
dev.off()



png(filename=paste0(img_save_dir, "HierModelBokPred.png"))
bok_pred_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOryxPred.png"))
oryx_pred_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelZebPred.png"))
zeb_pred_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelKuduPred.png"))
kudu_pred_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOstPred.png"))
ost_pred_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOthPred.png"))
oth_pred_comb_bic
dev.off()

png(filename=paste0(img_save_dir, "HierModelBokLow.png"))
bok_low_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOryxLow.png"))
oryx_low_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelZebLow.png"))
zeb_low_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelKuduLow.png"))
kudu_low_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOstLow.png"))
ost_low_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOthLow.png"))
oth_low_comb_bic
dev.off()

png(filename=paste0(img_save_dir, "HierModelBokHigh.png"))
bok_high_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOryxHigh.png"))
oryx_high_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelZebHigh.png"))
zeb_high_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelKuduHigh.png"))
kudu_high_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOstHigh.png"))
ost_high_comb_bic
dev.off()
png(filename=paste0(img_save_dir, "HierModelOthHigh.png"))
oth_high_comb_bic
dev.off()


#### Not included in thesis graphs not updated


### Model by forward selection

var_list74 <- c("Elev", "dist_river_clip", "NDVI250", "dist_fence_clip", "dist_minor_clip", "dist_track_clip", "dist_wat_clip", "dist_elec_clip", "any_water_clip", "any_road_clip")
mod_2d_74 <- mn_function(var_list74)

pred_mn_fsel <- predict.vglmMRSea(newdata=sectionz_ndvi[,c("x.pos","y.pos", var_list74)], object=mod_2d_74$mod_2d$bestModel, newdists=sectionz_ndvi_dists$dataDist, type="response", conf_int=T)

latlong <- sectionz_ndvi[,c("Lat", "Long")]

mn_fsel_pred <- cbind(latlong, pred_mn_fsel$predictions)
mn_fsel_low <- cbind(latlong, pred_mn_fsel$lower_limit)
mn_fsel_high <- cbind(latlong, pred_mn_fsel$higher_limit)

# link scale

pred_mn_fsel_link <- predict.vglmMRSea(newdata=sectionz_ndvi[,c("x.pos","y.pos", var_list74)], object=mod_2d_74$mod_2d$bestModel, newdists=sectionz_ndvi_dists$dataDist, type="link", conf_int=T)

mn_fsel_pred_link <- cbind(latlong, pred_mn_fsel_link$predictions)
mn_fsel_low_link <- cbind(latlong, pred_mn_fsel_link$lower_limit)
mn_fsel_high_link <- cbind(latlong, pred_mn_fsel_link$higher_limit)

colnames(mn_fsel_pred_link) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")

# 2d smooth only

# filter out any NAs
testcoefs_mn_fsel <- mod_2d_74$mod_2d$bestModel@coefficients
cfmsk_mn_fw <- !is.na(testcoefs_mn_fsel)
tstcfs_mn_fsel <- testcoefs_mn_fsel[cfmsk_mn_fw]

sp_col_mn_fw <- mod_2d_74$mod_2d$bestModel@splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_mn_fsel <- str_sub(names(tstcfs_mn_fsel),1,strlenclnam)
coefmask_mn_fsel <- coefnamsplt_mn_fsel == nam2d

# create radial gaussian bases
radii_mn_fw <- sp_col_mn_fw$radii
radiiInd_mn_fw <- sp_col_mn_fw$radiusIndices
aR_mn_fw <- sp_col_mn_fw$knotPos
lrf_mn_fw <- LRF.g(radiiInd_mn_fw, distMats$dataDist, radii_mn_fw, aR_mn_fw)

# combine coefmask and facts
coefz_mn_fw <- tstcfs_mn_fsel[coefmask_mn_fsel]
# get predicted values on link scale
dim(coefz_mn_fw) <- c(5,6)
coefz_mn_fw <- t(coefz_mn_fw)
predtm_mn_fw <- lrf_mn_fw %*% coefz_mn_fw
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

summ_mn_fw <- summary(mod_2d_74$mod_2d$bestModel)
covmat_mn_fw <- summ_mn_fw@cov.unscaled
rcoefs_mn_fw <- rmvnorm(1000, testcoefs_mn_fsel, sigma=covmat_mn_fw)


bootcoefz_mn_fw <- rcoefs_mn_fw[1,coefmask_mn_fsel]
dim(bootcoefz_mn_fw) <- c(5,6)
bootcoefz_mn_fw <- t(bootcoefz_mn_fw)
predboot_mn_fw_2d <- lrf_mn_fw %*% bootcoefz_mn_fw
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  print(bt)
  bootcf <- rcoefs_mn_fw[bt,coefmask_mn_fsel]
  dim(bootcf) <- c(5,6)
  bootcf <- t(bootcf)
  predbt <- lrf_mn_fw %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_mn_fw_2d <- cbind(predboot_mn_fw_2d, predbt)
}

predboot_mn_fw_2d_1 <- predboot_mn_fw_2d[,seq(1,5000,5)]
predboot_mn_fw_2d_2 <- predboot_mn_fw_2d[,seq(2,5000,5)]
predboot_mn_fw_2d_3 <- predboot_mn_fw_2d[,seq(3,5000,5)]
predboot_mn_fw_2d_4 <- predboot_mn_fw_2d[,seq(4,5000,5)]
predboot_mn_fw_2d_5 <- predboot_mn_fw_2d[,seq(5,5000,5)]

pred_2d_ci_mn_fw_1 <- t(apply(predboot_mn_fw_2d_1, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_fw_2 <- t(apply(predboot_mn_fw_2d_2, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_fw_3 <- t(apply(predboot_mn_fw_2d_3, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_fw_4 <- t(apply(predboot_mn_fw_2d_4, 1, quantile, probs=c(0.025, 0.975)))
pred_2d_ci_mn_fw_5 <- t(apply(predboot_mn_fw_2d_5, 1, quantile, probs=c(0.025, 0.975)))

data_2d_mn_fw_1 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_fw[,1], pred_2d_ci_mn_fw_1))
colnames(data_2d_mn_fw_1) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_fw_2 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_fw[,2], pred_2d_ci_mn_fw_2))
colnames(data_2d_mn_fw_2) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_fw_3 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_fw[,3], pred_2d_ci_mn_fw_3))
colnames(data_2d_mn_fw_3) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_fw_4 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_fw[,4], pred_2d_ci_mn_fw_4))
colnames(data_2d_mn_fw_4) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_2d_mn_fw_5 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_mn_fw[,5], pred_2d_ci_mn_fw_5))
colnames(data_2d_mn_fw_5) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

# calc values for 1d variables
data_1d_mn_pred <- mn_fsel_pred_link[,3:7] - cbind(data_2d_mn_fw_1$Prediction, data_2d_mn_fw_2$Prediction, 
                                                  data_2d_mn_fw_3$Prediction, data_2d_mn_fw_4$Prediction, 
                                                  data_2d_mn_fw_5$Prediction)
data_1d_mn_pred <- cbind(mn_fsel_pred_link[,1:2], data_1d_mn_pred)
colnames(data_1d_mn_pred) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")


data_1d_mn_low <- mn_fsel_low_link[,3:7] - cbind(data_2d_mn_fw_1$LowerCI, data_2d_mn_fw_2$LowerCI, 
                                                data_2d_mn_fw_3$LowerCI, data_2d_mn_fw_4$LowerCI, 
                                                data_2d_mn_fw_5$LowerCI)
data_1d_mn_low <- cbind(mn_fsel_low_link[,1:2], data_1d_mn_low)
colnames(data_1d_mn_low) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")


data_1d_mn_high <- mn_fsel_high_link[,3:7] - cbind(data_2d_mn_fw_1$UpperCI, data_2d_mn_fw_2$UpperCI, 
                                                  data_2d_mn_fw_3$UpperCI, data_2d_mn_fw_4$UpperCI, 
                                                  data_2d_mn_fw_5$UpperCI)
data_1d_mn_high <- cbind(mn_fsel_high_link[,1:2], data_1d_mn_high)
colnames(data_1d_mn_high) <- c("Lat", "Long", "First", "Second", "Third", "Fourth", "Fifth")


# combine models


combined_mn_fsel_pred <- pred_mn_fsel$predictions * cbind(data_48$Prediction, data_48$Prediction, data_48$Prediction, data_48$Prediction, data_48$Prediction, data_48$Prediction)
combined_mn_fsel_low <- pred_mn_fsel$lower_limit * cbind(data_48$LowerCI, data_48$LowerCI, data_48$LowerCI, data_48$LowerCI, data_48$LowerCI, data_48$LowerCI)
combined_mn_fsel_high <- pred_mn_fsel$higher_limit * cbind(data_48$UpperCI, data_48$UpperCI, data_48$UpperCI, data_48$UpperCI, data_48$UpperCI, data_48$UpperCI)

combined_mn_fsel_pred <- cbind(latlong, combined_mn_fsel_pred)
combined_mn_fsel_low <- cbind(latlong, combined_mn_fsel_low)
combined_mn_fsel_high <- cbind(latlong, combined_mn_fsel_high)

bok_pred_mn_fsel <- ggplot(mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Springbok predictions from multinomial model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

oryx_pred_mn_fsel <- ggplot(mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Oryx predictions from multinomial model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

zeb_pred_mn_fsel <- ggplot(mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Zebra predictions from multinomial model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

kudu_pred_mn_fsel <- ggplot(mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Kudu predictions from multinomial model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

ost_pred_mn_fsel <- ggplot(mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Ostrich predictions from multinomial model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

oth_pred_mn_fsel <- ggplot(mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Other predictions from multinomial model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

# Confidence interval

bok_low_mn_fsel <- ggplot(mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Springbok predictions from multinomial model, lower limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

bok_high_mn_fsel <- ggplot(mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Springbok predictions from multinomial model, upper limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

oryx_low_mn_fsel <- ggplot(mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Oryx predictions from multinomial model, lower limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

oryx_high_mn_fsel <- ggplot(mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Oryx predictions from multinomial model, upper limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

zeb_low_mn_fsel <- ggplot(mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Zebra predictions from multinomial model, lower limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

zeb_high_mn_fsel <- ggplot(mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Zebra predictions from multinomial model, upper limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

kudu_low_mn_fsel <- ggplot(mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Kudu predictions from multinomial model, lower limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

kudu_high_mn_fsel <- ggplot(mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Kudu predictions from multinomial model, upper limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

ost_low_mn_fsel <- ggplot(mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Ostrich predictions from multinomial model, lower limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

ost_high_mn_fsel <- ggplot(mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Ostrich predictions from multinomial model, upper limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

oth_low_mn_fsel <- ggplot(mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Other predictions from multinomial model, lower limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

oth_high_mn_fsel <- ggplot(mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Other predictions from multinomial model, upper limit of CI")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Other), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(0,1)) +
  theme(legend.position="right")

## Whole model just smooth

pred_mn_fsel1 <- ggplot(mn_fsel_pred_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for first linear predictor, link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=First), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_fsel2 <- ggplot(mn_fsel_pred_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for second linear predictor, link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Second), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_fsel3 <- ggplot(mn_fsel_pred_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for third linear predictor, link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Third), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_fsel4 <- ggplot(mn_fsel_pred_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for fourth linear predictor, link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fourth), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_fsel5 <- ggplot(mn_fsel_pred_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for fifth linear predictor, link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fifth), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_2d_fsel1 <- ggplot(data_2d_mn_fw_1, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for first linear predictor, 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_2d_fsel2 <- ggplot(data_2d_mn_fw_2, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for second linear predictor, 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_2d_fsel3 <- ggplot(data_2d_mn_fw_3, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for third linear predictor, 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_2d_fsel4 <- ggplot(data_2d_mn_fw_4, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for fourth linear predictor, 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

pred_mn_2d_fsel5 <- ggplot(data_2d_mn_fw_5, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for fifth linear predictor, 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-31,25)) +
  theme(legend.position="right")

###

pred_mn_1d_fsel1 <- ggplot(data_1d_mn_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for first linear predictor, 1d variables")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=First), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-50,100)) +
  theme(legend.position="right")

pred_mn_1d_fsel2 <- ggplot(data_1d_mn_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for second linear predictor, 1d variables")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Second), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-50,100)) +
  theme(legend.position="right")

pred_mn_1d_fsel3 <- ggplot(data_1d_mn_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for third linear predictor, 1d variables")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Third), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-50,100)) +
  theme(legend.position="right")

pred_mn_1d_fsel4 <- ggplot(data_1d_mn_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for fourth linear predictor, 1d variables")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fourth), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-50,100)) +
  theme(legend.position="right")

pred_mn_1d_fsel5 <- ggplot(data_1d_mn_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Predictions from multinomial model for fifth linear predictor, 1d variables")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Fifth), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", limits=c(-50,100)) +
  theme(legend.position="right")

## 1d variables
#c("Elev", "dist_river_clip", "NDVI250", "dist_fence_clip", "dist_minor_clip", 
# "dist_track_clip", "dist_wat_clip", "dist_elec_clip", "any_water_clip", "any_road_clip")

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
ndvi_predvals <- seq(min(sectionz_ndvi$NDVI250), max(sectionz_ndvi$NDVI250), length.out=200)
fence_predvals <- seq(min(sectionz_ndvi$dist_fence_clip), max(sectionz_ndvi$dist_fence_clip), length.out=200)
minor_predvals <- seq(min(sectionz_ndvi$dist_minor_clip), max(sectionz_ndvi$dist_minor_clip), length.out=200)
river_predvals <- seq(min(sectionz_ndvi$dist_river_clip), max(sectionz_ndvi$dist_river_clip), length.out=200)
track_predvals <- seq(min(sectionz_ndvi$dist_track_clip), max(sectionz_ndvi$dist_track_clip), length.out=200)
wat_predvals <- seq(min(sectionz_ndvi$dist_wat_clip), max(sectionz_ndvi$dist_wat_clip), length.out=200)
elec_predvals <- seq(min(sectionz_ndvi$dist_elec_clip), max(sectionz_ndvi$dist_elec_clip), length.out=200)
water_predvals <- seq(min(sectionz_ndvi$any_water_clip), max(sectionz_ndvi$any_water_clip), length.out=200)
anyroad_predvals <- seq(min(sectionz_ndvi$any_road_clip), max(sectionz_ndvi$any_road_clip), length.out=200)


# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
elev_normal <- (sectionz_ndvi$Elev - mean(sectionz_ndvi$Elev)) / sd(sectionz_ndvi$Elev)
ndvi_normal <- (sectionz_ndvi$NDVI250 - mean(sectionz_ndvi$NDVI250)) / sd(sectionz_ndvi$NDVI250)
fence_normal <- (sectionz_ndvi$dist_fence_clip - mean(sectionz_ndvi$dist_fence_clip)) / sd(sectionz_ndvi$dist_fence_clip)
minor_normal <- (sectionz_ndvi$dist_minor_clip - mean(sectionz_ndvi$dist_minor_clip)) / sd(sectionz_ndvi$dist_minor_clip)
river_normal <- (sectionz_ndvi$dist_river_clip - mean(sectionz_ndvi$dist_river_clip)) / sd(sectionz_ndvi$dist_river_clip)
track_normal <- (sectionz_ndvi$dist_track_clip - mean(sectionz_ndvi$dist_track_clip)) / sd(sectionz_ndvi$dist_track_clip)
wat_normal <- (sectionz_ndvi$dist_wat_clip - mean(sectionz_ndvi$dist_wat_clip)) / sd(sectionz_ndvi$dist_wat_clip)
elec_normal <- (sectionz_ndvi$dist_elec_clip - mean(sectionz_ndvi$dist_elec_clip)) / sd(sectionz_ndvi$dist_elec_clip)
water_normal <- (sectionz_ndvi$any_water_clip - mean(sectionz_ndvi$any_water_clip)) / sd(sectionz_ndvi$any_water_clip)
anyroad_normal <- (sectionz_ndvi$any_road_clip - mean(sectionz_ndvi$any_road_clip)) / sd(sectionz_ndvi$any_road_clip)


dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_elev <- (elev_normal - median(elev_normal))^2
dist_med_ndvi <- (ndvi_normal - median(ndvi_normal))^2
dist_med_fence <- (fence_normal - median(fence_normal))^2
dist_med_minor <- (minor_normal - median(minor_normal))^2
dist_med_river <- (river_normal - median(river_normal))^2
dist_med_track <- (track_normal - median(track_normal))^2
dist_med_wat <- (wat_normal - median(wat_normal))^2
dist_med_elec <- (elec_normal - median(elec_normal))^2
dist_med_water <- (water_normal - median(water_normal))^2
dist_med_anyroad <- (anyroad_normal - median(anyroad_normal))^2


dist_med_mn_fw <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_ndvi + dist_med_fence + 
  dist_med_minor + dist_med_river + dist_med_track + dist_med_wat + dist_med_elec + 
  dist_med_water + dist_med_anyroad
index_med_mn_fw <- which.min(dist_med_mn_fw)
median_vals_mn_fw <- sectionz_ndvi[index_med_mn_fw,]
distMats_pred_mn_fw <- makeDists(
  cbind(rep(median_vals_mn_fw$x.pos, 200), rep(median_vals_mn_fw$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_mn_fw <- as.data.frame(cbind(pred_mn_fsel$predictions, sectionz_ndvi$Elev, sectionz_ndvi$NDVI250, 
                                      sectionz_ndvi$dist_fence_clip, sectionz_ndvi$dist_minor_clip, 
                                      sectionz_ndvi$dist_river_clip, sectionz_ndvi$dist_track_clip,
                                      sectionz_ndvi$dist_wat_clip, sectionz_ndvi$dist_elec_clip,
                                      sectionz_ndvi$any_water_clip, sectionz_ndvi$any_road_clip))
colnames(rug_vals_mn_fw) <- c("Springbok", "Oryx", "Zebra", "Kudu", "Ostrich", "Other", "Elev", "NDVI250", 
                              "dist_fence", "dist_minor", "dist_river", "dist_track", "dist_wat",
                              "dist_elec", "dist_water", "dist_anyroad")


predict_mn_fw <- as.data.frame(cbind(rep(median_vals_mn_fw$x.pos, 200), 
                                     rep(median_vals_mn_fw$y.pos, 200), 
                                     rep(median_vals_mn_fw$Elev, 200),
                                     rep(median_vals_mn_fw$NDVI250, 200),
                                     rep(median_vals_mn_fw$dist_fence_clip, 200),
                                     rep(median_vals_mn_fw$dist_minor_clip, 200),
                                     rep(median_vals_mn_fw$dist_river_clip, 200),
                                     rep(median_vals_mn_fw$dist_track_clip, 200),
                                     rep(median_vals_mn_fw$dist_wat_clip, 200),
                                     rep(median_vals_mn_fw$dist_elec_clip, 200),
                                     rep(median_vals_mn_fw$any_water_clip, 200),
                                     rep(median_vals_mn_fw$any_road_clip, 200)))
colnames(predict_mn_fw) <- c("x.pos", "y.pos", "Elev", "NDVI250", "dist_fence_clip", 
                             "dist_minor_clip", "dist_river_clip","dist_track_clip",
                             "dist_wat_clip", "dist_elec_clip", "any_water_clip", "any_road_clip")

predict_elev_mn_fw <- predict_mn_fw
predict_elev_mn_fw$Elev <- elev_predvals
pred_bootz_mn_fw_elev <- predict.vglmMRSea(newdata=predict_elev_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)


pred_df_elev_mn_fw_bok <- cbind(predict_elev_mn_fw, pred_bootz_mn_fw_elev$predictions[,1], pred_bootz_mn_fw_elev$lower_limit[,1], pred_bootz_mn_fw_elev$higher_limit[,1])
colnames(pred_df_elev_mn_fw_bok) <- c(colnames(predict_elev_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_fw_oryx <- cbind(predict_elev_mn_fw, pred_bootz_mn_fw_elev$predictions[,2], pred_bootz_mn_fw_elev$lower_limit[,2], pred_bootz_mn_fw_elev$higher_limit[,2])
colnames(pred_df_elev_mn_fw_oryx) <- c(colnames(predict_elev_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_fw_zeb <- cbind(predict_elev_mn_fw, pred_bootz_mn_fw_elev$predictions[,3], pred_bootz_mn_fw_elev$lower_limit[,3], pred_bootz_mn_fw_elev$higher_limit[,3])
colnames(pred_df_elev_mn_fw_zeb) <- c(colnames(predict_elev_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_fw_kudu <- cbind(predict_elev_mn_fw, pred_bootz_mn_fw_elev$predictions[,4], pred_bootz_mn_fw_elev$lower_limit[,4], pred_bootz_mn_fw_elev$higher_limit[,4])
colnames(pred_df_elev_mn_fw_kudu) <- c(colnames(predict_elev_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_fw_ost <- cbind(predict_elev_mn_fw, pred_bootz_mn_fw_elev$predictions[,5], pred_bootz_mn_fw_elev$lower_limit[,5], pred_bootz_mn_fw_elev$higher_limit[,5])
colnames(pred_df_elev_mn_fw_ost) <- c(colnames(predict_elev_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elev_mn_fw_oth <- cbind(predict_elev_mn_fw, pred_bootz_mn_fw_elev$predictions[,6], pred_bootz_mn_fw_elev$lower_limit[,6], pred_bootz_mn_fw_elev$higher_limit[,6])
colnames(pred_df_elev_mn_fw_oth) <- c(colnames(predict_elev_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_elev_mn_fw_bok <- ggplot(pred_df_elev_mn_fw_bok, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_fw_oryx <- ggplot(pred_df_elev_mn_fw_oryx, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_fw_zeb <- ggplot(pred_df_elev_mn_fw_zeb, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_fw_kudu <- ggplot(pred_df_elev_mn_fw_kudu, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_fw_ost <- ggplot(pred_df_elev_mn_fw_ost, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elev_mn_fw_oth <- ggplot(pred_df_elev_mn_fw_oth, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_elev_mn_fw_bok, rep("Springbok", nrow(pred_df_elev_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_elev_mn_fw_oryx, rep("Oryx", nrow(pred_df_elev_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_elev_mn_fw_zeb, rep("Zebra", nrow(pred_df_elev_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_elev_mn_fw_kudu, rep("Kudu", nrow(pred_df_elev_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_elev_mn_fw_ost, rep("Ostrich", nrow(pred_df_elev_mn_fw_ost)))
oth_tmp <- cbind(pred_df_elev_mn_fw_oth, rep("Other", nrow(pred_df_elev_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_elev_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_elev_mn_fw <- ggplot(pred_df_elev_mn_fw_comb, aes(x=Elev, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_ndvi_mn_fw <- predict_mn_fw
predict_ndvi_mn_fw$NDVI250 <- ndvi_predvals
pred_bootz_mn_fw_ndvi <- predict.vglmMRSea(newdata=predict_ndvi_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_ndvi_mn_fw_bok <- cbind(predict_ndvi_mn_fw, pred_bootz_mn_fw_ndvi$predictions[,1], pred_bootz_mn_fw_ndvi$lower_limit[,1], pred_bootz_mn_fw_ndvi$higher_limit[,1])
colnames(pred_df_ndvi_mn_fw_bok) <- c(colnames(predict_ndvi_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_fw_oryx <- cbind(predict_ndvi_mn_fw, pred_bootz_mn_fw_ndvi$predictions[,2], pred_bootz_mn_fw_ndvi$lower_limit[,2], pred_bootz_mn_fw_ndvi$higher_limit[,2])
colnames(pred_df_ndvi_mn_fw_oryx) <- c(colnames(predict_ndvi_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_fw_zeb <- cbind(predict_ndvi_mn_fw, pred_bootz_mn_fw_ndvi$predictions[,3], pred_bootz_mn_fw_ndvi$lower_limit[,3], pred_bootz_mn_fw_ndvi$higher_limit[,3])
colnames(pred_df_ndvi_mn_fw_zeb) <- c(colnames(predict_ndvi_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_fw_kudu <- cbind(predict_ndvi_mn_fw, pred_bootz_mn_fw_ndvi$predictions[,4], pred_bootz_mn_fw_ndvi$lower_limit[,4], pred_bootz_mn_fw_ndvi$higher_limit[,4])
colnames(pred_df_ndvi_mn_fw_kudu) <- c(colnames(predict_ndvi_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_fw_ost <- cbind(predict_ndvi_mn_fw, pred_bootz_mn_fw_ndvi$predictions[,5], pred_bootz_mn_fw_ndvi$lower_limit[,5], pred_bootz_mn_fw_ndvi$higher_limit[,5])
colnames(pred_df_ndvi_mn_fw_ost) <- c(colnames(predict_ndvi_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_ndvi_mn_fw_oth <- cbind(predict_ndvi_mn_fw, pred_bootz_mn_fw_ndvi$predictions[,6], pred_bootz_mn_fw_ndvi$lower_limit[,6], pred_bootz_mn_fw_ndvi$higher_limit[,6])
colnames(pred_df_ndvi_mn_fw_oth) <- c(colnames(predict_ndvi_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_ndvi_mn_fw_bok <- ggplot(pred_df_ndvi_mn_fw_bok, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_fw_oryx <- ggplot(pred_df_ndvi_mn_fw_oryx, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_fw_zeb <- ggplot(pred_df_ndvi_mn_fw_zeb, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_fw_kudu <- ggplot(pred_df_ndvi_mn_fw_kudu, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_fw_ost <- ggplot(pred_df_ndvi_mn_fw_ost, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_ndvi_mn_fw_oth <- ggplot(pred_df_ndvi_mn_fw_oth, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_ndvi_mn_fw_bok, rep("Springbok", nrow(pred_df_ndvi_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_ndvi_mn_fw_oryx, rep("Oryx", nrow(pred_df_ndvi_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_ndvi_mn_fw_zeb, rep("Zebra", nrow(pred_df_ndvi_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_ndvi_mn_fw_kudu, rep("Kudu", nrow(pred_df_ndvi_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_ndvi_mn_fw_ost, rep("Ostrich", nrow(pred_df_ndvi_mn_fw_ost)))
oth_tmp <- cbind(pred_df_ndvi_mn_fw_oth, rep("Other", nrow(pred_df_ndvi_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_ndvi_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_ndvi_mn_fw <- ggplot(pred_df_ndvi_mn_fw_comb, aes(x=NDVI250, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_fence_mn_fw <- predict_mn_fw
predict_fence_mn_fw$dist_fence_clip <- fence_predvals
pred_bootz_mn_fw_fence <- predict.vglmMRSea(newdata=predict_fence_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_fence_mn_fw_bok <- cbind(predict_fence_mn_fw, pred_bootz_mn_fw_fence$predictions[,1], pred_bootz_mn_fw_fence$lower_limit[,1], pred_bootz_mn_fw_fence$higher_limit[,1])
colnames(pred_df_fence_mn_fw_bok) <- c(colnames(predict_fence_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_fw_oryx <- cbind(predict_fence_mn_fw, pred_bootz_mn_fw_fence$predictions[,2], pred_bootz_mn_fw_fence$lower_limit[,2], pred_bootz_mn_fw_fence$higher_limit[,2])
colnames(pred_df_fence_mn_fw_oryx) <- c(colnames(predict_fence_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_fw_zeb <- cbind(predict_fence_mn_fw, pred_bootz_mn_fw_fence$predictions[,3], pred_bootz_mn_fw_fence$lower_limit[,3], pred_bootz_mn_fw_fence$higher_limit[,3])
colnames(pred_df_fence_mn_fw_zeb) <- c(colnames(predict_fence_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_fw_kudu <- cbind(predict_fence_mn_fw, pred_bootz_mn_fw_fence$predictions[,4], pred_bootz_mn_fw_fence$lower_limit[,4], pred_bootz_mn_fw_fence$higher_limit[,4])
colnames(pred_df_fence_mn_fw_kudu) <- c(colnames(predict_fence_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_fw_ost <- cbind(predict_fence_mn_fw, pred_bootz_mn_fw_fence$predictions[,5], pred_bootz_mn_fw_fence$lower_limit[,5], pred_bootz_mn_fw_fence$higher_limit[,5])
colnames(pred_df_fence_mn_fw_ost) <- c(colnames(predict_fence_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_fence_mn_fw_oth <- cbind(predict_fence_mn_fw, pred_bootz_mn_fw_fence$predictions[,6], pred_bootz_mn_fw_fence$lower_limit[,6], pred_bootz_mn_fw_fence$higher_limit[,6])
colnames(pred_df_fence_mn_fw_oth) <- c(colnames(predict_fence_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_fence_mn_fw_bok <- ggplot(pred_df_fence_mn_fw_bok, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_fw_oryx <- ggplot(pred_df_fence_mn_fw_oryx, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_fw_zeb <- ggplot(pred_df_fence_mn_fw_zeb, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_fw_kudu <- ggplot(pred_df_fence_mn_fw_kudu, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_fw_ost <- ggplot(pred_df_fence_mn_fw_ost, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_fence_mn_fw_oth <- ggplot(pred_df_fence_mn_fw_oth, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to Fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_fence_mn_fw_bok, rep("Springbok", nrow(pred_df_fence_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_fence_mn_fw_oryx, rep("Oryx", nrow(pred_df_fence_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_fence_mn_fw_zeb, rep("Zebra", nrow(pred_df_fence_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_fence_mn_fw_kudu, rep("Kudu", nrow(pred_df_fence_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_fence_mn_fw_ost, rep("Ostrich", nrow(pred_df_fence_mn_fw_ost)))
oth_tmp <- cbind(pred_df_fence_mn_fw_oth, rep("Other", nrow(pred_df_fence_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_fence_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_fence_mn_fw <- ggplot(pred_df_fence_mn_fw_comb, aes(x=dist_fence_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to fence') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_minor_mn_fw <- predict_mn_fw
predict_minor_mn_fw$dist_minor_clip <- minor_predvals
pred_bootz_mn_fw_minor <- predict.vglmMRSea(newdata=predict_minor_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_minor_mn_fw_bok <- cbind(predict_minor_mn_fw, pred_bootz_mn_fw_minor$predictions[,1], pred_bootz_mn_fw_minor$lower_limit[,1], pred_bootz_mn_fw_minor$higher_limit[,1])
colnames(pred_df_minor_mn_fw_bok) <- c(colnames(predict_minor_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_fw_oryx <- cbind(predict_minor_mn_fw, pred_bootz_mn_fw_minor$predictions[,2], pred_bootz_mn_fw_minor$lower_limit[,2], pred_bootz_mn_fw_minor$higher_limit[,2])
colnames(pred_df_minor_mn_fw_oryx) <- c(colnames(predict_minor_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_fw_zeb <- cbind(predict_minor_mn_fw, pred_bootz_mn_fw_minor$predictions[,3], pred_bootz_mn_fw_minor$lower_limit[,3], pred_bootz_mn_fw_minor$higher_limit[,3])
colnames(pred_df_minor_mn_fw_zeb) <- c(colnames(predict_minor_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_fw_kudu <- cbind(predict_minor_mn_fw, pred_bootz_mn_fw_minor$predictions[,4], pred_bootz_mn_fw_minor$lower_limit[,4], pred_bootz_mn_fw_minor$higher_limit[,4])
colnames(pred_df_minor_mn_fw_kudu) <- c(colnames(predict_minor_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_fw_ost <- cbind(predict_minor_mn_fw, pred_bootz_mn_fw_minor$predictions[,5], pred_bootz_mn_fw_minor$lower_limit[,5], pred_bootz_mn_fw_minor$higher_limit[,5])
colnames(pred_df_minor_mn_fw_ost) <- c(colnames(predict_minor_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_minor_mn_fw_oth <- cbind(predict_minor_mn_fw, pred_bootz_mn_fw_minor$predictions[,6], pred_bootz_mn_fw_minor$lower_limit[,6], pred_bootz_mn_fw_minor$higher_limit[,6])
colnames(pred_df_minor_mn_fw_oth) <- c(colnames(predict_minor_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_minor_mn_fw_bok <- ggplot(pred_df_minor_mn_fw_bok, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_fw_oryx <- ggplot(pred_df_minor_mn_fw_oryx, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_fw_zeb <- ggplot(pred_df_minor_mn_fw_zeb, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_fw_kudu <- ggplot(pred_df_minor_mn_fw_kudu, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_fw_ost <- ggplot(pred_df_minor_mn_fw_ost, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_minor_mn_fw_oth <- ggplot(pred_df_minor_mn_fw_oth, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_minor_mn_fw_bok, rep("Springbok", nrow(pred_df_minor_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_minor_mn_fw_oryx, rep("Oryx", nrow(pred_df_minor_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_minor_mn_fw_zeb, rep("Zebra", nrow(pred_df_minor_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_minor_mn_fw_kudu, rep("Kudu", nrow(pred_df_minor_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_minor_mn_fw_ost, rep("Ostrich", nrow(pred_df_minor_mn_fw_ost)))
oth_tmp <- cbind(pred_df_minor_mn_fw_oth, rep("Other", nrow(pred_df_minor_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_minor_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_minor_mn_fw <- ggplot(pred_df_minor_mn_fw_comb, aes(x=dist_minor_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to minor river') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_river_mn_fw <- predict_mn_fw
predict_river_mn_fw$dist_river_clip <- river_predvals
pred_bootz_mn_fw_river <- predict.vglmMRSea(newdata=predict_river_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_river_mn_fw_bok <- cbind(predict_river_mn_fw, pred_bootz_mn_fw_river$predictions[,1], pred_bootz_mn_fw_river$lower_limit[,1], pred_bootz_mn_fw_river$higher_limit[,1])
colnames(pred_df_river_mn_fw_bok) <- c(colnames(predict_river_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_fw_oryx <- cbind(predict_river_mn_fw, pred_bootz_mn_fw_river$predictions[,2], pred_bootz_mn_fw_river$lower_limit[,2], pred_bootz_mn_fw_river$higher_limit[,2])
colnames(pred_df_river_mn_fw_oryx) <- c(colnames(predict_river_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_fw_zeb <- cbind(predict_river_mn_fw, pred_bootz_mn_fw_river$predictions[,3], pred_bootz_mn_fw_river$lower_limit[,3], pred_bootz_mn_fw_river$higher_limit[,3])
colnames(pred_df_river_mn_fw_zeb) <- c(colnames(predict_river_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_fw_kudu <- cbind(predict_river_mn_fw, pred_bootz_mn_fw_river$predictions[,4], pred_bootz_mn_fw_river$lower_limit[,4], pred_bootz_mn_fw_river$higher_limit[,4])
colnames(pred_df_river_mn_fw_kudu) <- c(colnames(predict_river_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_fw_ost <- cbind(predict_river_mn_fw, pred_bootz_mn_fw_river$predictions[,5], pred_bootz_mn_fw_river$lower_limit[,5], pred_bootz_mn_fw_river$higher_limit[,5])
colnames(pred_df_river_mn_fw_ost) <- c(colnames(predict_river_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_river_mn_fw_oth <- cbind(predict_river_mn_fw, pred_bootz_mn_fw_river$predictions[,6], pred_bootz_mn_fw_river$lower_limit[,6], pred_bootz_mn_fw_river$higher_limit[,6])
colnames(pred_df_river_mn_fw_oth) <- c(colnames(predict_river_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_river_mn_fw_bok <- ggplot(pred_df_river_mn_fw_bok, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_fw_oryx <- ggplot(pred_df_river_mn_fw_oryx, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_fw_zeb <- ggplot(pred_df_river_mn_fw_zeb, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_fw_kudu <- ggplot(pred_df_river_mn_fw_kudu, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_fw_ost <- ggplot(pred_df_river_mn_fw_ost, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_river_mn_fw_oth <- ggplot(pred_df_river_mn_fw_oth, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_river_mn_fw_bok, rep("Springbok", nrow(pred_df_river_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_river_mn_fw_oryx, rep("Oryx", nrow(pred_df_river_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_river_mn_fw_zeb, rep("Zebra", nrow(pred_df_river_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_river_mn_fw_kudu, rep("Kudu", nrow(pred_df_river_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_river_mn_fw_ost, rep("Ostrich", nrow(pred_df_river_mn_fw_ost)))
oth_tmp <- cbind(pred_df_river_mn_fw_oth, rep("Other", nrow(pred_df_river_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_river_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_river_mn_fw <- ggplot(pred_df_river_mn_fw_comb, aes(x=dist_river_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to major river') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))



predict_track_mn_fw <- predict_mn_fw
predict_track_mn_fw$dist_track_clip <- track_predvals
pred_bootz_mn_fw_track <- predict.vglmMRSea(newdata=predict_track_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_track_mn_fw_bok <- cbind(predict_track_mn_fw, pred_bootz_mn_fw_track$predictions[,1], pred_bootz_mn_fw_track$lower_limit[,1], pred_bootz_mn_fw_track$higher_limit[,1])
colnames(pred_df_track_mn_fw_bok) <- c(colnames(predict_track_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_track_mn_fw_oryx <- cbind(predict_track_mn_fw, pred_bootz_mn_fw_track$predictions[,2], pred_bootz_mn_fw_track$lower_limit[,2], pred_bootz_mn_fw_track$higher_limit[,2])
colnames(pred_df_track_mn_fw_oryx) <- c(colnames(predict_track_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_track_mn_fw_zeb <- cbind(predict_track_mn_fw, pred_bootz_mn_fw_track$predictions[,3], pred_bootz_mn_fw_track$lower_limit[,3], pred_bootz_mn_fw_track$higher_limit[,3])
colnames(pred_df_track_mn_fw_zeb) <- c(colnames(predict_track_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_track_mn_fw_kudu <- cbind(predict_track_mn_fw, pred_bootz_mn_fw_track$predictions[,4], pred_bootz_mn_fw_track$lower_limit[,4], pred_bootz_mn_fw_track$higher_limit[,4])
colnames(pred_df_track_mn_fw_kudu) <- c(colnames(predict_track_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_track_mn_fw_ost <- cbind(predict_track_mn_fw, pred_bootz_mn_fw_track$predictions[,5], pred_bootz_mn_fw_track$lower_limit[,5], pred_bootz_mn_fw_track$higher_limit[,5])
colnames(pred_df_track_mn_fw_ost) <- c(colnames(predict_track_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_track_mn_fw_oth <- cbind(predict_track_mn_fw, pred_bootz_mn_fw_track$predictions[,6], pred_bootz_mn_fw_track$lower_limit[,6], pred_bootz_mn_fw_track$higher_limit[,6])
colnames(pred_df_track_mn_fw_oth) <- c(colnames(predict_track_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_track_mn_fw_bok <- ggplot(pred_df_track_mn_fw_bok, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_track_mn_fw_oryx <- ggplot(pred_df_track_mn_fw_oryx, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_track_mn_fw_zeb <- ggplot(pred_df_track_mn_fw_zeb, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_track_mn_fw_kudu <- ggplot(pred_df_track_mn_fw_kudu, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_track_mn_fw_ost <- ggplot(pred_df_track_mn_fw_ost, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_track_mn_fw_oth <- ggplot(pred_df_track_mn_fw_oth, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_track_mn_fw_bok, rep("Springbok", nrow(pred_df_track_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_track_mn_fw_oryx, rep("Oryx", nrow(pred_df_track_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_track_mn_fw_zeb, rep("Zebra", nrow(pred_df_track_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_track_mn_fw_kudu, rep("Kudu", nrow(pred_df_track_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_track_mn_fw_ost, rep("Ostrich", nrow(pred_df_track_mn_fw_ost)))
oth_tmp <- cbind(pred_df_track_mn_fw_oth, rep("Other", nrow(pred_df_track_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_track_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_track_mn_fw <- ggplot(pred_df_track_mn_fw_comb, aes(x=dist_track_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_wat_mn_fw <- predict_mn_fw
predict_wat_mn_fw$dist_wat_clip <- wat_predvals
pred_bootz_mn_fw_wat <- predict.vglmMRSea(newdata=predict_wat_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_wat_mn_fw_bok <- cbind(predict_wat_mn_fw, pred_bootz_mn_fw_wat$predictions[,1], pred_bootz_mn_fw_wat$lower_limit[,1], pred_bootz_mn_fw_wat$higher_limit[,1])
colnames(pred_df_wat_mn_fw_bok) <- c(colnames(predict_wat_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_wat_mn_fw_oryx <- cbind(predict_wat_mn_fw, pred_bootz_mn_fw_wat$predictions[,2], pred_bootz_mn_fw_wat$lower_limit[,2], pred_bootz_mn_fw_wat$higher_limit[,2])
colnames(pred_df_wat_mn_fw_oryx) <- c(colnames(predict_wat_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_wat_mn_fw_zeb <- cbind(predict_wat_mn_fw, pred_bootz_mn_fw_wat$predictions[,3], pred_bootz_mn_fw_wat$lower_limit[,3], pred_bootz_mn_fw_wat$higher_limit[,3])
colnames(pred_df_wat_mn_fw_zeb) <- c(colnames(predict_wat_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_wat_mn_fw_kudu <- cbind(predict_wat_mn_fw, pred_bootz_mn_fw_wat$predictions[,4], pred_bootz_mn_fw_wat$lower_limit[,4], pred_bootz_mn_fw_wat$higher_limit[,4])
colnames(pred_df_wat_mn_fw_kudu) <- c(colnames(predict_wat_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_wat_mn_fw_ost <- cbind(predict_wat_mn_fw, pred_bootz_mn_fw_wat$predictions[,5], pred_bootz_mn_fw_wat$lower_limit[,5], pred_bootz_mn_fw_wat$higher_limit[,5])
colnames(pred_df_wat_mn_fw_ost) <- c(colnames(predict_wat_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_wat_mn_fw_oth <- cbind(predict_wat_mn_fw, pred_bootz_mn_fw_wat$predictions[,6], pred_bootz_mn_fw_wat$lower_limit[,6], pred_bootz_mn_fw_wat$higher_limit[,6])
colnames(pred_df_wat_mn_fw_oth) <- c(colnames(predict_wat_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_wat_mn_fw_bok <- ggplot(pred_df_wat_mn_fw_bok, aes(x=dist_wat_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to wat') + geom_ribbon(aes(x=dist_wat_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_wat_mn_fw_oryx <- ggplot(pred_df_wat_mn_fw_oryx, aes(x=dist_wat_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to wat') + geom_ribbon(aes(x=dist_wat_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_wat_mn_fw_zeb <- ggplot(pred_df_wat_mn_fw_zeb, aes(x=dist_wat_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to wat') + geom_ribbon(aes(x=dist_wat_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_wat_mn_fw_kudu <- ggplot(pred_df_wat_mn_fw_kudu, aes(x=dist_wat_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to wat') + geom_ribbon(aes(x=dist_wat_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_wat_mn_fw_ost <- ggplot(pred_df_wat_mn_fw_ost, aes(x=dist_wat_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to wat') + geom_ribbon(aes(x=dist_wat_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_wat_mn_fw_oth <- ggplot(pred_df_wat_mn_fw_oth, aes(x=dist_wat_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to wat') + geom_ribbon(aes(x=dist_wat_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_wat_mn_fw_bok, rep("Springbok", nrow(pred_df_wat_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_wat_mn_fw_oryx, rep("Oryx", nrow(pred_df_wat_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_wat_mn_fw_zeb, rep("Zebra", nrow(pred_df_wat_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_wat_mn_fw_kudu, rep("Kudu", nrow(pred_df_wat_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_wat_mn_fw_ost, rep("Ostrich", nrow(pred_df_wat_mn_fw_ost)))
oth_tmp <- cbind(pred_df_wat_mn_fw_oth, rep("Other", nrow(pred_df_wat_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_wat_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_wat_mn_fw <- ggplot(pred_df_wat_mn_fw_comb, aes(x=dist_wat_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to major wat') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_elec_mn_fw <- predict_mn_fw
predict_elec_mn_fw$dist_elec_clip <- elec_predvals
pred_bootz_mn_fw_elec <- predict.vglmMRSea(newdata=predict_elec_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_elec_mn_fw_bok <- cbind(predict_elec_mn_fw, pred_bootz_mn_fw_elec$predictions[,1], pred_bootz_mn_fw_elec$lower_limit[,1], pred_bootz_mn_fw_elec$higher_limit[,1])
colnames(pred_df_elec_mn_fw_bok) <- c(colnames(predict_elec_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elec_mn_fw_oryx <- cbind(predict_elec_mn_fw, pred_bootz_mn_fw_elec$predictions[,2], pred_bootz_mn_fw_elec$lower_limit[,2], pred_bootz_mn_fw_elec$higher_limit[,2])
colnames(pred_df_elec_mn_fw_oryx) <- c(colnames(predict_elec_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elec_mn_fw_zeb <- cbind(predict_elec_mn_fw, pred_bootz_mn_fw_elec$predictions[,3], pred_bootz_mn_fw_elec$lower_limit[,3], pred_bootz_mn_fw_elec$higher_limit[,3])
colnames(pred_df_elec_mn_fw_zeb) <- c(colnames(predict_elec_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elec_mn_fw_kudu <- cbind(predict_elec_mn_fw, pred_bootz_mn_fw_elec$predictions[,4], pred_bootz_mn_fw_elec$lower_limit[,4], pred_bootz_mn_fw_elec$higher_limit[,4])
colnames(pred_df_elec_mn_fw_kudu) <- c(colnames(predict_elec_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elec_mn_fw_ost <- cbind(predict_elec_mn_fw, pred_bootz_mn_fw_elec$predictions[,5], pred_bootz_mn_fw_elec$lower_limit[,5], pred_bootz_mn_fw_elec$higher_limit[,5])
colnames(pred_df_elec_mn_fw_ost) <- c(colnames(predict_elec_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_elec_mn_fw_oth <- cbind(predict_elec_mn_fw, pred_bootz_mn_fw_elec$predictions[,6], pred_bootz_mn_fw_elec$lower_limit[,6], pred_bootz_mn_fw_elec$higher_limit[,6])
colnames(pred_df_elec_mn_fw_oth) <- c(colnames(predict_elec_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_elec_mn_fw_bok <- ggplot(pred_df_elec_mn_fw_bok, aes(x=dist_elec_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to elec') + geom_ribbon(aes(x=dist_elec_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elec_mn_fw_oryx <- ggplot(pred_df_elec_mn_fw_oryx, aes(x=dist_elec_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to elec') + geom_ribbon(aes(x=dist_elec_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elec_mn_fw_zeb <- ggplot(pred_df_elec_mn_fw_zeb, aes(x=dist_elec_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to elec') + geom_ribbon(aes(x=dist_elec_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elec_mn_fw_kudu <- ggplot(pred_df_elec_mn_fw_kudu, aes(x=dist_elec_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to elec') + geom_ribbon(aes(x=dist_elec_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elec_mn_fw_ost <- ggplot(pred_df_elec_mn_fw_ost, aes(x=dist_elec_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to elec') + geom_ribbon(aes(x=dist_elec_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_elec_mn_fw_oth <- ggplot(pred_df_elec_mn_fw_oth, aes(x=dist_elec_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to elec') + geom_ribbon(aes(x=dist_elec_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_elec_mn_fw_bok, rep("Springbok", nrow(pred_df_elec_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_elec_mn_fw_oryx, rep("Oryx", nrow(pred_df_elec_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_elec_mn_fw_zeb, rep("Zebra", nrow(pred_df_elec_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_elec_mn_fw_kudu, rep("Kudu", nrow(pred_df_elec_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_elec_mn_fw_ost, rep("Ostrich", nrow(pred_df_elec_mn_fw_ost)))
oth_tmp <- cbind(pred_df_elec_mn_fw_oth, rep("Other", nrow(pred_df_elec_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_elec_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_elec_mn_fw <- ggplot(pred_df_elec_mn_fw_comb, aes(x=dist_elec_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to electricity line') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_water_mn_fw <- predict_mn_fw
predict_water_mn_fw$any_water_clip <- water_predvals
pred_bootz_mn_fw_water <- predict.vglmMRSea(newdata=predict_water_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_water_mn_fw_bok <- cbind(predict_water_mn_fw, pred_bootz_mn_fw_water$predictions[,1], pred_bootz_mn_fw_water$lower_limit[,1], pred_bootz_mn_fw_water$higher_limit[,1])
colnames(pred_df_water_mn_fw_bok) <- c(colnames(predict_water_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_water_mn_fw_oryx <- cbind(predict_water_mn_fw, pred_bootz_mn_fw_water$predictions[,2], pred_bootz_mn_fw_water$lower_limit[,2], pred_bootz_mn_fw_water$higher_limit[,2])
colnames(pred_df_water_mn_fw_oryx) <- c(colnames(predict_water_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_water_mn_fw_zeb <- cbind(predict_water_mn_fw, pred_bootz_mn_fw_water$predictions[,3], pred_bootz_mn_fw_water$lower_limit[,3], pred_bootz_mn_fw_water$higher_limit[,3])
colnames(pred_df_water_mn_fw_zeb) <- c(colnames(predict_water_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_water_mn_fw_kudu <- cbind(predict_water_mn_fw, pred_bootz_mn_fw_water$predictions[,4], pred_bootz_mn_fw_water$lower_limit[,4], pred_bootz_mn_fw_water$higher_limit[,4])
colnames(pred_df_water_mn_fw_kudu) <- c(colnames(predict_water_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_water_mn_fw_ost <- cbind(predict_water_mn_fw, pred_bootz_mn_fw_water$predictions[,5], pred_bootz_mn_fw_water$lower_limit[,5], pred_bootz_mn_fw_water$higher_limit[,5])
colnames(pred_df_water_mn_fw_ost) <- c(colnames(predict_water_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_water_mn_fw_oth <- cbind(predict_water_mn_fw, pred_bootz_mn_fw_water$predictions[,6], pred_bootz_mn_fw_water$lower_limit[,6], pred_bootz_mn_fw_water$higher_limit[,6])
colnames(pred_df_water_mn_fw_oth) <- c(colnames(predict_water_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_water_mn_fw_bok <- ggplot(pred_df_water_mn_fw_bok, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to water') + geom_ribbon(aes(x=dist_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_water_mn_fw_oryx <- ggplot(pred_df_water_mn_fw_oryx, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to water') + geom_ribbon(aes(x=dist_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_water_mn_fw_zeb <- ggplot(pred_df_water_mn_fw_zeb, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to water') + geom_ribbon(aes(x=dist_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_water_mn_fw_kudu <- ggplot(pred_df_water_mn_fw_kudu, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to water') + geom_ribbon(aes(x=dist_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_water_mn_fw_ost <- ggplot(pred_df_water_mn_fw_ost, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to water') + geom_ribbon(aes(x=dist_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_water_mn_fw_oth <- ggplot(pred_df_water_mn_fw_oth, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to water') + geom_ribbon(aes(x=dist_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_water_mn_fw_bok, rep("Springbok", nrow(pred_df_water_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_water_mn_fw_oryx, rep("Oryx", nrow(pred_df_water_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_water_mn_fw_zeb, rep("Zebra", nrow(pred_df_water_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_water_mn_fw_kudu, rep("Kudu", nrow(pred_df_water_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_water_mn_fw_ost, rep("Ostrich", nrow(pred_df_water_mn_fw_ost)))
oth_tmp <- cbind(pred_df_water_mn_fw_oth, rep("Other", nrow(pred_df_water_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_water_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_water_mn_fw <- ggplot(pred_df_water_mn_fw_comb, aes(x=any_water_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to any water') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))


predict_anyroad_mn_fw <- predict_mn_fw
predict_anyroad_mn_fw$any_road_clip <- anyroad_predvals
pred_bootz_mn_fw_anyroad <- predict.vglmMRSea(newdata=predict_anyroad_mn_fw, object=mod_2d_74$mod_2d$bestModel, newdists=distMats_pred_mn_fw$dataDist, conf_int=T)

pred_df_anyroad_mn_fw_bok <- cbind(predict_anyroad_mn_fw, pred_bootz_mn_fw_anyroad$predictions[,1], pred_bootz_mn_fw_anyroad$lower_limit[,1], pred_bootz_mn_fw_anyroad$higher_limit[,1])
colnames(pred_df_anyroad_mn_fw_bok) <- c(colnames(predict_anyroad_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_anyroad_mn_fw_oryx <- cbind(predict_anyroad_mn_fw, pred_bootz_mn_fw_anyroad$predictions[,2], pred_bootz_mn_fw_anyroad$lower_limit[,2], pred_bootz_mn_fw_anyroad$higher_limit[,2])
colnames(pred_df_anyroad_mn_fw_oryx) <- c(colnames(predict_anyroad_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_anyroad_mn_fw_zeb <- cbind(predict_anyroad_mn_fw, pred_bootz_mn_fw_anyroad$predictions[,3], pred_bootz_mn_fw_anyroad$lower_limit[,3], pred_bootz_mn_fw_anyroad$higher_limit[,3])
colnames(pred_df_anyroad_mn_fw_zeb) <- c(colnames(predict_anyroad_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_anyroad_mn_fw_kudu <- cbind(predict_anyroad_mn_fw, pred_bootz_mn_fw_anyroad$predictions[,4], pred_bootz_mn_fw_anyroad$lower_limit[,4], pred_bootz_mn_fw_anyroad$higher_limit[,4])
colnames(pred_df_anyroad_mn_fw_kudu) <- c(colnames(predict_anyroad_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_anyroad_mn_fw_ost <- cbind(predict_anyroad_mn_fw, pred_bootz_mn_fw_anyroad$predictions[,5], pred_bootz_mn_fw_anyroad$lower_limit[,5], pred_bootz_mn_fw_anyroad$higher_limit[,5])
colnames(pred_df_anyroad_mn_fw_ost) <- c(colnames(predict_anyroad_mn_fw), "predictions", "LowerCI", "UpperCI")
pred_df_anyroad_mn_fw_oth <- cbind(predict_anyroad_mn_fw, pred_bootz_mn_fw_anyroad$predictions[,6], pred_bootz_mn_fw_anyroad$lower_limit[,6], pred_bootz_mn_fw_anyroad$higher_limit[,6])
colnames(pred_df_anyroad_mn_fw_oth) <- c(colnames(predict_anyroad_mn_fw), "predictions", "LowerCI", "UpperCI")

plotout_anyroad_mn_fw_bok <- ggplot(pred_df_anyroad_mn_fw_bok, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to anyroad') + geom_ribbon(aes(x=dist_anyroad_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_anyroad_mn_fw_oryx <- ggplot(pred_df_anyroad_mn_fw_oryx, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to anyroad') + geom_ribbon(aes(x=dist_anyroad_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_anyroad_mn_fw_zeb <- ggplot(pred_df_anyroad_mn_fw_zeb, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to anyroad') + geom_ribbon(aes(x=dist_anyroad_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_anyroad_mn_fw_kudu <- ggplot(pred_df_anyroad_mn_fw_kudu, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to anyroad') + geom_ribbon(aes(x=dist_anyroad_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_anyroad_mn_fw_ost <- ggplot(pred_df_anyroad_mn_fw_ost, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to anyroad') + geom_ribbon(aes(x=dist_anyroad_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))
plotout_anyroad_mn_fw_oth <- ggplot(pred_df_anyroad_mn_fw_oth, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to anyroad') + geom_ribbon(aes(x=dist_anyroad_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') # + geom_rug(data=rug_vals_mn_fw, sides='b', size=0.5, length=unit(0.02,'npc'))

bok_tmp <- cbind(pred_df_anyroad_mn_fw_bok, rep("Springbok", nrow(pred_df_anyroad_mn_fw_bok)))
oryx_tmp <- cbind(pred_df_anyroad_mn_fw_oryx, rep("Oryx", nrow(pred_df_anyroad_mn_fw_oryx)))
zeb_tmp <- cbind(pred_df_anyroad_mn_fw_zeb, rep("Zebra", nrow(pred_df_anyroad_mn_fw_zeb)))
kudu_tmp <- cbind(pred_df_anyroad_mn_fw_kudu, rep("Kudu", nrow(pred_df_anyroad_mn_fw_kudu)))
ost_tmp <- cbind(pred_df_anyroad_mn_fw_ost, rep("Ostrich", nrow(pred_df_anyroad_mn_fw_ost)))
oth_tmp <- cbind(pred_df_anyroad_mn_fw_oth, rep("Other", nrow(pred_df_anyroad_mn_fw_oth)))
colnames(bok_tmp)[16] <- "Animal"
colnames(oryx_tmp)[16] <- "Animal"
colnames(zeb_tmp)[16] <- "Animal"
colnames(kudu_tmp)[16] <- "Animal"
colnames(ost_tmp)[16] <- "Animal"
colnames(oth_tmp)[16] <- "Animal"
pred_df_anyroad_mn_fw_comb <- rbind(bok_tmp,oryx_tmp,zeb_tmp,kudu_tmp,ost_tmp,oth_tmp)
plotout_anyroad_mn_fw <- ggplot(pred_df_anyroad_mn_fw_comb, aes(x=any_road_clip, y=predictions, color=Animal)) + 
  theme_bw() + ylab('predicted probability') + xlab('Distance to any road') + geom_line() +
  scale_colour_manual(values=cbbPalette) # + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_rug(data=rug_vals_mn_bc, sides='b', size=0.5, length=unit(0.02,'npc'))



# Combined models

bok_pred_comb_fsel <- ggplot(combined_mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Springbok predictions from hierarchical model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

oryx_pred_comb_fsel <- ggplot(combined_mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Oryx predictions from hierarchical model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

zeb_pred_comb_fsel <- ggplot(combined_mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Zebra predictions from hierarchical model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

kudu_pred_comb_fsel <- ggplot(combined_mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Kudu predictions from hierarchical model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

ost_pred_comb_fsel <- ggplot(combined_mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Ostrich predictions from hierarchical model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

oth_pred_comb_fsel <- ggplot(combined_mn_fsel_pred, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Other predictions from hierarchical model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")


bok_low_comb_fsel <- ggplot(combined_mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Springbok predictions from hierarchical model, lower limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

oryx_low_comb_fsel <- ggplot(combined_mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Oryx predictions from hierarchical model, lower limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

zeb_low_comb_fsel <- ggplot(combined_mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Zebra predictions from hierarchical model, lower limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

kudu_low_comb_fsel <- ggplot(combined_mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Kudu predictions from hierarchical model, lower limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

ost_low_comb_fsel <- ggplot(combined_mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Ostrich predictions from hierarchical model, lower limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

oth_low_comb_fsel <- ggplot(combined_mn_fsel_low, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Other predictions from hierarchical model, lower limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")


bok_high_comb_fsel <- ggplot(combined_mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Springbok predictions from hierarchical model, upper limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

oryx_high_comb_fsel <- ggplot(combined_mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Oryx predictions from hierarchical model, upper limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Oryx), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

zeb_high_comb_fsel <- ggplot(combined_mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Zebra predictions from hierarchical model, upper limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Zebra), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

kudu_high_comb_fsel <- ggplot(combined_mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Kudu predictions from hierarchical model, upper limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Kudu), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

ost_high_comb_fsel <- ggplot(combined_mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Ostrich predictions from hierarchical model, upper limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Ostrich), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")

oth_high_comb_fsel <- ggplot(combined_mn_fsel_high, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  ggtitle(paste("Other predictions from hierarchical model, upper limit")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Springbok), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.01,0.02,0.03,0.04,0.05,0.1,1.0), limits=c(0,1)) +
  theme(legend.position="right")


# write out graphs
png(filename=paste0(img_save_dir, "MnModelBokPred2.png"))
bok_pred_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOryxPred2.png"))
oryx_pred_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelZebPred2.png"))
zeb_pred_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelKuduPred2.png"))
kudu_pred_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOstPred2.png"))
ost_pred_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOthPred2.png"))
oth_pred_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelBokLow2.png"))
bok_low_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOryxLow2.png"))
oryx_low_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelZebLow2.png"))
zeb_low_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelKuduLow2.png"))
kudu_low_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOstLow2.png"))
ost_low_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOthLow2.png"))
oth_low_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelBokHigh2.png"))
bok_high_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOryxHigh2.png"))
oryx_high_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelZebHigh2.png"))
zeb_high_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelKuduHigh2.png"))
kudu_high_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOstHigh2.png"))
ost_high_mn_fsel
dev.off()
png(filename=paste0(img_save_dir, "MnModelOthHigh2.png"))
oth_high_mn_fsel
dev.off()

png(filename=paste0(img_save_dir, "MnModelLink1_2.png"))
pred_mn_fsel1
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink2_2.png"))
pred_mn_fsel2
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink3_2.png"))
pred_mn_fsel3
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink4_2.png"))
pred_mn_fsel4
dev.off()
png(filename=paste0(img_save_dir, "MnModelLink5_2.png"))
pred_mn_fsel5
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth1_2.png"))
pred_mn_2d_fsel1
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth2_2.png"))
pred_mn_2d_fsel2
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth3_2.png"))
pred_mn_2d_fsel3
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth4_2.png"))
pred_mn_2d_fsel4
dev.off()
png(filename=paste0(img_save_dir, "MnModel2dSmooth5_2.png"))
pred_mn_2d_fsel5
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth1_2.png"))
pred_mn_1d_fsel1
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth2_2.png"))
pred_mn_1d_fsel2
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth3_2.png"))
pred_mn_1d_fsel3
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth4_2.png"))
pred_mn_1d_fsel4
dev.off()
png(filename=paste0(img_save_dir, "MnModel1dSmooth5_2.png"))
pred_mn_1d_fsel5
dev.off()

png(filename=paste0(img_save_dir, "MnModelRiver1d_2.png"))
plotout_river_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelMinor1d_2.png"))
plotout_minor_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelFence1d_2.png"))
plotout_fence_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelNDVI1d_2.png"))
plotout_ndvi_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelElev1d_2.png"))
plotout_elev_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelTrack1d_2.png"))
plotout_track_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelWat1d_2.png"))
plotout_wat_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelElec1d_2.png"))
plotout_elec_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelWater1d_2.png"))
plotout_water_mn_fw
dev.off()
png(filename=paste0(img_save_dir, "MnModelAnyroad1d_2.png"))
plotout_anyroad_mn_fw
dev.off()

png(filename=paste0(img_save_dir, "HierModelBokPred_2.png"))
bok_pred_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOryxPred_2.png"))
oryx_pred_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelZebPred_2.png"))
zeb_pred_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelKuduPred_2.png"))
kudu_pred_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOstPred_2.png"))
ost_pred_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOthPred_2.png"))
oth_pred_comb_fsel
dev.off()

png(filename=paste0(img_save_dir, "HierModelBokLow_2.png"))
bok_low_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOryxLow_2.png"))
oryx_low_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelZebLow_2.png"))
zeb_low_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelKuduLow_2.png"))
kudu_low_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOstLow_2.png"))
ost_low_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOthLow_2.png"))
oth_low_comb_fsel
dev.off()

png(filename=paste0(img_save_dir, "HierModelBokHigh_2.png"))
bok_high_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOryxHigh_2.png"))
oryx_high_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelZebHigh_2.png"))
zeb_high_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelKuduHigh_2.png"))
kudu_high_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOstHigh_2.png"))
ost_high_comb_fsel
dev.off()
png(filename=paste0(img_save_dir, "HierModelOthHigh_2.png"))
oth_high_comb_fsel
dev.off()
















