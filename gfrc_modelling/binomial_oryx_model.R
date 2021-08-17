# This file fits the binomial oryx model and creates output images
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
library(scales)


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

sectionz_ndvi$response <- sectionz_ndvi$OryxPA
PA_mod <- glm(
  response~1, 
  family='binomial', 
  data=sectionz_ndvi
)
set.seed(3004)
knotgrid <- getKnotgrid(coordData = cbind(sectionz_ndvi$x.pos, sectionz_ndvi$y.pos), numKnots = 300, plot=FALSE)
distMats <- makeDists(
  cbind(sectionz_ndvi$x.pos, sectionz_ndvi$y.pos),
  na.omit(knotgrid)
)

# First model just smooth

salsa2dlist_in <- list(
  fitnessMeasure = 'BIC', 
  knotgrid = knotgrid, 
  startKnots=6, 
  minKnots=2,
  maxKnots=12, 
  gap=0
)
oryx_2d_mod <- runSALSA2D(
  PA_mod, 
  salsa2dlist_in, 
  d2k=distMats$dataDist, 
  k2k=distMats$knotDist,
  suppress.printout = TRUE
)
an_res_oryx <- Anova(oryx_2d_mod$bestModel, test='F')
rownames(an_res_oryx) <- c(oryx_2d_mod$bestModel$varshortnames, "2dsmooth", "Residuals")

coefz_oryx_sm <- oryx_2d_mod$bestModel$coefficients
covmat_oryx_sm <- summary(oryx_2d_mod$bestModel)$cov.robust
rcoefs_oryx_sm <- rmvnorm(1000, coefz_oryx_sm, sigma=covmat_oryx_sm)

pred_bootz_oryx_sm_all <- predict.gamMRSea(object=oryx_2d_mod$bestModel, coeff=rcoefs_oryx_sm[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=oryx_2d_mod$bestModel, coeff=rcoefs_oryx_sm[bt,])
  pred_bootz_oryx_sm_all <- cbind(pred_bootz_oryx_sm_all, pred_boot)
}
pred_ci_oryx_sm <- t(apply(pred_bootz_oryx_sm_all, 1, quantile, probs=c(0.025, 0.975)))

pred_oryx_sm <- predict.gamMRSea(object=oryx_2d_mod$bestModel)
data_oryx_sm <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_oryx_sm, pred_ci_oryx_sm))
colnames(data_oryx_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_oryx_sm <- cbind(data_oryx_sm, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_oryx_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")

oryx_pred_sm <- ggplot(data_oryx_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette =  "YlOrBr", direction=1, limits=c(0, 0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
oryx_pred_sm <- ggplot(data_oryx_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_sm_lw <- ggplot(data_oryx_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0, 0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
oryx_pred_sm_lw <- ggplot(data_oryx_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_sm_hi <- ggplot(data_oryx_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0, 0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
oryx_pred_sm_hi <- ggplot(data_oryx_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OryxModel1Pred.png"))
oryx_pred_sm
dev.off()

png(filename=paste0(img_save_dir, "OryxModel1Lw.png"))
oryx_pred_sm_lw
dev.off()

png(filename=paste0(img_save_dir, "OryxModel1Hi.png"))
oryx_pred_sm_hi
dev.off()

# Second model backward selection (forward selection gives the same)
var_list_back <- c("NDVI250", "Elev", "dist_fence_clip", "dist_track_clip", "any_road_clip", "dist_wat_clip", 
                   "dist_dam_clip", "any_water_clip")
oryx_mod_bc <- rez_for_var_list(var_list_back, PA_mod, distMats)


coefz_oryx_bc <- oryx_mod_bc$mod_2d$bestModel$coefficients
covmat_oryx_bc <- summary(oryx_mod_bc$mod_2d$bestModel)$cov.robust
rcoefs_oryx_bc <- rmvnorm(1000, coefz_oryx_bc, sigma=covmat_oryx_bc)

pred_bootz_oryx_bc_all <- predict.gamMRSea(object=oryx_mod_bc$mod_2d$bestModel, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=oryx_mod_bc$mod_2d$bestModel, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_all <- cbind(pred_bootz_oryx_bc_all, pred_boot)
}
pred_ci_oryx_bc <- t(apply(pred_bootz_oryx_bc_all, 1, quantile, probs=c(0.025, 0.975)))

pred_oryx_bc <- predict.gamMRSea(object=oryx_mod_bc$mod_2d$bestModel)
data_oryx_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_oryx_bc, pred_ci_oryx_bc))
colnames(data_oryx_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_oryx_bc <- cbind(data_oryx_bc, (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_oryx_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY")

graph_pred_oryx_bc <- ggplot(data_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0, 0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_oryx_bc <- ggplot(data_oryx_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlab("distance (km)") + ylab("distance (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_oryx_bc_lw <- ggplot(data_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0, 0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_oryx_bc_lw <- ggplot(data_oryx_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlab("distance (km)") + ylab("distance (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_oryx_bc_hi <- ggplot(data_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0, 0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_oryx_bc_hi <- ggplot(data_oryx_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlab("distance (km)") + ylab("distance (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OryxModel2Pred.png"))
graph_pred_oryx_bc
dev.off()

png(filename=paste0(img_save_dir, "OryxModel2Lw.png"))
graph_pred_oryx_bc_lw
dev.off()

png(filename=paste0(img_save_dir, "OryxModel2Hi.png"))
graph_pred_oryx_bc_hi
dev.off()

### Not included in final thesis not updated from lat long

pred_bootz_oryx_bc_all_link <- predict.gamMRSea(object=oryx_mod_bc$mod_2d$bestModel, coeff=rcoefs_oryx_bc[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=oryx_mod_bc$mod_2d$bestModel, coeff=rcoefs_oryx_bc[bt,], type='link')
  pred_bootz_oryx_bc_all_link <- cbind(pred_bootz_oryx_bc_all, pred_boot)
}
pred_ci_oryx_bc_link <- t(apply(pred_bootz_oryx_bc_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_oryx_bc_link <- predict.gamMRSea(object=oryx_mod_bc$mod_2d$bestModel)
data_oryx_bc_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_oryx_bc_link, pred_ci_oryx_bc_link))
colnames(data_oryx_bc_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")




# filter out any NAs
testcoefs_oryx_bc <- oryx_mod_bc$mod_2d$bestModel$coefficients
cfmsk_oryx_bc <- !is.na(testcoefs_oryx_bc)
tstcfs_oryx_bc <- testcoefs_oryx_bc[cfmsk_oryx_bc]

sp_col_oryx_bc <- oryx_mod_bc$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_oryx_bc <- str_sub(names(tstcfs_oryx_bc),1,strlenclnam)
coefmask_oryx_bc <- coefnamsplt_oryx_bc == nam2d

# create radial gaussian bases
radii_oryx_bc <- sp_col_oryx_bc$radii
radiiInd_oryx_bc <- sp_col_oryx_bc$radiusIndices
aR_oryx_bc <- sp_col_oryx_bc$knotPos
lrf_oryx_bc <- LRF.g(radiiInd_oryx_bc, distMats$dataDist, radii_oryx_bc, aR_oryx_bc)

# combine coefmask and facts
coefz_oryx_bc <- tstcfs_oryx_bc[coefmask_oryx_bc]
# get predicted values on link scale
predtm_oryx_bc <- lrf_oryx_bc %*% coefz_oryx_bc
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_oryx_bc <- rcoefs_oryx_bc[1,][cfmsk_oryx_bc][coefmask_oryx_bc]
predboot_oryx_bc_2d <- lrf_oryx_bc %*% bootcoefz_oryx_bc
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_oryx_bc[bt,][cfmsk_oryx_bc][coefmask_oryx_bc]
  predbt <- lrf_oryx_bc %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_oryx_bc_2d <- cbind(predboot_oryx_bc_2d, predbt)
}
pred_2d_ci_oryx_bc <- t(apply(predboot_oryx_bc_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_oryx_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_oryx_bc, pred_2d_ci_oryx_bc))
colnames(data_2d_oryx_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


oryx_2d_bc <- ggplot(data_2d_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_2d_bc_lw <- ggplot(data_2d_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_2d_bc_hi <- ggplot(data_2d_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette ="YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_oryx_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_oryx_bc_link - predtm_oryx_bc, pred_ci_oryx_bc_link - pred_2d_ci_oryx_bc))
colnames(data_not2d_oryx_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


oryx_not2d_bc <- ggplot(data_not2d_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_not2d_bc_lw <- ggplot(data_not2d_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_not2d_bc_hi <- ggplot(data_not2d_oryx_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))


oryx_pred_bc_link <- ggplot(data_oryx_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_bc_lw_link <- ggplot(data_oryx_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_bc_hi_link <- ggplot(data_oryx_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-16, -4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OryxModel22d.png"))
  oryx_2d_bc
dev.off()
png(filename=paste0(img_save_dir, "OryxModel22dLw.png"))
  oryx_2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "OryxModel22dHi.png"))
  oryx_2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "OryxModel2NotPrd.png"))
  oryx_not2d_bc
dev.off()
png(filename=paste0(img_save_dir, "OryxModel2NotLw.png"))
  oryx_not2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "OryxModel2NotHi.png"))
  oryx_not2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "OryxModel2Lnk.png"))
oryx_pred_bc_link
dev.off()
png(filename=paste0(img_save_dir, "OryxModel2LnkLw.png"))
oryx_pred_bc_lw_link
dev.off()
png(filename=paste0(img_save_dir, "OryxModel2LnkHi.png"))
oryx_pred_bc_hi_link
dev.off()

### Included in thesis below

## 1d calculations

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
ndvi_predvals <- seq(min(sectionz_ndvi$NDVI250), max(sectionz_ndvi$NDVI250), length.out=200)
fence_predvals <- seq(min(sectionz_ndvi$dist_fence_clip), max(sectionz_ndvi$dist_fence_clip), length.out=200)
track_predvals <- seq(min(sectionz_ndvi$dist_track_clip), max(sectionz_ndvi$dist_track_clip), length.out=200)
anyrd_predvals <- seq(min(sectionz_ndvi$any_road_clip), max(sectionz_ndvi$any_road_clip), length.out=200)
wat_predvals <- seq(min(sectionz_ndvi$dist_wat_clip), max(sectionz_ndvi$dist_wat_clip), length.out=200)
dam_predvals <- seq(min(sectionz_ndvi$dist_dam_clip), max(sectionz_ndvi$dist_dam_clip), length.out=200)
anywater_predvals <- seq(min(sectionz_ndvi$any_water_clip), max(sectionz_ndvi$any_water_clip), length.out=200)
#ndvi, elev, fence, track, any road, wat, dam, any water

# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
elev_normal <- (sectionz_ndvi$Elev - mean(sectionz_ndvi$Elev)) / sd(sectionz_ndvi$Elev)
ndvi_normal <- (sectionz_ndvi$NDVI250 - mean(sectionz_ndvi$NDVI250)) / sd(sectionz_ndvi$NDVI250)
fence_normal <- (sectionz_ndvi$dist_fence_clip - mean(sectionz_ndvi$dist_fence_clip)) / sd(sectionz_ndvi$dist_fence_clip)
track_normal <- (sectionz_ndvi$dist_track_clip - mean(sectionz_ndvi$dist_track_clip)) / sd(sectionz_ndvi$dist_track_clip)
anyrd_normal <- (sectionz_ndvi$any_road_clip - mean(sectionz_ndvi$any_road_clip)) / sd(sectionz_ndvi$any_road_clip)
wat_normal <- (sectionz_ndvi$dist_wat_clip - mean(sectionz_ndvi$dist_wat_clip)) / sd(sectionz_ndvi$dist_wat_clip)
dam_normal <- (sectionz_ndvi$dist_dam_clip - mean(sectionz_ndvi$dist_dam_clip)) / sd(sectionz_ndvi$dist_dam_clip)
anywater_normal <- (sectionz_ndvi$any_water_clip - mean(sectionz_ndvi$any_water_clip)) / sd(sectionz_ndvi$any_water_clip)


dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_elev <- (elev_normal - median(elev_normal))^2
dist_med_ndvi <- (ndvi_normal - median(ndvi_normal))^2
dist_med_fence <- (fence_normal - median(fence_normal))^2
dist_med_track <- (track_normal - median(track_normal))^2
dist_med_anyrd <- (anyrd_normal - median(anyrd_normal))^2
dist_med_wat <- (wat_normal - median(wat_normal))^2
dist_med_dam <- (dam_normal - median(dam_normal))^2
dist_med_anywater <- (anywater_normal - median(anywater_normal))^2



dist_med_oryx_bc <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_ndvi + dist_med_fence + 
  dist_med_track + dist_med_anyrd + dist_med_wat + dist_med_dam + dist_med_anywater
index_med_oryx_bc <- which.min(dist_med_oryx_bc)
median_vals_oryx_bc <- sectionz_ndvi[index_med_oryx_bc,]
distMats_pred_oryx_bc <- makeDists(
  cbind(rep(median_vals_oryx_bc$x.pos, 200), rep(median_vals_oryx_bc$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_oryx_bc <- as.data.frame(cbind(pred_oryx_bc, sectionz_ndvi$Elev, sectionz_ndvi$NDVI250, 
                                   sectionz_ndvi$dist_fence_clip, sectionz_ndvi$dist_track_clip, 
                                   sectionz_ndvi$any_road_clip,sectionz_ndvi$dist_wat_clip, 
                                   sectionz_ndvi$dist_dam_clip, sectionz_ndvi$any_water_clip))
colnames(rug_vals_oryx_bc) <- c("predictions", "Elev", "NDVI250", "dist_fence", "dist_track",
                               "any_road", "dist_wat", "dist_dam", "any_water")
rug_vals_oryx_bc$dist_fence_km <- rug_vals_oryx_bc$dist_fence / 1000
rug_vals_oryx_bc$dist_track_km <- rug_vals_oryx_bc$dist_track / 1000
rug_vals_oryx_bc$any_road_km <- rug_vals_oryx_bc$any_road / 1000
rug_vals_oryx_bc$dist_wat_km <- rug_vals_oryx_bc$dist_wat / 1000
rug_vals_oryx_bc$dist_dam_km <- rug_vals_oryx_bc$dist_dam / 1000
rug_vals_oryx_bc$any_water_km <- rug_vals_oryx_bc$any_water / 1000


predict_oryx_bc <- as.data.frame(cbind(rep(median_vals_oryx_bc$x.pos, 200), 
                                          rep(median_vals_oryx_bc$y.pos, 200), 
                                          rep(median_vals_oryx_bc$Elev, 200),
                                          rep(median_vals_oryx_bc$NDVI250, 200),
                                          rep(median_vals_oryx_bc$dist_fence_clip, 200),
                                          rep(median_vals_oryx_bc$dist_track_clip, 200),
                                          rep(median_vals_oryx_bc$any_road_clip, 200),
                                          rep(median_vals_oryx_bc$dist_wat_clip, 200),
                                          rep(median_vals_oryx_bc$dist_dam_clip, 200),
                                          rep(median_vals_oryx_bc$any_water_clip, 200)))
colnames(predict_oryx_bc) <- c("x.pos", "y.pos", "Elev", "NDVI250", "dist_fence_clip", 
                                  "dist_track_clip", "any_road_clip",
                                  "dist_wat_clip", "dist_dam_clip", "any_water_clip")

predict_elev_oryx_bc <- predict_oryx_bc
predict_elev_oryx_bc$Elev <- elev_predvals
pred_elev_oryx_bc <- predict.gamMRSea(newdata=predict_elev_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_elev <- predict.gamMRSea(newdata=predict_elev_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_elev <- cbind(pred_bootz_oryx_bc_elev, pred_boot)
}
pred_ci_elev_oryx_bc <- t(apply(pred_bootz_oryx_bc_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_oryx_bc <- cbind(predict_elev_oryx_bc, pred_elev_oryx_bc, pred_ci_elev_oryx_bc)
colnames(pred_df_elev_oryx_bc) <- c(colnames(predict_elev_oryx_bc), "predictions", "LowerCI", "UpperCI")

predict_ndvi_oryx_bc <- predict_oryx_bc
predict_ndvi_oryx_bc$NDVI250 <- ndvi_predvals
pred_ndvi_oryx_bc <- predict.gamMRSea(newdata=predict_ndvi_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_ndvi <- predict.gamMRSea(newdata=predict_ndvi_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_ndvi_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_ndvi <- cbind(pred_bootz_oryx_bc_ndvi, pred_boot)
}
pred_ci_ndvi_oryx_bc <- t(apply(pred_bootz_oryx_bc_ndvi, 1, quantile, probs=c(0.025, 0.975)))
pred_df_ndvi_oryx_bc <- cbind(predict_ndvi_oryx_bc, pred_ndvi_oryx_bc, pred_ci_ndvi_oryx_bc)
colnames(pred_df_ndvi_oryx_bc) <- c(colnames(predict_ndvi_oryx_bc), "predictions", "LowerCI", "UpperCI")

predict_fence_oryx_bc <- predict_oryx_bc
predict_fence_oryx_bc$dist_fence_clip <- fence_predvals
pred_fence_oryx_bc <- predict.gamMRSea(newdata=predict_fence_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_fence <- predict.gamMRSea(newdata=predict_fence_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_fence_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_fence <- cbind(pred_bootz_oryx_bc_fence, pred_boot)
}
pred_ci_fence_oryx_bc <- t(apply(pred_bootz_oryx_bc_fence, 1, quantile, probs=c(0.025, 0.975)))
pred_df_fence_oryx_bc <- cbind(predict_fence_oryx_bc, pred_fence_oryx_bc, pred_ci_fence_oryx_bc)
colnames(pred_df_fence_oryx_bc) <- c(colnames(predict_fence_oryx_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_oryx_bc$dist_fence_km <- pred_df_fence_oryx_bc$dist_fence_clip / 1000

predict_track_oryx_bc <- predict_oryx_bc
predict_track_oryx_bc$dist_track_clip <- track_predvals
pred_track_oryx_bc <- predict.gamMRSea(newdata=predict_track_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_track <- predict.gamMRSea(newdata=predict_track_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_track <- cbind(pred_bootz_oryx_bc_track, pred_boot)
}
pred_ci_track_oryx_bc <- t(apply(pred_bootz_oryx_bc_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_oryx_bc <- cbind(predict_track_oryx_bc, pred_track_oryx_bc, pred_ci_track_oryx_bc)
colnames(pred_df_track_oryx_bc) <- c(colnames(predict_track_oryx_bc), "predictions", "LowerCI", "UpperCI")
pred_df_track_oryx_bc$dist_track_km <- pred_df_track_oryx_bc$dist_track_clip / 1000

predict_anyrd_oryx_bc <- predict_oryx_bc
predict_anyrd_oryx_bc$any_road_clip <- anyrd_predvals
pred_anyrd_oryx_bc <- predict.gamMRSea(newdata=predict_anyrd_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_anyrd <- predict.gamMRSea(newdata=predict_anyrd_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_anyrd <- cbind(pred_bootz_oryx_bc_anyrd, pred_boot)
}
pred_ci_anyrd_oryx_bc <- t(apply(pred_bootz_oryx_bc_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_oryx_bc <- cbind(predict_anyrd_oryx_bc, pred_anyrd_oryx_bc, pred_ci_anyrd_oryx_bc)
colnames(pred_df_anyrd_oryx_bc) <- c(colnames(predict_anyrd_oryx_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anyrd_oryx_bc$any_road_km <- pred_df_anyrd_oryx_bc$any_road_clip / 1000

predict_wat_oryx_bc <- predict_oryx_bc
predict_wat_oryx_bc$dist_wat_clip <- wat_predvals
pred_wat_oryx_bc <- predict.gamMRSea(newdata=predict_wat_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_wat <- predict.gamMRSea(newdata=predict_wat_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_wat_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_wat <- cbind(pred_bootz_oryx_bc_wat, pred_boot)
}
pred_ci_wat_oryx_bc <- t(apply(pred_bootz_oryx_bc_wat, 1, quantile, probs=c(0.025, 0.975)))
pred_df_wat_oryx_bc <- cbind(predict_wat_oryx_bc, pred_wat_oryx_bc, pred_ci_wat_oryx_bc)
colnames(pred_df_wat_oryx_bc) <- c(colnames(predict_wat_oryx_bc), "predictions", "LowerCI", "UpperCI")
pred_df_wat_oryx_bc$dist_wat_km <- pred_df_wat_oryx_bc$dist_wat_clip / 1000

predict_dam_oryx_bc <- predict_oryx_bc
predict_dam_oryx_bc$dist_dam_clip <- dam_predvals
pred_dam_oryx_bc <- predict.gamMRSea(newdata=predict_dam_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_dam <- predict.gamMRSea(newdata=predict_dam_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_dam_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_dam <- cbind(pred_bootz_oryx_bc_dam, pred_boot)
}
pred_ci_dam_oryx_bc <- t(apply(pred_bootz_oryx_bc_dam, 1, quantile, probs=c(0.025, 0.975)))
pred_df_dam_oryx_bc <- cbind(predict_dam_oryx_bc, pred_dam_oryx_bc, pred_ci_dam_oryx_bc)
colnames(pred_df_dam_oryx_bc) <- c(colnames(predict_dam_oryx_bc), "predictions", "LowerCI", "UpperCI")
pred_df_dam_oryx_bc$dist_dam_km <- pred_df_dam_oryx_bc$dist_dam_clip / 1000

predict_anywt_oryx_bc <- predict_oryx_bc
predict_anywt_oryx_bc$any_water_clip <- anywater_predvals
pred_anywt_oryx_bc <- predict.gamMRSea(newdata=predict_anywt_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist)
pred_bootz_oryx_bc_anywt <- predict.gamMRSea(newdata=predict_anywt_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_oryx_bc, object=oryx_mod_bc$mod_2d$bestModel, g2k=distMats_pred_oryx_bc$dataDist, coeff=rcoefs_oryx_bc[bt,])
  pred_bootz_oryx_bc_anywt <- cbind(pred_bootz_oryx_bc_anywt, pred_boot)
}
pred_ci_anywt_oryx_bc <- t(apply(pred_bootz_oryx_bc_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_oryx_bc <- cbind(predict_anywt_oryx_bc, pred_anywt_oryx_bc, pred_ci_anywt_oryx_bc)
colnames(pred_df_anywt_oryx_bc) <- c(colnames(predict_anywt_oryx_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anywt_oryx_bc$any_water_km <- pred_df_anywt_oryx_bc$any_water_clip / 1000


pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_oryx_bc] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")
df_pointz <- cbind(df_pointz, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(df_pointz) <- c("Lat", "Long", "Points", "UTMX", "UTMY", "distX", "distY")

pt_med <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of median point used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.05, 0.05)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))
pt_med <- ggplot(df_pointz, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(6, 6)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))


plotout_elev_oryx_bc <- ggplot(pred_df_elev_oryx_bc, aes(x=Elev, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Elevation (m)') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc')) +  theme(text = element_text(size=20))
plotout_ndvi_oryx_bc <- ggplot(pred_df_ndvi_oryx_bc, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc')) +  theme(text = element_text(size=20))
plotout_fence_oryx_bc <- ggplot(pred_df_fence_oryx_bc, aes(x=dist_fence_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to fence (km)') + geom_ribbon(aes(x=dist_fence_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc')) +  theme(text = element_text(size=20))
plotout_track_oryx_bc <- ggplot(pred_df_track_oryx_bc, aes(x=dist_track_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to track (km)') + geom_ribbon(aes(x=dist_track_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc')) +  theme(text = element_text(size=20))
plotout_anyrd_oryx_bc <- ggplot(pred_df_anyrd_oryx_bc, aes(x=any_road_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any road (km)') + geom_ribbon(aes(x=any_road_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc')) +  theme(text = element_text(size=20))
plotout_wat_oryx_bc <- ggplot(pred_df_wat_oryx_bc, aes(x=dist_wat_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to waterhole (km)') + geom_ribbon(aes(x=dist_wat_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc'))+  theme(text = element_text(size=20))
plotout_dam_oryx_bc <- ggplot(pred_df_dam_oryx_bc, aes(x=dist_dam_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to dam (km)') + geom_ribbon(aes(x=dist_dam_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc'))+  theme(text = element_text(size=20))
plotout_anywt_oryx_bc <- ggplot(pred_df_anywt_oryx_bc, aes(x=any_water_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any water (km)') + geom_ribbon(aes(x=any_water_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bc, sides='b', size=0.5, length=unit(0.02,'npc')) +  theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed2_oryx.png"))
pt_med
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Elev_oryx.png"))
plotout_elev_oryx_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2NDVI_oryx.png"))
plotout_ndvi_oryx_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Fence_oryx.png"))
plotout_fence_oryx_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Track_oryx.png"))
plotout_track_oryx_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Anyrd_oryx.png"))
plotout_anyrd_oryx_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2wat_oryx.png"))
plotout_wat_oryx_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2dam_oryx.png"))
plotout_dam_oryx_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2anywater_oryx.png"))
plotout_anywt_oryx_bc
dev.off()






## best BIC model
var_list_bic <- c("Elev", "any_water_clip")
oryx_mod_bic <- rez_for_var_list(var_list_bic, PA_mod, distMats)


coefz_oryx_bic <- oryx_mod_bic$mod_2d$bestModel$coefficients
covmat_oryx_bic <- summary(oryx_mod_bic$mod_2d$bestModel)$cov.robust
rcoefs_oryx_bic <- rmvnorm(1000, coefz_oryx_bic, sigma=covmat_oryx_bic)

pred_bootz_oryx_bic_all <- predict.gamMRSea(object=oryx_mod_bic$mod_2d$bestModel, coeff=rcoefs_oryx_bic[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=oryx_mod_bic$mod_2d$bestModel, coeff=rcoefs_oryx_bic[bt,])
  pred_bootz_oryx_bic_all <- cbind(pred_bootz_oryx_bic_all, pred_boot)
}
pred_ci_oryx_bic <- t(apply(pred_bootz_oryx_bic_all, 1, quantile, probs=c(0.025, 0.975)))


pred_oryx_bic <- predict.gamMRSea(object=oryx_mod_bic$mod_2d$bestModel)
data_oryx_bic <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_oryx_bic, pred_ci_oryx_bic))
colnames(data_oryx_bic) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

graph_pred_oryx_bic <- ggplot(data_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_oryx_bic_lw <- ggplot(data_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_oryx_bic_hi <- ggplot(data_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OryxModel3Pred.png"))
graph_pred_oryx_bic
dev.off()

png(filename=paste0(img_save_dir, "OryxModel3Lw.png"))
graph_pred_oryx_bic_lw
dev.off()

png(filename=paste0(img_save_dir, "OryxModel3Hi.png"))
graph_pred_oryx_bic_hi
dev.off()


pred_bootz_oryx_bic_all_link <- predict.gamMRSea(object=oryx_mod_bic$mod_2d$bestModel, coeff=rcoefs_oryx_bic[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=oryx_mod_bic$mod_2d$bestModel, coeff=rcoefs_oryx_bic[bt,], type='link')
  pred_bootz_oryx_bic_all_link <- cbind(pred_bootz_oryx_bic_all, pred_boot)
}
pred_ci_oryx_bic_link <- t(apply(pred_bootz_oryx_bic_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_oryx_bic_link <- predict.gamMRSea(object=oryx_mod_bic$mod_2d$bestModel)
data_oryx_bic_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_oryx_bic_link, pred_ci_oryx_bic_link))
colnames(data_oryx_bic_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")




# filter out any NAs
testcoefs_oryx_bic <- oryx_mod_bic$mod_2d$bestModel$coefficients
cfmsk_oryx_bic <- !is.na(testcoefs_oryx_bic)
tstcfs_oryx_bic <- testcoefs_oryx_bic[cfmsk_oryx_bic]

sp_col_oryx_bic <- oryx_mod_bic$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_oryx_bic <- str_sub(names(tstcfs_oryx_bic),1,strlenclnam)
coefmask_oryx_bic <- coefnamsplt_oryx_bic == nam2d

# create radial gaussian bases
radii_oryx_bic <- sp_col_oryx_bic$radii
radiiInd_oryx_bic <- sp_col_oryx_bic$radiusIndices
aR_oryx_bic <- sp_col_oryx_bic$knotPos
lrf_oryx_bic <- LRF.g(radiiInd_oryx_bic, distMats$dataDist, radii_oryx_bic, aR_oryx_bic)

# combine coefmask and facts
coefz_oryx_bic <- tstcfs_oryx_bic[coefmask_oryx_bic]
# get predicted values on link scale
predtm_oryx_bic <- lrf_oryx_bic %*% coefz_oryx_bic
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_oryx_bic <- rcoefs_oryx_bic[1,][cfmsk_oryx_bic][coefmask_oryx_bic]
predboot_oryx_bic_2d <- lrf_oryx_bic %*% bootcoefz_oryx_bic
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_oryx_bic[bt,][cfmsk_oryx_bic][coefmask_oryx_bic]
  predbt <- lrf_oryx_bic %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_oryx_bic_2d <- cbind(predboot_oryx_bic_2d, predbt)
}
pred_2d_ci_oryx_bic <- t(apply(predboot_oryx_bic_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_oryx_bic <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_oryx_bic, pred_2d_ci_oryx_bic))
colnames(data_2d_oryx_bic) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


oryx_2d_bic <- ggplot(data_2d_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette =  "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_2d_bic_lw <- ggplot(data_2d_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_2d_bic_hi <- ggplot(data_2d_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_oryx_bic <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_oryx_bic_link - predtm_oryx_bic, pred_ci_oryx_bic_link - pred_2d_ci_oryx_bic))
colnames(data_not2d_oryx_bic) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


oryx_not2d_bic <- ggplot(data_not2d_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_not2d_bic_lw <- ggplot(data_not2d_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_not2d_bic_hi <- ggplot(data_not2d_oryx_bic, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))


oryx_pred_bic_link <- ggplot(data_oryx_bic_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_bic_lw_link <- ggplot(data_oryx_bic_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

oryx_pred_bic_hi_link <- ggplot(data_oryx_bic_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OryxModel32d.png"))
oryx_2d_bic
dev.off()
png(filename=paste0(img_save_dir, "OryxModel32dLw.png"))
oryx_2d_bic_lw
dev.off()
png(filename=paste0(img_save_dir, "OryxModel32dHi.png"))
oryx_2d_bic_hi
dev.off()
png(filename=paste0(img_save_dir, "OryxModel3NotPrd.png"))
oryx_not2d_bic
dev.off()
png(filename=paste0(img_save_dir, "OryxModel3NotLw.png"))
oryx_not2d_bic_lw
dev.off()
png(filename=paste0(img_save_dir, "OryxModel3NotHi.png"))
oryx_not2d_bic_hi
dev.off()
png(filename=paste0(img_save_dir, "OryxModel3Lnk.png"))
oryx_pred_bic_link
dev.off()
png(filename=paste0(img_save_dir, "OryxModel3LnkLw.png"))
oryx_pred_bic_lw_link
dev.off()
png(filename=paste0(img_save_dir, "OryxModel3LnkHi.png"))
oryx_pred_bic_hi_link
dev.off()

## 1d calculations

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
anywater_predvals <- seq(min(sectionz_ndvi$any_water_clip), max(sectionz_ndvi$any_water_clip), length.out=200)
#ndvi, elev, fence, track, any road, wat, dam, any water

dist_med_oryx_bic <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_anywater
index_med_oryx_bic <- which.min(dist_med_oryx_bic)
median_vals_oryx_bic <- sectionz_ndvi[index_med_oryx_bic,]
distMats_pred_oryx_bic <- makeDists(
  cbind(rep(median_vals_oryx_bic$x.pos, 200), rep(median_vals_oryx_bic$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_oryx_bic <- as.data.frame(cbind(pred_oryx_bic, sectionz_ndvi$Elev, sectionz_ndvi$any_water_clip))
colnames(rug_vals_oryx_bic) <- c("predictions", "Elev", "any_water")


predict_oryx_bic <- as.data.frame(cbind(rep(median_vals_oryx_bic$x.pos, 200), 
                                       rep(median_vals_oryx_bic$y.pos, 200), 
                                       rep(median_vals_oryx_bic$Elev, 200),
                                       rep(median_vals_oryx_bic$any_water_clip, 200)))
colnames(predict_oryx_bic) <- c("x.pos", "y.pos", "Elev", "any_water_clip")

predict_elev_oryx_bic <- predict_oryx_bic
predict_elev_oryx_bic$Elev <- elev_predvals
pred_elev_oryx_bic <- predict.gamMRSea(newdata=predict_elev_oryx_bic, object=oryx_mod_bic$mod_2d$bestModel, g2k=distMats_pred_oryx_bic$dataDist)
pred_bootz_oryx_bic_elev <- predict.gamMRSea(newdata=predict_elev_oryx_bic, object=oryx_mod_bic$mod_2d$bestModel, g2k=distMats_pred_oryx_bic$dataDist, coeff=rcoefs_oryx_bic[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_oryx_bic, object=oryx_mod_bic$mod_2d$bestModel, g2k=distMats_pred_oryx_bic$dataDist, coeff=rcoefs_oryx_bic[bt,])
  pred_bootz_oryx_bic_elev <- cbind(pred_bootz_oryx_bic_elev, pred_boot)
}
pred_ci_elev_oryx_bic <- t(apply(pred_bootz_oryx_bic_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_oryx_bic <- cbind(predict_elev_oryx_bic, pred_elev_oryx_bic, pred_ci_elev_oryx_bic)
colnames(pred_df_elev_oryx_bic) <- c(colnames(predict_elev_oryx_bic), "predictions", "LowerCI", "UpperCI")

predict_anywt_oryx_bic <- predict_oryx_bic
predict_anywt_oryx_bic$any_water_clip <- anywater_predvals
pred_anywt_oryx_bic <- predict.gamMRSea(newdata=predict_anywt_oryx_bic, object=oryx_mod_bic$mod_2d$bestModel, g2k=distMats_pred_oryx_bic$dataDist)
pred_bootz_oryx_bic_anywt <- predict.gamMRSea(newdata=predict_anywt_oryx_bic, object=oryx_mod_bic$mod_2d$bestModel, g2k=distMats_pred_oryx_bic$dataDist, coeff=rcoefs_oryx_bic[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_oryx_bic, object=oryx_mod_bic$mod_2d$bestModel, g2k=distMats_pred_oryx_bic$dataDist, coeff=rcoefs_oryx_bic[bt,])
  pred_bootz_oryx_bic_anywt <- cbind(pred_bootz_oryx_bic_anywt, pred_boot)
}
pred_ci_anywt_oryx_bic <- t(apply(pred_bootz_oryx_bic_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_oryx_bic <- cbind(predict_anywt_oryx_bic, pred_anywt_oryx_bic, pred_ci_anywt_oryx_bic)
colnames(pred_df_anywt_oryx_bic) <- c(colnames(predict_anywt_oryx_bic), "predictions", "LowerCI", "UpperCI")




pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_oryx_bc] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")

pt_med_bic <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of median point used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.05, 0.05)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))

plotout_elev_oryx_bic <- ggplot(pred_df_elev_oryx_bic, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bic, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anywt_oryx_bic <- ggplot(pred_df_anywt_oryx_bic, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to any water') + geom_ribbon(aes(x=any_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_oryx_bic, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed3_oryx.png"))
pt_med_bic
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Elev_oryx.png"))
plotout_elev_oryx_bic
dev.off()
png(filename=paste0(img_save_dir, "PtMod3anywater_oryx.png"))
plotout_anywt_oryx_bic
dev.off()




















