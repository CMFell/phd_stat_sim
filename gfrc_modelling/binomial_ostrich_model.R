# This file fits the binomial ostrich model and creates output images
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

sectionz_ndvi$response <- sectionz_ndvi$OstrichPA
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
ostr_2d_mod <- runSALSA2D(
  PA_mod, 
  salsa2dlist_in, 
  d2k=distMats$dataDist, 
  k2k=distMats$knotDist,
  suppress.printout = TRUE
)
an_res_ostr <- Anova(ostr_2d_mod$bestModel, test='F')
rownames(an_res_ostr) <- c(ostr_2d_mod$bestModel$varshortnames, "2dsmooth", "Residuals")

coefz_ostr_sm <- ostr_2d_mod$bestModel$coefficients
covmat_ostr_sm <- summary(ostr_2d_mod$bestModel)$cov.robust
rcoefs_ostr_sm <- rmvnorm(1000, coefz_ostr_sm, sigma=covmat_ostr_sm)

pred_bootz_ostr_sm_all <- predict.gamMRSea(object=ostr_2d_mod$bestModel, coeff=rcoefs_ostr_sm[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=ostr_2d_mod$bestModel, coeff=rcoefs_ostr_sm[bt,])
  pred_bootz_ostr_sm_all <- cbind(pred_bootz_ostr_sm_all, pred_boot)
}
pred_ci_ostr_sm <- t(apply(pred_bootz_ostr_sm_all, 1, quantile, probs=c(0.025, 0.975)))

pred_ostr_sm <- predict.gamMRSea(object=ostr_2d_mod$bestModel)
data_ostr_sm <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_ostr_sm, pred_ci_ostr_sm))
colnames(data_ostr_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_ostr_sm <- cbind(data_ostr_sm, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_ostr_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")

ostr_pred_sm <- ggplot(data_ostr_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
ostr_pred_sm <- ggplot(data_ostr_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ostr_pred_sm_lw <- ggplot(data_ostr_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
ostr_pred_sm_lw <- ggplot(data_ostr_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

ostr_pred_sm_hi <- ggplot(data_ostr_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
ostr_pred_sm_hi  <- ggplot(data_ostr_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OstrModel1Pred.png"))
ostr_pred_sm
dev.off()

png(filename=paste0(img_save_dir, "OstrModel1Lw.png"))
ostr_pred_sm_lw
dev.off()

png(filename=paste0(img_save_dir, "OstrModel1Hi.png"))
ostr_pred_sm_hi
dev.off()






# Second model backward selection
var_list_back <- c("NDVI250", "Elev", "dist_fence_clip", "dist_elec_clip", "dist_road_clip", "dist_track_clip", 
                   "any_road_clip", "dist_wat_clip", "dist_dam_clip", "dist_river_clip", "dist_minor_clip", "any_water_clip")
ostr_mod_bc <- rez_for_var_list(var_list_back, PA_mod, distMats)


coefz_ostr_bc <- ostr_mod_bc$mod_2d$bestModel$coefficients
covmat_ostr_bc <- summary(ostr_mod_bc$mod_2d$bestModel)$cov.robust
rcoefs_ostr_bc <- rmvnorm(1000, coefz_ostr_bc, sigma=covmat_ostr_bc)

pred_bootz_ostr_bc_all <- predict.gamMRSea(object=ostr_mod_bc$mod_2d$bestModel, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=ostr_mod_bc$mod_2d$bestModel, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_all <- cbind(pred_bootz_ostr_bc_all, pred_boot)
}
pred_ci_ostr_bc <- t(apply(pred_bootz_ostr_bc_all, 1, quantile, probs=c(0.025, 0.975)))

pred_ostr_bc <- predict.gamMRSea(object=ostr_mod_bc$mod_2d$bestModel)
data_ostr_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_ostr_bc, pred_ci_ostr_bc))
colnames(data_ostr_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_ostr_bc <- cbind(data_ostr_bc, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_ostr_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")

graph_pred_ostr_bc <- ggplot(data_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_ostr_bc <- ggplot(data_ostr_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_ostr_bc_lw <- ggplot(data_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_ostr_bc_lw <- ggplot(data_ostr_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_ostr_bc_hi <- ggplot(data_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_ostr_bc_hi <- ggplot(data_ostr_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OstrModel2Pred.png"))
graph_pred_ostr_bc
dev.off()

png(filename=paste0(img_save_dir, "OstrModel2Lw.png"))
graph_pred_ostr_bc_lw
dev.off()

png(filename=paste0(img_save_dir, "OstrModel2Hi.png"))
graph_pred_ostr_bc_hi
dev.off()

### Not included in thesis graphs not updated to utm coords

pred_bootz_ostr_bc_all_link <- predict.gamMRSea(object=ostr_mod_bc$mod_2d$bestModel, coeff=rcoefs_ostr_bc[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=ostr_mod_bc$mod_2d$bestModel, coeff=rcoefs_ostr_bc[bt,], type='link')
  pred_bootz_ostr_bc_all_link <- cbind(pred_bootz_ostr_bc_all, pred_boot)
}
pred_ci_ostr_bc_link <- t(apply(pred_bootz_ostr_bc_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_ostr_bc_link <- predict.gamMRSea(object=ostr_mod_bc$mod_2d$bestModel)
data_ostr_bc_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_ostr_bc_link, pred_ci_ostr_bc_link))
colnames(data_ostr_bc_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


# filter out any NAs
testcoefs_ostr_bc <- ostr_mod_bc$mod_2d$bestModel$coefficients
cfmsk_ostr_bc <- !is.na(testcoefs_ostr_bc)
tstcfs_ostr_bc <- testcoefs_ostr_bc[cfmsk_ostr_bc]

sp_col_ostr_bc <- ostr_mod_bc$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_ostr_bc <- str_sub(names(tstcfs_ostr_bc),1,strlenclnam)
coefmask_ostr_bc <- coefnamsplt_ostr_bc == nam2d

# create radial gaussian bases
radii_ostr_bc <- sp_col_ostr_bc$radii
radiiInd_ostr_bc <- sp_col_ostr_bc$radiusIndices
aR_ostr_bc <- sp_col_ostr_bc$knotPos
lrf_ostr_bc <- LRF.g(radiiInd_ostr_bc, distMats$dataDist, radii_ostr_bc, aR_ostr_bc)

# combine coefmask and facts
coefz_ostr_bc <- tstcfs_ostr_bc[coefmask_ostr_bc]
# get predicted values on link scale
predtm_ostr_bc <- lrf_ostr_bc %*% coefz_ostr_bc
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_ostr_bc <- rcoefs_ostr_bc[1,][cfmsk_ostr_bc][coefmask_ostr_bc]
predboot_ostr_bc_2d <- lrf_ostr_bc %*% bootcoefz_ostr_bc
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_ostr_bc[bt,][cfmsk_ostr_bc][coefmask_ostr_bc]
  predbt <- lrf_ostr_bc %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_ostr_bc_2d <- cbind(predboot_ostr_bc_2d, predbt)
}
pred_2d_ci_ostr_bc <- t(apply(predboot_ostr_bc_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_ostr_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_ostr_bc, pred_2d_ci_ostr_bc))
colnames(data_2d_ostr_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


ostr_2d_bc <- ggplot(data_2d_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

ostr_2d_bc_lw <- ggplot(data_2d_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

ostr_2d_bc_hi <- ggplot(data_2d_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_ostr_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_ostr_bc_link - predtm_ostr_bc, pred_ci_ostr_bc_link - pred_2d_ci_ostr_bc))
colnames(data_not2d_ostr_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


ostr_not2d_bc <- ggplot(data_not2d_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

ostr_not2d_bc_lw <- ggplot(data_not2d_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

ostr_not2d_bc_hi <- ggplot(data_not2d_ostr_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))


ostr_pred_bc_link <- ggplot(data_ostr_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

ostr_pred_bc_lw_link <- ggplot(data_ostr_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

ostr_pred_bc_hi_link <- ggplot(data_ostr_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "OstrModel22d.png"))
  ostr_2d_bc
dev.off()
png(filename=paste0(img_save_dir, "OstrModel22dLw.png"))
  ostr_2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "OstrModel22dHi.png"))
  ostr_2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "OstrModel2NotPrd.png"))
  ostr_not2d_bc
dev.off()
png(filename=paste0(img_save_dir, "OstrModel2NotLw.png"))
  ostr_not2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "OstrModel2NotHi.png"))
  ostr_not2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "OstrModel2Lnk.png"))
ostr_pred_bc_link
dev.off()
png(filename=paste0(img_save_dir, "OstrModel2LnkLw.png"))
ostr_pred_bc_lw_link
dev.off()
png(filename=paste0(img_save_dir, "OstrModel2LnkHi.png"))
ostr_pred_bc_hi_link
dev.off()

### Thesis included below

## 1d calculations

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
ndvi_predvals <- seq(min(sectionz_ndvi$NDVI250), max(sectionz_ndvi$NDVI250), length.out=200)
fence_predvals <- seq(min(sectionz_ndvi$dist_fence_clip), max(sectionz_ndvi$dist_fence_clip), length.out=200)
elec_predvals <- seq(min(sectionz_ndvi$dist_elec_clip), max(sectionz_ndvi$dist_elec_clip), length.out=200)
road_predvals <- seq(min(sectionz_ndvi$dist_road_clip), max(sectionz_ndvi$dist_road_clip), length.out=200)
track_predvals <- seq(min(sectionz_ndvi$dist_track_clip), max(sectionz_ndvi$dist_track_clip), length.out=200)
anyrd_predvals <- seq(min(sectionz_ndvi$any_road_clip), max(sectionz_ndvi$any_road_clip), length.out=200)
wat_predvals <- seq(min(sectionz_ndvi$dist_wat_clip), max(sectionz_ndvi$dist_wat_clip), length.out=200)
dam_predvals <- seq(min(sectionz_ndvi$dist_dam_clip), max(sectionz_ndvi$dist_dam_clip), length.out=200)
river_predvals <- seq(min(sectionz_ndvi$dist_river_clip), max(sectionz_ndvi$dist_river_clip), length.out=200)
minor_predvals <- seq(min(sectionz_ndvi$dist_minor_clip), max(sectionz_ndvi$dist_minor_clip), length.out=200)
anywater_predvals <- seq(min(sectionz_ndvi$any_water_clip), max(sectionz_ndvi$any_water_clip), length.out=200)


# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
elev_normal <- (sectionz_ndvi$Elev - mean(sectionz_ndvi$Elev)) / sd(sectionz_ndvi$Elev)
ndvi_normal <- (sectionz_ndvi$NDVI250 - mean(sectionz_ndvi$NDVI250)) / sd(sectionz_ndvi$NDVI250)
fence_normal <- (sectionz_ndvi$dist_fence_clip - mean(sectionz_ndvi$dist_fence_clip)) / sd(sectionz_ndvi$dist_fence_clip)
elec_normal <- (sectionz_ndvi$dist_elec_clip - mean(sectionz_ndvi$dist_elec_clip)) / sd(sectionz_ndvi$dist_elec_clip)
road_normal <- (sectionz_ndvi$dist_road_clip - mean(sectionz_ndvi$dist_road_clip)) / sd(sectionz_ndvi$dist_road_clip)
track_normal <- (sectionz_ndvi$dist_track_clip - mean(sectionz_ndvi$dist_track_clip)) / sd(sectionz_ndvi$dist_track_clip)
anyrd_normal <- (sectionz_ndvi$any_road_clip - mean(sectionz_ndvi$any_road_clip)) / sd(sectionz_ndvi$any_road_clip)
wat_normal <- (sectionz_ndvi$dist_wat_clip - mean(sectionz_ndvi$dist_wat_clip)) / sd(sectionz_ndvi$dist_wat_clip)
dam_normal <- (sectionz_ndvi$dist_dam_clip - mean(sectionz_ndvi$dist_dam_clip)) / sd(sectionz_ndvi$dist_dam_clip)
river_normal <- (sectionz_ndvi$dist_river_clip - mean(sectionz_ndvi$dist_river_clip)) / sd(sectionz_ndvi$dist_river_clip)
minor_normal <- (sectionz_ndvi$dist_minor_clip - mean(sectionz_ndvi$dist_minor_clip)) / sd(sectionz_ndvi$dist_minor_clip)
anywater_normal <- (sectionz_ndvi$any_water_clip - mean(sectionz_ndvi$any_water_clip)) / sd(sectionz_ndvi$any_water_clip)


dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_elev <- (elev_normal - median(elev_normal))^2
dist_med_ndvi <- (ndvi_normal - median(ndvi_normal))^2
dist_med_fence <- (fence_normal - median(fence_normal))^2
dist_med_elec <- (elec_normal - median(elec_normal))^2
dist_med_road <- (road_normal - median(road_normal))^2
dist_med_track <- (track_normal - median(track_normal))^2
dist_med_anyrd <- (anyrd_normal - median(anyrd_normal))^2
dist_med_wat <- (wat_normal - median(wat_normal))^2
dist_med_dam <- (dam_normal - median(dam_normal))^2
dist_med_river <- (river_normal - median(river_normal))^2
dist_med_minor <- (minor_normal - median(minor_normal))^2
dist_med_anywater <- (anywater_normal - median(anywater_normal))^2


dist_med_ostr_bc <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_ndvi + dist_med_fence + dist_med_elec + 
  dist_med_road + dist_med_track + dist_med_anyrd + dist_med_wat + dist_med_dam + dist_med_river + 
  dist_med_minor + dist_med_anywater
index_med_ostr_bc <- which.min(dist_med_ostr_bc)
median_vals_ostr_bc <- sectionz_ndvi[index_med_ostr_bc,]
distMats_pred_ostr_bc <- makeDists(
  cbind(rep(median_vals_ostr_bc$x.pos, 200), rep(median_vals_ostr_bc$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_ostr_bc <- as.data.frame(cbind(pred_ostr_bc, sectionz_ndvi$Elev, sectionz_ndvi$NDVI250, 
                                   sectionz_ndvi$dist_elec_clip, sectionz_ndvi$dist_fence_clip, 
                                   sectionz_ndvi$dist_road_clip, sectionz_ndvi$dist_track_clip, 
                                   sectionz_ndvi$any_road_clip, sectionz_ndvi$dist_wat_clip,
                                   sectionz_ndvi$dist_dam_clip, sectionz_ndvi$any_water_clip,
                                   sectionz_ndvi$dist_minor_clip, sectionz_ndvi$dist_river_clip))
colnames(rug_vals_ostr_bc) <- c("predictions", "Elev", "NDVI250", "dist_elec", "dist_fence", "dist_road", "dist_track",
                               "any_road", "dist_wat", "dist_dam", "any_water", "dist_minor", "dist_river")
rug_vals_ostr_bc$dist_elec_km <- rug_vals_ostr_bc$dist_elec / 1000
rug_vals_ostr_bc$dist_fence_km <- rug_vals_ostr_bc$dist_fence / 1000
rug_vals_ostr_bc$dist_road_km <- rug_vals_ostr_bc$dist_road / 1000
rug_vals_ostr_bc$dist_track_km <- rug_vals_ostr_bc$dist_track / 1000
rug_vals_ostr_bc$any_road_km <- rug_vals_ostr_bc$any_road / 1000
rug_vals_ostr_bc$dist_wat_km <- rug_vals_ostr_bc$dist_wat / 1000
rug_vals_ostr_bc$dist_dam_km <- rug_vals_ostr_bc$dist_dam / 1000
rug_vals_ostr_bc$any_water_km <- rug_vals_ostr_bc$any_water / 1000
rug_vals_ostr_bc$dist_minor_km <- rug_vals_ostr_bc$dist_minor / 1000
rug_vals_ostr_bc$dist_river_km <- rug_vals_ostr_bc$dist_river / 1000


predict_ostr_bc <- as.data.frame(cbind(rep(median_vals_ostr_bc$x.pos, 200), 
                                          rep(median_vals_ostr_bc$y.pos, 200), 
                                          rep(median_vals_ostr_bc$Elev, 200),
                                          rep(median_vals_ostr_bc$NDVI250, 200),
                                          rep(median_vals_ostr_bc$dist_elec_clip, 200),
                                          rep(median_vals_ostr_bc$dist_fence_clip, 200),
                                          rep(median_vals_ostr_bc$dist_road_clip, 200),
                                          rep(median_vals_ostr_bc$dist_track_clip, 200),
                                          rep(median_vals_ostr_bc$any_road_clip, 200),
                                          rep(median_vals_ostr_bc$dist_wat_clip, 200),
                                          rep(median_vals_ostr_bc$dist_dam_clip, 200),
                                          rep(median_vals_ostr_bc$any_water_clip, 200),
                                          rep(median_vals_ostr_bc$dist_minor_clip, 200),
                                          rep(median_vals_ostr_bc$dist_river_clip, 200)))
colnames(predict_ostr_bc) <- c("x.pos", "y.pos", "Elev", "NDVI250", "dist_elec_clip", "dist_fence_clip", 
                               "dist_road_clip", "dist_track_clip", "any_road_clip",
                               "dist_wat_clip", "dist_dam_clip", "any_water_clip", "dist_minor_clip", "dist_river_clip")

predict_elev_ostr_bc <- predict_ostr_bc
predict_elev_ostr_bc$Elev <- elev_predvals
pred_elev_ostr_bc <- predict.gamMRSea(newdata=predict_elev_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_elev <- predict.gamMRSea(newdata=predict_elev_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_elev <- cbind(pred_bootz_ostr_bc_elev, pred_boot)
}
pred_ci_elev_ostr_bc <- t(apply(pred_bootz_ostr_bc_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_ostr_bc <- cbind(predict_elev_ostr_bc, pred_elev_ostr_bc, pred_ci_elev_ostr_bc)
colnames(pred_df_elev_ostr_bc) <- c(colnames(predict_elev_ostr_bc), "predictions", "LowerCI", "UpperCI")

predict_ndvi_ostr_bc <- predict_ostr_bc
predict_ndvi_ostr_bc$NDVI250 <- ndvi_predvals
pred_ndvi_ostr_bc <- predict.gamMRSea(newdata=predict_ndvi_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_ndvi <- predict.gamMRSea(newdata=predict_ndvi_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_ndvi_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_ndvi <- cbind(pred_bootz_ostr_bc_ndvi, pred_boot)
}
pred_ci_ndvi_ostr_bc <- t(apply(pred_bootz_ostr_bc_ndvi, 1, quantile, probs=c(0.025, 0.975)))
pred_df_ndvi_ostr_bc <- cbind(predict_ndvi_ostr_bc, pred_ndvi_ostr_bc, pred_ci_ndvi_ostr_bc)
colnames(pred_df_ndvi_ostr_bc) <- c(colnames(predict_ndvi_ostr_bc), "predictions", "LowerCI", "UpperCI")

predict_fence_ostr_bc <- predict_ostr_bc
predict_fence_ostr_bc$dist_fence_clip <- fence_predvals
pred_fence_ostr_bc <- predict.gamMRSea(newdata=predict_fence_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_fence <- predict.gamMRSea(newdata=predict_fence_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_fence_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_fence <- cbind(pred_bootz_ostr_bc_fence, pred_boot)
}
pred_ci_fence_ostr_bc <- t(apply(pred_bootz_ostr_bc_fence, 1, quantile, probs=c(0.025, 0.975)))
pred_df_fence_ostr_bc <- cbind(predict_fence_ostr_bc, pred_fence_ostr_bc, pred_ci_fence_ostr_bc)
colnames(pred_df_fence_ostr_bc) <- c(colnames(predict_fence_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_ostr_bc$dist_fence_km <- pred_df_fence_ostr_bc$dist_fence_clip / 1000

predict_elec_ostr_bc <- predict_ostr_bc
predict_elec_ostr_bc$dist_elec_clip <- elec_predvals
pred_elec_ostr_bc <- predict.gamMRSea(newdata=predict_elec_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_elec <- predict.gamMRSea(newdata=predict_elec_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elec_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_elec <- cbind(pred_bootz_ostr_bc_elec, pred_boot)
}
pred_ci_elec_ostr_bc <- t(apply(pred_bootz_ostr_bc_elec, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elec_ostr_bc <- cbind(predict_elec_ostr_bc, pred_elec_ostr_bc, pred_ci_elec_ostr_bc)
colnames(pred_df_elec_ostr_bc) <- c(colnames(predict_elec_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_elec_ostr_bc$dist_elec_km <- pred_df_elec_ostr_bc$dist_elec_clip / 1000

predict_road_ostr_bc <- predict_ostr_bc
predict_road_ostr_bc$dist_road_clip <- road_predvals
pred_road_ostr_bc <- predict.gamMRSea(newdata=predict_road_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_road <- predict.gamMRSea(newdata=predict_road_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_road_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_road <- cbind(pred_bootz_ostr_bc_road, pred_boot)
}
pred_ci_road_ostr_bc <- t(apply(pred_bootz_ostr_bc_road, 1, quantile, probs=c(0.025, 0.975)))
pred_df_road_ostr_bc <- cbind(predict_road_ostr_bc, pred_road_ostr_bc, pred_ci_road_ostr_bc)
colnames(pred_df_road_ostr_bc) <- c(colnames(predict_road_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_road_ostr_bc$dist_road_km <- pred_df_road_ostr_bc$dist_road_clip / 1000

predict_track_ostr_bc <- predict_ostr_bc
predict_track_ostr_bc$dist_track_clip <- track_predvals
pred_track_ostr_bc <- predict.gamMRSea(newdata=predict_track_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_track <- predict.gamMRSea(newdata=predict_track_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_track <- cbind(pred_bootz_ostr_bc_track, pred_boot)
}
pred_ci_track_ostr_bc <- t(apply(pred_bootz_ostr_bc_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_ostr_bc <- cbind(predict_track_ostr_bc, pred_track_ostr_bc, pred_ci_track_ostr_bc)
colnames(pred_df_track_ostr_bc) <- c(colnames(predict_track_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_track_ostr_bc$dist_track_km <- pred_df_track_ostr_bc$dist_track_clip / 1000

predict_anyrd_ostr_bc <- predict_ostr_bc
predict_anyrd_ostr_bc$any_road_clip <- anyrd_predvals
pred_anyrd_ostr_bc <- predict.gamMRSea(newdata=predict_anyrd_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_anyrd <- predict.gamMRSea(newdata=predict_anyrd_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_anyrd <- cbind(pred_bootz_ostr_bc_anyrd, pred_boot)
}
pred_ci_anyrd_ostr_bc <- t(apply(pred_bootz_ostr_bc_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_ostr_bc <- cbind(predict_anyrd_ostr_bc, pred_anyrd_ostr_bc, pred_ci_anyrd_ostr_bc)
colnames(pred_df_anyrd_ostr_bc) <- c(colnames(predict_anyrd_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anyrd_ostr_bc$any_road_km <- pred_df_anyrd_ostr_bc$any_road_clip / 1000

predict_wat_ostr_bc <- predict_ostr_bc
predict_wat_ostr_bc$dist_wat_clip <- wat_predvals
pred_wat_ostr_bc <- predict.gamMRSea(newdata=predict_wat_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_wat <- predict.gamMRSea(newdata=predict_wat_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_wat_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_wat <- cbind(pred_bootz_ostr_bc_wat, pred_boot)
}
pred_ci_wat_ostr_bc <- t(apply(pred_bootz_ostr_bc_wat, 1, quantile, probs=c(0.025, 0.975)))
pred_df_wat_ostr_bc <- cbind(predict_wat_ostr_bc, pred_wat_ostr_bc, pred_ci_wat_ostr_bc)
colnames(pred_df_wat_ostr_bc) <- c(colnames(predict_wat_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_wat_ostr_bc$dist_wat_km <- pred_df_wat_ostr_bc$dist_wat_clip / 1000

predict_dam_ostr_bc <- predict_ostr_bc
predict_dam_ostr_bc$dist_dam_clip <- dam_predvals
pred_dam_ostr_bc <- predict.gamMRSea(newdata=predict_dam_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_dam <- predict.gamMRSea(newdata=predict_dam_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_dam_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_dam <- cbind(pred_bootz_ostr_bc_dam, pred_boot)
}
pred_ci_dam_ostr_bc <- t(apply(pred_bootz_ostr_bc_dam, 1, quantile, probs=c(0.025, 0.975)))
pred_df_dam_ostr_bc <- cbind(predict_dam_ostr_bc, pred_dam_ostr_bc, pred_ci_dam_ostr_bc)
colnames(pred_df_dam_ostr_bc) <- c(colnames(predict_dam_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_dam_ostr_bc$dist_dam_km <- pred_df_dam_ostr_bc$dist_dam_clip / 1000

predict_anywt_ostr_bc <- predict_ostr_bc
predict_anywt_ostr_bc$any_water_clip <- anywater_predvals
pred_anywt_ostr_bc <- predict.gamMRSea(newdata=predict_anywt_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_anywt <- predict.gamMRSea(newdata=predict_anywt_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_anywt <- cbind(pred_bootz_ostr_bc_anywt, pred_boot)
}
pred_ci_anywt_ostr_bc <- t(apply(pred_bootz_ostr_bc_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_ostr_bc <- cbind(predict_anywt_ostr_bc, pred_anywt_ostr_bc, pred_ci_anywt_ostr_bc)
colnames(pred_df_anywt_ostr_bc) <- c(colnames(predict_anywt_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anywt_ostr_bc$any_water_km <- pred_df_anywt_ostr_bc$any_water_clip / 1000

predict_minor_ostr_bc <- predict_ostr_bc
predict_minor_ostr_bc$dist_minor_clip <- minor_predvals
pred_minor_ostr_bc <- predict.gamMRSea(newdata=predict_minor_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_minor <- predict.gamMRSea(newdata=predict_minor_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_minor_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_minor <- cbind(pred_bootz_ostr_bc_minor, pred_boot)
}
pred_ci_minor_ostr_bc <- t(apply(pred_bootz_ostr_bc_minor, 1, quantile, probs=c(0.025, 0.975)))
pred_df_minor_ostr_bc <- cbind(predict_minor_ostr_bc, pred_minor_ostr_bc, pred_ci_minor_ostr_bc)
colnames(pred_df_minor_ostr_bc) <- c(colnames(predict_minor_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_minor_ostr_bc$dist_minor_km <- pred_df_minor_ostr_bc$dist_minor_clip / 1000

predict_river_ostr_bc <- predict_ostr_bc
predict_river_ostr_bc$dist_river_clip <- river_predvals
pred_river_ostr_bc <- predict.gamMRSea(newdata=predict_river_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist)
pred_bootz_ostr_bc_river <- predict.gamMRSea(newdata=predict_river_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_river_ostr_bc, object=ostr_mod_bc$mod_2d$bestModel, g2k=distMats_pred_ostr_bc$dataDist, coeff=rcoefs_ostr_bc[bt,])
  pred_bootz_ostr_bc_river <- cbind(pred_bootz_ostr_bc_river, pred_boot)
}
pred_ci_river_ostr_bc <- t(apply(pred_bootz_ostr_bc_river, 1, quantile, probs=c(0.025, 0.975)))
pred_df_river_ostr_bc <- cbind(predict_river_ostr_bc, pred_river_ostr_bc, pred_ci_river_ostr_bc)
colnames(pred_df_river_ostr_bc) <- c(colnames(predict_river_ostr_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_ostr_bc$dist_river_km <- pred_df_river_ostr_bc$dist_river_clip / 1000


pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_ostr_bc] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")
df_pointz <- cbind(df_pointz, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(df_pointz) <- c("Lat", "Long", "Points", "UTMX", "UTMY", "distX", "distY")

pt_med2 <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of median point used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.05, 0.05)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))
pt_med2 <- ggplot(df_pointz, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(6, 6)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))

plotout_elev_ostr_bc <- ggplot(pred_df_elev_ostr_bc, aes(x=Elev, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Elevation (m)') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_ndvi_ostr_bc <- ggplot(pred_df_ndvi_ostr_bc, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_fence_ostr_bc <- ggplot(pred_df_fence_ostr_bc, aes(x=dist_fence_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to fence (km)') + geom_ribbon(aes(x=dist_fence_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_elec_ostr_bc <- ggplot(pred_df_elec_ostr_bc, aes(x=dist_elec_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to electricity line (km)') + geom_ribbon(aes(x=dist_elec_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc'))  + theme(text = element_text(size=20))
plotout_road_ostr_bc <- ggplot(pred_df_road_ostr_bc, aes(x=dist_road_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to road (km)') + geom_ribbon(aes(x=dist_road_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc'))  + theme(text = element_text(size=20))
plotout_track_ostr_bc <- ggplot(pred_df_track_ostr_bc, aes(x=dist_track_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to track (km)') + geom_ribbon(aes(x=dist_track_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc'))  + theme(text = element_text(size=20))
plotout_anyrd_ostr_bc <- ggplot(pred_df_anyrd_ostr_bc, aes(x=any_road_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any road (km)') + geom_ribbon(aes(x=any_road_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_wat_ostr_bc <- ggplot(pred_df_wat_ostr_bc, aes(x=dist_wat_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to water hole (km)') + geom_ribbon(aes(x=dist_wat_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_dam_ostr_bc <- ggplot(pred_df_dam_ostr_bc, aes(x=dist_dam_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to dam (km)') + geom_ribbon(aes(x=dist_dam_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anywt_ostr_bc <- ggplot(pred_df_anywt_ostr_bc, aes(x=any_water_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any water (km)') + geom_ribbon(aes(x=any_water_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_minor_ostr_bc <- ggplot(pred_df_minor_ostr_bc, aes(x=dist_minor_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to minor river (km)') + geom_ribbon(aes(x=dist_minor_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_river_ostr_bc <- ggplot(pred_df_river_ostr_bc, aes(x=dist_river_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to river (km)') + geom_ribbon(aes(x=dist_river_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_ostr_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed2_ostr.png"))
pt_med2
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Elev_ostr.png"))
plotout_elev_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2NDVI_ostr.png"))
plotout_ndvi_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Fence_ostr.png"))
plotout_fence_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Elec_ostr.png"))
plotout_elec_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Road_ostr.png"))
plotout_road_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Track_ostr.png"))
plotout_track_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Anyrd_ostr.png"))
plotout_anyrd_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2wat_ostr.png"))
plotout_wat_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2dam_ostr.png"))
plotout_dam_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2anywater_ostr.png"))
plotout_anywt_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2minor_ostr.png"))
plotout_minor_ostr_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2river_ostr.png"))
plotout_river_ostr_bc
dev.off()




