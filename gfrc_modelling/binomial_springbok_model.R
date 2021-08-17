# This file fits the binomial springbok model and creates output images
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

sectionz_ndvi$response <- sectionz_ndvi$SpringbokPA
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
bok_2d_mod <- runSALSA2D(
  PA_mod, 
  salsa2dlist_in, 
  d2k=distMats$dataDist, 
  k2k=distMats$knotDist,
  suppress.printout = TRUE
)
an_res_bok <- Anova(bok_2d_mod$bestModel, test='F')
rownames(an_res_bok) <- c(bok_2d_mod$bestModel$varshortnames, "2dsmooth", "Residuals")

coefz_bok_sm <- bok_2d_mod$bestModel$coefficients
covmat_bok_sm <- summary(bok_2d_mod$bestModel)$cov.robust
rcoefs_bok_sm <- rmvnorm(1000, coefz_bok_sm, sigma=covmat_bok_sm)

pred_bootz_bok_sm_all <- predict.gamMRSea(object=bok_2d_mod$bestModel, coeff=rcoefs_bok_sm[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=bok_2d_mod$bestModel, coeff=rcoefs_bok_sm[bt,])
  pred_bootz_bok_sm_all <- cbind(pred_bootz_bok_sm_all, pred_boot)
}
pred_ci_bok_sm <- t(apply(pred_bootz_bok_sm_all, 1, quantile, probs=c(0.025, 0.975)))

pred_bok_sm <- predict.gamMRSea(object=bok_2d_mod$bestModel)
data_bok_sm <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_bok_sm, pred_ci_bok_sm))
colnames(data_bok_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

data_bok_sm <- as.data.frame(cbind(data_bok_sm, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_bok_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")


bok_pred_sm <- ggplot(data_bok_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))
bok_pred_sm <- ggplot(data_bok_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

bok_pred_sm_lw <- ggplot(data_bok_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))
bok_pred_sm_lw <- ggplot(data_bok_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

bok_pred_sm_hi <- ggplot(data_bok_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))
bok_pred_sm_hi <- ggplot(data_bok_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "BokModel1Pred.png"))
bok_pred_sm
dev.off()

png(filename=paste0(img_save_dir, "BokModel1Lw.png"))
bok_pred_sm_lw
dev.off()

png(filename=paste0(img_save_dir, "BokModel1Hi.png"))
bok_pred_sm_hi
dev.off()

# Second model backward selection (forward selection gives the same)
var_list_back <- c("NDVI250", "Elev", "dist_fence_clip", "dist_elec_clip","dist_track_clip", 
                   "any_road_clip", "dist_wat_clip", "dist_dam_clip", "dist_river_clip", "any_water_clip")
bok_mod_bc <- rez_for_var_list(var_list_back, PA_mod, distMats)


coefz_bok_bc <- bok_mod_bc$mod_2d$bestModel$coefficients
covmat_bok_bc <- summary(bok_mod_bc$mod_2d$bestModel)$cov.robust
rcoefs_bok_bc <- rmvnorm(1000, coefz_bok_bc, sigma=covmat_bok_bc)

pred_bootz_bok_bc_all <- predict.gamMRSea(object=bok_mod_bc$mod_2d$bestModel, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=bok_mod_bc$mod_2d$bestModel, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_all <- cbind(pred_bootz_bok_bc_all, pred_boot)
}
pred_ci_bok_bc <- t(apply(pred_bootz_bok_bc_all, 1, quantile, probs=c(0.025, 0.975)))


pred_bok_bc <- predict.gamMRSea(object=bok_mod_bc$mod_2d$bestModel)
data_bok_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_bok_bc, pred_ci_bok_bc))
colnames(data_bok_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

data_bok_bc <- as.data.frame(cbind(data_bok_bc, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_bok_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")


graph_pred_bok_bc <- ggplot(data_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_bok_bc <- ggplot(data_bok_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_bok_bc_lw <- ggplot(data_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_bok_bc_lw <- ggplot(data_bok_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_bok_bc_hi <- ggplot(data_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_bok_bc_hi <- ggplot(data_bok_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "BokModel2Pred.png"))
graph_pred_bok_bc
dev.off()

png(filename=paste0(img_save_dir, "BokModel2Lw.png"))
graph_pred_bok_bc_lw
dev.off()

png(filename=paste0(img_save_dir, "BokModel2Hi.png"))
graph_pred_bok_bc_hi
dev.off()

### ot included in final thesis not updated

pred_bootz_bok_bc_all_link <- predict.gamMRSea(object=bok_mod_bc$mod_2d$bestModel, coeff=rcoefs_bok_bc[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=bok_mod_bc$mod_2d$bestModel, coeff=rcoefs_bok_bc[bt,], type='link')
  pred_bootz_bok_bc_all_link <- cbind(pred_bootz_bok_bc_all, pred_boot)
}
pred_ci_bok_bc_link <- t(apply(pred_bootz_bok_bc_all_link, 1, quantile, probs=c(0.025, 0.975)))

pred_bok_bc_link <- predict.gamMRSea(object=bok_mod_bc$mod_2d$bestModel, type='link')
data_bok_bc_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_bok_bc_link, pred_ci_bok_bc_link))
colnames(data_bok_bc_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


# filter out any NAs
testcoefs_bok_bc <- bok_mod_bc$mod_2d$bestModel$coefficients
cfmsk_bok_bc <- !is.na(testcoefs_bok_bc)
tstcfs_bok_bc <- testcoefs_bok_bc[cfmsk_bok_bc]

sp_col_bok_bc <- bok_mod_bc$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_bok_bc <- str_sub(names(tstcfs_bok_bc),1,strlenclnam)
coefmask_bok_bc <- coefnamsplt_bok_bc == nam2d

# create radial gaussian bases
radii_bok_bc <- sp_col_bok_bc$radii
radiiInd_bok_bc <- sp_col_bok_bc$radiusIndices
aR_bok_bc <- sp_col_bok_bc$knotPos
lrf_bok_bc <- LRF.g(radiiInd_bok_bc, distMats$dataDist, radii_bok_bc, aR_bok_bc)

# combine coefmask and facts
coefz_bok_bc <- tstcfs_bok_bc[coefmask_bok_bc]
# get predicted values on link scale
predtm_bok_bc <- lrf_bok_bc %*% coefz_bok_bc
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_bok_bc <- rcoefs_bok_bc[1,][cfmsk_bok_bc][coefmask_bok_bc]
predboot_bok_bc_2d <- lrf_bok_bc %*% bootcoefz_bok_bc
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_bok_bc[bt,][cfmsk_bok_bc][coefmask_bok_bc]
  predbt <- lrf_bok_bc %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_bok_bc_2d <- cbind(predboot_bok_bc_2d, predbt)
}
pred_2d_ci_bok_bc <- t(apply(predboot_bok_bc_2d, 1, quantile, probs=c(0.025, 0.975)))

data_lnk_bok_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_bok_bc_link, pred_ci_bok_bc_link))
colnames(data_lnk_bok_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

bok_lnk_bc <- ggplot(data_lnk_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10, 0), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))

bok_lnk_bc_lw <- ggplot(data_lnk_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))

bok_lnk_bc_hi <- ggplot(data_lnk_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))

data_2d_bok_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_bok_bc, pred_2d_ci_bok_bc))
colnames(data_2d_bok_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

bok_2d_bc <- ggplot(data_2d_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))

bok_2d_bc_lw <- ggplot(data_2d_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))

bok_2d_bc_hi <- ggplot(data_2d_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, 
                       guide=guide_colorbar(title=NULL)) +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_bok_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_bok_bc_link - predtm_bok_bc, pred_ci_bok_bc_link - pred_2d_ci_bok_bc))
colnames(data_not2d_bok_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


bok_not2d_bc <- ggplot(data_not2d_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

bok_not2d_bc_lw <- ggplot(data_not2d_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

bok_not2d_bc_hi <- ggplot(data_not2d_bok_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))


bok_pred_bc_link <- ggplot(data_bok_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

bok_pred_bc_lw_link <- ggplot(data_bok_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

bok_pred_bc_hi_link <- ggplot(data_bok_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,0), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "BokModel22d.png"))
  bok_2d_bc
dev.off()
png(filename=paste0(img_save_dir, "BokModel22dLw.png"))
  bok_2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "BokModel22dHi.png"))
  bok_2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "BokModel2NotPrd.png"))
  bok_not2d_bc
dev.off()
png(filename=paste0(img_save_dir, "BokModel2NotLw.png"))
  bok_not2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "BokModel2NotHi.png"))
  bok_not2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "BokModel2Lnk.png"))
bok_pred_bc_link
dev.off()
png(filename=paste0(img_save_dir, "BokModel2LnkLw.png"))
bok_pred_bc_lw_link
dev.off()
png(filename=paste0(img_save_dir, "BokModel2LnkHi.png"))
bok_pred_bc_hi_link
dev.off()

### included in thesis below

## 1d calculations

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
ndvi_predvals <- seq(min(sectionz_ndvi$NDVI250), max(sectionz_ndvi$NDVI250), length.out=200)
fence_predvals <- seq(min(sectionz_ndvi$dist_fence_clip), max(sectionz_ndvi$dist_fence_clip), length.out=200)
elec_predvals <- seq(min(sectionz_ndvi$dist_elec_clip), max(sectionz_ndvi$dist_elec_clip), length.out=200)
track_predvals <- seq(min(sectionz_ndvi$dist_track_clip), max(sectionz_ndvi$dist_track_clip), length.out=200)
anyrd_predvals <- seq(min(sectionz_ndvi$any_road_clip), max(sectionz_ndvi$any_road_clip), length.out=200)
wat_predvals <- seq(min(sectionz_ndvi$dist_wat_clip), max(sectionz_ndvi$dist_wat_clip), length.out=200)
dam_predvals <- seq(min(sectionz_ndvi$dist_dam_clip), max(sectionz_ndvi$dist_dam_clip), length.out=200)
river_predvals <- seq(min(sectionz_ndvi$dist_river_clip), max(sectionz_ndvi$dist_river_clip), length.out=200)
anywater_predvals <- seq(min(sectionz_ndvi$any_water_clip), max(sectionz_ndvi$any_water_clip), length.out=200)


# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
elev_normal <- (sectionz_ndvi$Elev - mean(sectionz_ndvi$Elev)) / sd(sectionz_ndvi$Elev)
ndvi_normal <- (sectionz_ndvi$NDVI250 - mean(sectionz_ndvi$NDVI250)) / sd(sectionz_ndvi$NDVI250)
fence_normal <- (sectionz_ndvi$dist_fence_clip - mean(sectionz_ndvi$dist_fence_clip)) / sd(sectionz_ndvi$dist_fence_clip)
elec_normal <- (sectionz_ndvi$dist_elec_clip - mean(sectionz_ndvi$dist_elec_clip)) / sd(sectionz_ndvi$dist_elec_clip)
track_normal <- (sectionz_ndvi$dist_track_clip - mean(sectionz_ndvi$dist_track_clip)) / sd(sectionz_ndvi$dist_track_clip)
anyrd_normal <- (sectionz_ndvi$any_road_clip - mean(sectionz_ndvi$any_road_clip)) / sd(sectionz_ndvi$any_road_clip)
wat_normal <- (sectionz_ndvi$dist_wat_clip - mean(sectionz_ndvi$dist_wat_clip)) / sd(sectionz_ndvi$dist_wat_clip)
dam_normal <- (sectionz_ndvi$dist_dam_clip - mean(sectionz_ndvi$dist_dam_clip)) / sd(sectionz_ndvi$dist_dam_clip)
river_normal <- (sectionz_ndvi$dist_river_clip - mean(sectionz_ndvi$dist_river_clip)) / sd(sectionz_ndvi$dist_river_clip)
anywater_normal <- (sectionz_ndvi$any_water_clip - mean(sectionz_ndvi$any_water_clip)) / sd(sectionz_ndvi$any_water_clip)


dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_elev <- (elev_normal - median(elev_normal))^2
dist_med_ndvi <- (ndvi_normal - median(ndvi_normal))^2
dist_med_fence <- (fence_normal - median(fence_normal))^2
dist_med_elec <- (elec_normal - median(elec_normal))^2
dist_med_track <- (track_normal - median(track_normal))^2
dist_med_anyrd <- (anyrd_normal - median(anyrd_normal))^2
dist_med_wat <- (wat_normal - median(wat_normal))^2
dist_med_dam <- (dam_normal - median(dam_normal))^2
dist_med_river <- (river_normal - median(river_normal))^2
dist_med_anywater <- (anywater_normal - median(anywater_normal))^2



dist_med_bok_bc <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_ndvi + dist_med_fence + 
  dist_med_elec + dist_med_track + dist_med_anyrd + dist_med_wat + dist_med_dam + dist_med_river + 
  dist_med_anywater
index_med_bok_bc <- which.min(dist_med_bok_bc)
median_vals_bok_bc <- sectionz_ndvi[index_med_bok_bc,]
distMats_pred_bok_bc <- makeDists(
  cbind(rep(median_vals_bok_bc$x.pos, 200), rep(median_vals_bok_bc$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_bok_bc <- as.data.frame(cbind(pred_bok_bc, sectionz_ndvi$Elev, sectionz_ndvi$NDVI250, 
                                   sectionz_ndvi$dist_fence_clip, sectionz_ndvi$dist_elec_clip, 
                                   sectionz_ndvi$dist_track_clip, sectionz_ndvi$any_road_clip,
                                   sectionz_ndvi$dist_wat_clip, sectionz_ndvi$dist_dam_clip,
                                   sectionz_ndvi$dist_river_clip, sectionz_ndvi$any_water_clip))
colnames(rug_vals_bok_bc) <- c("predictions", "Elev", "NDVI250", "dist_fence", "dist_elec", "dist_track",
                               "any_road", "dist_wat", "dist_dam", "dist_river", "any_water")
rug_vals_bok_bc$dist_fence_km <- rug_vals_bok_bc$dist_fence / 1000
rug_vals_bok_bc$dist_elec_km <- rug_vals_bok_bc$dist_elec / 1000
rug_vals_bok_bc$dist_track_km <- rug_vals_bok_bc$dist_track / 1000
rug_vals_bok_bc$any_road_km <- rug_vals_bok_bc$any_road / 1000
rug_vals_bok_bc$dist_wat_km <- rug_vals_bok_bc$dist_wat / 1000
rug_vals_bok_bc$dist_dam_km <- rug_vals_bok_bc$dist_dam / 1000
rug_vals_bok_bc$dist_river_km <- rug_vals_bok_bc$dist_river / 1000
rug_vals_bok_bc$any_water_km <- rug_vals_bok_bc$any_water / 1000

predict_bok_bc <- as.data.frame(cbind(rep(median_vals_bok_bc$x.pos, 200), 
                                          rep(median_vals_bok_bc$y.pos, 200), 
                                          rep(median_vals_bok_bc$Elev, 200),
                                          rep(median_vals_bok_bc$NDVI250, 200),
                                          rep(median_vals_bok_bc$dist_fence_clip, 200),
                                          rep(median_vals_bok_bc$dist_elec_clip, 200),
                                          rep(median_vals_bok_bc$dist_track_clip, 200),
                                          rep(median_vals_bok_bc$any_road_clip, 200),
                                          rep(median_vals_bok_bc$dist_wat_clip, 200),
                                          rep(median_vals_bok_bc$dist_dam_clip, 200),
                                          rep(median_vals_bok_bc$dist_river_clip, 200),
                                          rep(median_vals_bok_bc$any_water_clip, 200)))
colnames(predict_bok_bc) <- c("x.pos", "y.pos", "Elev", "NDVI250", "dist_fence_clip", 
                                  "dist_elec_clip", "dist_track_clip", "any_road_clip",
                                  "dist_wat_clip", "dist_dam_clip", "dist_river_clip",
                                  "any_water_clip")

predict_elev_bok_bc <- predict_bok_bc
predict_elev_bok_bc$Elev <- elev_predvals
pred_elev_bok_bc <- predict.gamMRSea(newdata=predict_elev_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_elev <- predict.gamMRSea(newdata=predict_elev_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_elev <- cbind(pred_bootz_bok_bc_elev, pred_boot)
}
pred_ci_elev_bok_bc <- t(apply(pred_bootz_bok_bc_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_bok_bc <- cbind(predict_elev_bok_bc, pred_elev_bok_bc, pred_ci_elev_bok_bc)
colnames(pred_df_elev_bok_bc) <- c(colnames(predict_elev_bok_bc), "predictions", "LowerCI", "UpperCI")


predict_ndvi_bok_bc <- predict_bok_bc
predict_ndvi_bok_bc$NDVI250 <- ndvi_predvals
pred_ndvi_bok_bc <- predict.gamMRSea(newdata=predict_ndvi_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_ndvi <- predict.gamMRSea(newdata=predict_ndvi_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_ndvi_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_ndvi <- cbind(pred_bootz_bok_bc_ndvi, pred_boot)
}
pred_ci_ndvi_bok_bc <- t(apply(pred_bootz_bok_bc_ndvi, 1, quantile, probs=c(0.025, 0.975)))
pred_df_ndvi_bok_bc <- cbind(predict_ndvi_bok_bc, pred_ndvi_bok_bc, pred_ci_ndvi_bok_bc)
colnames(pred_df_ndvi_bok_bc) <- c(colnames(predict_ndvi_bok_bc), "predictions", "LowerCI", "UpperCI")

predict_fence_bok_bc <- predict_bok_bc
predict_fence_bok_bc$dist_fence_clip <- fence_predvals
pred_fence_bok_bc <- predict.gamMRSea(newdata=predict_fence_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_fence <- predict.gamMRSea(newdata=predict_fence_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_fence_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_fence <- cbind(pred_bootz_bok_bc_fence, pred_boot)
}
pred_ci_fence_bok_bc <- t(apply(pred_bootz_bok_bc_fence, 1, quantile, probs=c(0.025, 0.975)))
pred_df_fence_bok_bc <- cbind(predict_fence_bok_bc, pred_fence_bok_bc, pred_ci_fence_bok_bc)
colnames(pred_df_fence_bok_bc) <- c(colnames(predict_fence_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_bok_bc$dist_fence_km <- pred_df_fence_bok_bc$dist_fence_clip / 1000

predict_elec_bok_bc <- predict_bok_bc
predict_elec_bok_bc$dist_elec_clip <- elec_predvals
pred_elec_bok_bc <- predict.gamMRSea(newdata=predict_elec_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_elec <- predict.gamMRSea(newdata=predict_elec_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elec_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_elec <- cbind(pred_bootz_bok_bc_elec, pred_boot)
}
pred_ci_elec_bok_bc <- t(apply(pred_bootz_bok_bc_elec, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elec_bok_bc <- cbind(predict_elec_bok_bc, pred_elec_bok_bc, pred_ci_elec_bok_bc)
colnames(pred_df_elec_bok_bc) <- c(colnames(predict_elec_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_elec_bok_bc$dist_elec_km <- pred_df_elec_bok_bc$dist_elec_clip / 1000

predict_track_bok_bc <- predict_bok_bc
predict_track_bok_bc$dist_track_clip <- track_predvals
pred_track_bok_bc <- predict.gamMRSea(newdata=predict_track_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_track <- predict.gamMRSea(newdata=predict_track_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_track <- cbind(pred_bootz_bok_bc_track, pred_boot)
}
pred_ci_track_bok_bc <- t(apply(pred_bootz_bok_bc_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_bok_bc <- cbind(predict_track_bok_bc, pred_track_bok_bc, pred_ci_track_bok_bc)
colnames(pred_df_track_bok_bc) <- c(colnames(predict_track_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_track_bok_bc$dist_track_km <- pred_df_track_bok_bc$dist_track_clip / 1000

predict_anyrd_bok_bc <- predict_bok_bc
predict_anyrd_bok_bc$any_road_clip <- anyrd_predvals
pred_anyrd_bok_bc <- predict.gamMRSea(newdata=predict_anyrd_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_anyrd <- predict.gamMRSea(newdata=predict_anyrd_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_anyrd <- cbind(pred_bootz_bok_bc_anyrd, pred_boot)
}
pred_ci_anyrd_bok_bc <- t(apply(pred_bootz_bok_bc_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_bok_bc <- cbind(predict_anyrd_bok_bc, pred_anyrd_bok_bc, pred_ci_anyrd_bok_bc)
colnames(pred_df_anyrd_bok_bc) <- c(colnames(predict_anyrd_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anyrd_bok_bc$any_road_km <- pred_df_anyrd_bok_bc$any_road_clip / 1000

predict_wat_bok_bc <- predict_bok_bc
predict_wat_bok_bc$dist_wat_clip <- wat_predvals
pred_wat_bok_bc <- predict.gamMRSea(newdata=predict_wat_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_wat <- predict.gamMRSea(newdata=predict_wat_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_wat_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_wat <- cbind(pred_bootz_bok_bc_wat, pred_boot)
}
pred_ci_wat_bok_bc <- t(apply(pred_bootz_bok_bc_wat, 1, quantile, probs=c(0.025, 0.975)))
pred_df_wat_bok_bc <- cbind(predict_wat_bok_bc, pred_wat_bok_bc, pred_ci_wat_bok_bc)
colnames(pred_df_wat_bok_bc) <- c(colnames(predict_wat_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_wat_bok_bc$dist_wat_km <- pred_df_wat_bok_bc$dist_wat_clip / 1000

predict_dam_bok_bc <- predict_bok_bc
predict_dam_bok_bc$dist_dam_clip <- dam_predvals
pred_dam_bok_bc <- predict.gamMRSea(newdata=predict_dam_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_dam <- predict.gamMRSea(newdata=predict_dam_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_dam_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_dam <- cbind(pred_bootz_bok_bc_dam, pred_boot)
}
pred_ci_dam_bok_bc <- t(apply(pred_bootz_bok_bc_dam, 1, quantile, probs=c(0.025, 0.975)))
pred_df_dam_bok_bc <- cbind(predict_dam_bok_bc, pred_dam_bok_bc, pred_ci_dam_bok_bc)
colnames(pred_df_dam_bok_bc) <- c(colnames(predict_dam_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_dam_bok_bc$dist_dam_km <- pred_df_dam_bok_bc$dist_dam_clip / 1000

predict_river_bok_bc <- predict_bok_bc
predict_river_bok_bc$dist_river_clip <- river_predvals
pred_river_bok_bc <- predict.gamMRSea(newdata=predict_river_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_river <- predict.gamMRSea(newdata=predict_river_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_river_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_river <- cbind(pred_bootz_bok_bc_river, pred_boot)
}
pred_ci_river_bok_bc <- t(apply(pred_bootz_bok_bc_river, 1, quantile, probs=c(0.025, 0.975)))
pred_df_river_bok_bc <- cbind(predict_river_bok_bc, pred_river_bok_bc, pred_ci_river_bok_bc)
colnames(pred_df_river_bok_bc) <- c(colnames(predict_river_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_bok_bc$dist_river_km <- pred_df_river_bok_bc$dist_river_clip / 1000

predict_anywt_bok_bc <- predict_bok_bc
predict_anywt_bok_bc$any_water_clip <- anywater_predvals
pred_anywt_bok_bc <- predict.gamMRSea(newdata=predict_anywt_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist)
pred_bootz_bok_bc_anywt <- predict.gamMRSea(newdata=predict_anywt_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_bok_bc, object=bok_mod_bc$mod_2d$bestModel, g2k=distMats_pred_bok_bc$dataDist, coeff=rcoefs_bok_bc[bt,])
  pred_bootz_bok_bc_anywt <- cbind(pred_bootz_bok_bc_anywt, pred_boot)
}
pred_ci_anywt_bok_bc <- t(apply(pred_bootz_bok_bc_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_bok_bc <- cbind(predict_anywt_bok_bc, pred_anywt_bok_bc, pred_ci_anywt_bok_bc)
colnames(pred_df_anywt_bok_bc) <- c(colnames(predict_anywt_bok_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anywt_bok_bc$any_water_km <- pred_df_anywt_bok_bc$any_water_clip / 1000



pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_bok_bc] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")

df_pointz <- as.data.frame(cbind(df_pointz, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
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
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(6, 6)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))


plotout_elev_bok_bc <- ggplot(pred_df_elev_bok_bc, aes(x=Elev, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Elevation (m)') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_ndvi_bok_bc <- ggplot(pred_df_ndvi_bok_bc, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_fence_bok_bc <- ggplot(pred_df_fence_bok_bc, aes(x=dist_fence_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to fence (km)') + geom_ribbon(aes(x=dist_fence_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_elec_bok_bc <- ggplot(pred_df_elec_bok_bc, aes(x=dist_elec_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to electricity line (km)') + geom_ribbon(aes(x=dist_elec_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_track_bok_bc <- ggplot(pred_df_track_bok_bc, aes(x=dist_track_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to track (km)') + geom_ribbon(aes(x=dist_track_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anyrd_bok_bc <- ggplot(pred_df_anyrd_bok_bc, aes(x=any_road_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any road (km)') + geom_ribbon(aes(x=any_road_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_wat_bok_bc <- ggplot(pred_df_wat_bok_bc, aes(x=dist_wat_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to waterhole (km)') + geom_ribbon(aes(x=dist_wat_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_dam_bok_bc <- ggplot(pred_df_dam_bok_bc, aes(x=dist_dam_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to dam (km)') + geom_ribbon(aes(x=dist_dam_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_river_bok_bc <- ggplot(pred_df_river_bok_bc, aes(x=dist_river_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to river (km)') + geom_ribbon(aes(x=dist_river_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anywt_bok_bc <- ggplot(pred_df_anywt_bok_bc, aes(x=any_water_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any water (km)') + geom_ribbon(aes(x=any_water_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_bok_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed_bok.png"))
pt_med
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Elev_bok.png"))
plotout_elev_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2NDVI_bok.png"))
plotout_ndvi_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Fence_bok.png"))
plotout_fence_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Elec_bok.png"))
plotout_elec_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Track_bok.png"))
plotout_track_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Anyrd_bok.png"))
plotout_anyrd_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2wat_bok.png"))
plotout_wat_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2dam_bok.png"))
plotout_dam_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2river_bok.png"))
plotout_river_bok_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2anywater_bok.png"))
plotout_anywt_bok_bc
dev.off()























