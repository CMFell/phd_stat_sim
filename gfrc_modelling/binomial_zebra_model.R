# This file fits the binomial zebra model and creates output images
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

sectionz_ndvi$response <- sectionz_ndvi$ZebraPA
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
zebra_2d_mod <- runSALSA2D(
  PA_mod, 
  salsa2dlist_in, 
  d2k=distMats$dataDist, 
  k2k=distMats$knotDist,
  suppress.printout = TRUE
)
an_res_zebra <- Anova(zebra_2d_mod$bestModel, test='F')
rownames(an_res_zebra) <- c(zebra_2d_mod$bestModel$varshortnames, "2dsmooth", "Residuals")

coefz_zebra_sm <- zebra_2d_mod$bestModel$coefficients
covmat_zebra_sm <- summary(zebra_2d_mod$bestModel)$cov.robust
rcoefs_zebra_sm <- rmvnorm(1000, coefz_zebra_sm, sigma=covmat_zebra_sm)

pred_bootz_zebra_sm_all <- predict.gamMRSea(object=zebra_2d_mod$bestModel, coeff=rcoefs_zebra_sm[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=zebra_2d_mod$bestModel, coeff=rcoefs_zebra_sm[bt,])
  pred_bootz_zebra_sm_all <- cbind(pred_bootz_zebra_sm_all, pred_boot)
}
pred_ci_zebra_sm <- t(apply(pred_bootz_zebra_sm_all, 1, quantile, probs=c(0.025, 0.975)))

pred_zebra_sm <- predict.gamMRSea(object=zebra_2d_mod$bestModel)
data_zebra_sm <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_sm, pred_ci_zebra_sm))
colnames(data_zebra_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_zebra_sm <- cbind(data_zebra_sm, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_zebra_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")

zebra_pred_sm <- ggplot(data_zebra_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
zebra_pred_sm <- ggplot(data_zebra_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_sm_lw <- ggplot(data_zebra_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette ="YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
zebra_pred_sm_lw <- ggplot(data_zebra_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_sm_hi <- ggplot(data_zebra_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
zebra_pred_sm_hi <- ggplot(data_zebra_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "ZebraModel1Pred.png"))
zebra_pred_sm
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel1Lw.png"))
zebra_pred_sm_lw
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel1Hi.png"))
zebra_pred_sm_hi
dev.off()


### Not included in thesis


# Second model base selection
var_list_base <- c("Elev", "dist_track_clip", "any_road_clip", "dist_minor_clip", "any_water_clip")
zebra_mod_base <- rez_for_var_list(var_list_base, PA_mod, distMats)

coefz_zebra_base <- zebra_mod_base$mod_2d$bestModel$coefficients
covmat_zebra_base <- summary(zebra_mod_base$mod_2d$bestModel)$cov.robust
rcoefs_zebra_base <- rmvnorm(1000, coefz_zebra_base, sigma=covmat_zebra_base)

pred_bootz_zebra_base_all <- predict.gamMRSea(object=zebra_mod_base$mod_2d$bestModel, coeff=rcoefs_zebra_base[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=zebra_mod_base$mod_2d$bestModel, coeff=rcoefs_zebra_base[bt,])
  pred_bootz_zebra_base_all <- cbind(pred_bootz_zebra_base_all, pred_boot)
}
pred_ci_zebra_base <- t(apply(pred_bootz_zebra_base_all, 1, quantile, probs=c(0.025, 0.975)))


pred_zebra_base <- predict.gamMRSea(object=zebra_mod_base$mod_2d$bestModel)
data_zebra_base <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_base, pred_ci_zebra_base))
colnames(data_zebra_base) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

graph_pred_zebra_base <- ggplot(data_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette =  "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_zebra_base_lw <- ggplot(data_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_zebra_base_hi <- ggplot(data_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "ZebraModel2Pred.png"))
graph_pred_zebra_base
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel2Lw.png"))
graph_pred_zebra_base_lw
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel2Hi.png"))
graph_pred_zebra_base_hi
dev.off()


pred_bootz_zebra_base_all_link <- predict.gamMRSea(object=zebra_mod_base$mod_2d$bestModel, coeff=rcoefs_zebra_base[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=zebra_mod_base$mod_2d$bestModel, coeff=rcoefs_zebra_base[bt,], type='link')
  pred_bootz_zebra_base_all_link <- cbind(pred_bootz_zebra_base_all, pred_boot)
}
pred_ci_zebra_base_link <- t(apply(pred_bootz_zebra_base_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_zebra_base_link <- predict.gamMRSea(object=zebra_mod_base$mod_2d$bestModel)
data_zebra_base_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_base_link, pred_ci_zebra_base_link))
colnames(data_zebra_base_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


# filter out any NAs
testcoefs_zebra_base <- zebra_mod_base$mod_2d$bestModel$coefficients
cfmsk_zebra_base <- !is.na(testcoefs_zebra_base)
tstcfs_zebra_base <- testcoefs_zebra_base[cfmsk_zebra_base]

sp_col_zebra_base <- zebra_mod_base$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_zebra_base <- str_sub(names(tstcfs_zebra_base),1,strlenclnam)
coefmask_zebra_base <- coefnamsplt_zebra_base == nam2d

# create radial gaussian bases
radii_zebra_base <- sp_col_zebra_base$radii
radiiInd_zebra_base <- sp_col_zebra_base$radiusIndices
aR_zebra_base <- sp_col_zebra_base$knotPos
lrf_zebra_base <- LRF.g(radiiInd_zebra_base, distMats$dataDist, radii_zebra_base, aR_zebra_base)

# combine coefmask and facts
coefz_zebra_base <- tstcfs_zebra_base[coefmask_zebra_base]
# get predicted values on link scale
predtm_zebra_base <- lrf_zebra_base %*% coefz_zebra_base
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_zebra_base <- rcoefs_zebra_base[1,][cfmsk_zebra_base][coefmask_zebra_base]
predboot_zebra_base_2d <- lrf_zebra_base %*% bootcoefz_zebra_base
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_zebra_base[bt,][cfmsk_zebra_base][coefmask_zebra_base]
  predbt <- lrf_zebra_base %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_zebra_base_2d <- cbind(predboot_zebra_base_2d, predbt)
}
pred_2d_ci_zebra_base <- t(apply(predboot_zebra_base_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_zebra_base <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_zebra_base, pred_2d_ci_zebra_base))
colnames(data_2d_zebra_base) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

#
zebra_2d_base <- ggplot(data_2d_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_2d_base_lw <- ggplot(data_2d_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_2d_base_hi <- ggplot(data_2d_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_zebra_base <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_base_link - predtm_zebra_base, pred_ci_zebra_base_link - pred_2d_ci_zebra_base))
colnames(data_not2d_zebra_base) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


zebra_not2d_base <- ggplot(data_not2d_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_not2d_base_lw <- ggplot(data_not2d_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_not2d_base_hi <- ggplot(data_not2d_zebra_base, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))


zebra_pred_base_link <- ggplot(data_zebra_base_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_base_lw_link <- ggplot(data_zebra_base_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_base_hi_link <- ggplot(data_zebra_base_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "ZebraModel22d.png"))
  zebra_2d_base
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel22dLw.png"))
  zebra_2d_base_lw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel22dHi.png"))
  zebra_2d_base_hi
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel2NotPrd.png"))
  zebra_not2d_base
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel2NotLw.png"))
  zebra_not2d_base_lw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel2NotHi.png"))
  zebra_not2d_base_hi
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel2Lnk.png"))
zebra_pred_base_link
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel2LnkLw.png"))
zebra_pred_base_lw_link
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel2LnkHi.png"))
zebra_pred_base_hi_link
dev.off()

## 1d calculations
# Elev, dist track, any road, dist minor, any water dam road
# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
track_predvals <- seq(min(sectionz_ndvi$dist_track_clip), max(sectionz_ndvi$dist_track_clip), length.out=200)
anyrd_predvals <- seq(min(sectionz_ndvi$any_road_clip), max(sectionz_ndvi$any_road_clip), length.out=200)
minor_predvals <- seq(min(sectionz_ndvi$dist_minor_clip), max(sectionz_ndvi$dist_minor_clip), length.out=200)
dam_predvals <- seq(min(sectionz_ndvi$dist_dam_clip), max(sectionz_ndvi$dist_dam_clip), length.out=200)
anywater_predvals <- seq(min(sectionz_ndvi$any_water_clip), max(sectionz_ndvi$any_water_clip), length.out=200)
road_predvals <- seq(min(sectionz_ndvi$dist_road_clip), max(sectionz_ndvi$dist_road_clip), length.out=200)
#ndvi, elev, fence, track, any road, wat, dam, any water

# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
elev_normal <- (sectionz_ndvi$Elev - mean(sectionz_ndvi$Elev)) / sd(sectionz_ndvi$Elev)
track_normal <- (sectionz_ndvi$dist_track_clip - mean(sectionz_ndvi$dist_track_clip)) / sd(sectionz_ndvi$dist_track_clip)
anyrd_normal <- (sectionz_ndvi$any_road_clip - mean(sectionz_ndvi$any_road_clip)) / sd(sectionz_ndvi$any_road_clip)
minor_normal <- (sectionz_ndvi$dist_minor_clip - mean(sectionz_ndvi$dist_minor_clip)) / sd(sectionz_ndvi$dist_minor_clip)
dam_normal <- (sectionz_ndvi$dist_dam_clip - mean(sectionz_ndvi$dist_dam_clip)) / sd(sectionz_ndvi$dist_dam_clip)
anywater_normal <- (sectionz_ndvi$any_water_clip - mean(sectionz_ndvi$any_water_clip)) / sd(sectionz_ndvi$any_water_clip)
road_normal <- (sectionz_ndvi$dist_road_clip - mean(sectionz_ndvi$dist_road_clip)) / sd(sectionz_ndvi$dist_road_clip)

dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_elev <- (elev_normal - median(elev_normal))^2
dist_med_track <- (track_normal - median(track_normal))^2
dist_med_anyrd <- (anyrd_normal - median(anyrd_normal))^2
dist_med_minor <- (minor_normal - median(minor_normal))^2
dist_med_dam <- (dam_normal - median(dam_normal))^2
dist_med_anywater <- (anywater_normal - median(anywater_normal))^2
dist_med_road <- (road_normal - median(road_normal))^2


dist_med_zebra_base <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_track + dist_med_anyrd + 
  dist_med_minor + dist_med_anywater
index_med_zebra_base <- which.min(dist_med_zebra_base)
median_vals_zebra_base <- sectionz_ndvi[index_med_zebra_base,]
distMats_pred_zebra_base <- makeDists(
  cbind(rep(median_vals_zebra_base$x.pos, 200), rep(median_vals_zebra_base$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_zebra_base <- as.data.frame(cbind(pred_zebra_base, sectionz_ndvi$Elev, sectionz_ndvi$dist_track_clip, 
                                   sectionz_ndvi$any_road_clip,sectionz_ndvi$dist_minor_clip, 
                                   sectionz_ndvi$any_water_clip))
colnames(rug_vals_zebra_base) <- c("predictions", "Elev", "dist_track", "any_road", "dist_minor", "any_water")

predict_zebra_base <- as.data.frame(cbind(rep(median_vals_zebra_base$x.pos, 200), 
                                          rep(median_vals_zebra_base$y.pos, 200), 
                                          rep(median_vals_zebra_base$Elev, 200),
                                          rep(median_vals_zebra_base$dist_track_clip, 200),
                                          rep(median_vals_zebra_base$any_road_clip, 200),
                                          rep(median_vals_zebra_base$dist_minor_clip, 200),
                                          rep(median_vals_zebra_base$any_water_clip, 200)))
colnames(predict_zebra_base) <- c("x.pos", "y.pos", "Elev", "dist_track_clip", "any_road_clip", "dist_minor_clip", "any_water_clip")

predict_elev_zebra_base <- predict_zebra_base
predict_elev_zebra_base$Elev <- elev_predvals
pred_elev_zebra_base <- predict.gamMRSea(newdata=predict_elev_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist)
pred_bootz_zebra_base_elev <- predict.gamMRSea(newdata=predict_elev_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[bt,])
  pred_bootz_zebra_base_elev <- cbind(pred_bootz_zebra_base_elev, pred_boot)
}
pred_ci_elev_zebra_base <- t(apply(pred_bootz_zebra_base_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_zebra_base <- cbind(predict_elev_zebra_base, pred_elev_zebra_base, pred_ci_elev_zebra_base)
colnames(pred_df_elev_zebra_base) <- c(colnames(predict_elev_zebra_base), "predictions", "LowerCI", "UpperCI")

predict_track_zebra_base <- predict_zebra_base
predict_track_zebra_base$dist_track_clip <- track_predvals
pred_track_zebra_base <- predict.gamMRSea(newdata=predict_track_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist)
pred_bootz_zebra_base_track <- predict.gamMRSea(newdata=predict_track_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[bt,])
  pred_bootz_zebra_base_track <- cbind(pred_bootz_zebra_base_track, pred_boot)
}
pred_ci_track_zebra_base <- t(apply(pred_bootz_zebra_base_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_zebra_base <- cbind(predict_track_zebra_base, pred_track_zebra_base, pred_ci_track_zebra_base)
colnames(pred_df_track_zebra_base) <- c(colnames(predict_track_zebra_base), "predictions", "LowerCI", "UpperCI")

predict_anyrd_zebra_base <- predict_zebra_base
predict_anyrd_zebra_base$any_road_clip <- anyrd_predvals
pred_anyrd_zebra_base <- predict.gamMRSea(newdata=predict_anyrd_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist)
pred_bootz_zebra_base_anyrd <- predict.gamMRSea(newdata=predict_anyrd_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[bt,])
  pred_bootz_zebra_base_anyrd <- cbind(pred_bootz_zebra_base_anyrd, pred_boot)
}
pred_ci_anyrd_zebra_base <- t(apply(pred_bootz_zebra_base_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_zebra_base <- cbind(predict_anyrd_zebra_base, pred_anyrd_zebra_base, pred_ci_anyrd_zebra_base)
colnames(pred_df_anyrd_zebra_base) <- c(colnames(predict_anyrd_zebra_base), "predictions", "LowerCI", "UpperCI")

predict_minor_zebra_base <- predict_zebra_base
predict_minor_zebra_base$dist_minor_clip <- minor_predvals
pred_minor_zebra_base <- predict.gamMRSea(newdata=predict_minor_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist)
pred_bootz_zebra_base_minor <- predict.gamMRSea(newdata=predict_minor_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_minor_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[bt,])
  pred_bootz_zebra_base_minor <- cbind(pred_bootz_zebra_base_minor, pred_boot)
}
pred_ci_minor_zebra_base <- t(apply(pred_bootz_zebra_base_minor, 1, quantile, probs=c(0.025, 0.975)))
pred_df_minor_zebra_base <- cbind(predict_minor_zebra_base, pred_minor_zebra_base, pred_ci_minor_zebra_base)
colnames(pred_df_minor_zebra_base) <- c(colnames(predict_minor_zebra_base), "predictions", "LowerCI", "UpperCI")

predict_anywt_zebra_base <- predict_zebra_base
predict_anywt_zebra_base$any_water_clip <- anywater_predvals
pred_anywt_zebra_base <- predict.gamMRSea(newdata=predict_anywt_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist)
pred_bootz_zebra_base_anywt <- predict.gamMRSea(newdata=predict_anywt_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_zebra_base, object=zebra_mod_base$mod_2d$bestModel, g2k=distMats_pred_zebra_base$dataDist, coeff=rcoefs_zebra_base[bt,])
  pred_bootz_zebra_base_anywt <- cbind(pred_bootz_zebra_base_anywt, pred_boot)
}
pred_ci_anywt_zebra_base <- t(apply(pred_bootz_zebra_base_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_zebra_base <- cbind(predict_anywt_zebra_base, pred_anywt_zebra_base, pred_ci_anywt_zebra_base)
colnames(pred_df_anywt_zebra_base) <- c(colnames(predict_anywt_zebra_base), "predictions", "LowerCI", "UpperCI")


pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_zebra_base] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")

pt_med2 <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of median point used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.05, 0.05)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))

plotout_elev_zebra_base <- ggplot(pred_df_elev_zebra_base, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_base, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_track_zebra_base <- ggplot(pred_df_track_zebra_base, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_base, sides='b', size=0.5, length=unit(0.02,'npc'))+ theme(text = element_text(size=20))
plotout_anyrd_zebra_base <- ggplot(pred_df_anyrd_zebra_base, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to any road') + geom_ribbon(aes(x=any_road_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_base, sides='b', size=0.5, length=unit(0.02,'npc'))+ theme(text = element_text(size=20))
plotout_minor_zebra_base <- ggplot(pred_df_minor_zebra_base, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor river') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_base, sides='b', size=0.5, length=unit(0.02,'npc'))+ theme(text = element_text(size=20))
plotout_anywt_zebra_base <- ggplot(pred_df_anywt_zebra_base, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to any water') + geom_ribbon(aes(x=any_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_base, sides='b', size=0.5, length=unit(0.02,'npc'))+ theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed2_zebra_base.png"))
pt_med2
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Elev_zebra.png"))
plotout_elev_zebra_base
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Track_zebra.png"))
plotout_track_zebra_base
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Anyrd_zebra.png"))
plotout_anyrd_zebra_base
dev.off()
png(filename=paste0(img_save_dir, "PtMod2minor_zebra.png"))
plotout_minor_zebra_base
dev.off()
png(filename=paste0(img_save_dir, "PtMod2anywater_zebra.png"))
plotout_anywt_zebra_base
dev.off()


### Included in thesis below



# Third model backward selection
var_list_bc <- c("Elev", "dist_track_clip", "any_road_clip", "dist_minor_clip", "any_water_clip", "dist_dam_clip")
zebra_mod_bc <- rez_for_var_list(var_list_bc, PA_mod, distMats)


coefz_zebra_bc <- zebra_mod_bc$mod_2d$bestModel$coefficients
covmat_zebra_bc <- summary(zebra_mod_bc$mod_2d$bestModel)$cov.robust
rcoefs_zebra_bc <- rmvnorm(1000, coefz_zebra_bc, sigma=covmat_zebra_bc)

pred_bootz_zebra_bc_all <- predict.gamMRSea(object=zebra_mod_bc$mod_2d$bestModel, coeff=rcoefs_zebra_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=zebra_mod_bc$mod_2d$bestModel, coeff=rcoefs_zebra_bc[bt,])
  pred_bootz_zebra_bc_all <- cbind(pred_bootz_zebra_bc_all, pred_boot)
}
pred_ci_zebra_bc <- t(apply(pred_bootz_zebra_bc_all, 1, quantile, probs=c(0.025, 0.975)))

pred_zebra_bc <- predict.gamMRSea(object=zebra_mod_bc$mod_2d$bestModel)
data_zebra_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_bc, pred_ci_zebra_bc))
colnames(data_zebra_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_zebra_bc <- cbind(data_zebra_bc, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_zebra_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY","distX","distY")

graph_pred_zebra_bc <- ggplot(data_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_zebra_bc <- ggplot(data_zebra_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_zebra_bc_lw <- ggplot(data_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_zebra_bc_lw <- ggplot(data_zebra_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_zebra_bc_hi <- ggplot(data_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_zebra_bc_hi <- ggplot(data_zebra_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "ZebraModel3Pred.png"))
graph_pred_zebra_bc
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel3Lw.png"))
graph_pred_zebra_bc_lw
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel3Hi.png"))
graph_pred_zebra_bc_hi
dev.off()


### Not included in thesis

pred_bootz_zebra_bc_all_link <- predict.gamMRSea(object=zebra_mod_bc$mod_2d$bestModel, coeff=rcoefs_zebra_bc[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=zebra_mod_bc$mod_2d$bestModel, coeff=rcoefs_zebra_bc[bt,], type='link')
  pred_bootz_zebra_bc_all_link <- cbind(pred_bootz_zebra_bc_all, pred_boot)
}
pred_ci_zebra_bc_link <- t(apply(pred_bootz_zebra_bc_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_zebra_bc_link <- predict.gamMRSea(object=zebra_mod_bc$mod_2d$bestModel)
data_zebra_bc_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_bc_link, pred_ci_zebra_bc_link))
colnames(data_zebra_bc_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


# filter out any NAs
testcoefs_zebra_bc <- zebra_mod_bc$mod_2d$bestModel$coefficients
cfmsk_zebra_bc <- !is.na(testcoefs_zebra_bc)
tstcfs_zebra_bc <- testcoefs_zebra_bc[cfmsk_zebra_bc]

sp_col_zebra_bc <- zebra_mod_bc$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_zebra_bc <- str_sub(names(tstcfs_zebra_bc),1,strlenclnam)
coefmask_zebra_bc <- coefnamsplt_zebra_bc == nam2d

# create radial gaussian bases
radii_zebra_bc <- sp_col_zebra_bc$radii
radiiInd_zebra_bc <- sp_col_zebra_bc$radiusIndices
aR_zebra_bc <- sp_col_zebra_bc$knotPos
lrf_zebra_bc <- LRF.g(radiiInd_zebra_bc, distMats$dataDist, radii_zebra_bc, aR_zebra_bc)

# combine coefmask and facts
coefz_zebra_bc <- tstcfs_zebra_bc[coefmask_zebra_bc]
# get predicted values on link scale
predtm_zebra_bc <- lrf_zebra_bc %*% coefz_zebra_bc
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_zebra_bc <- rcoefs_zebra_bc[1,][cfmsk_zebra_bc][coefmask_zebra_bc]
predboot_zebra_bc_2d <- lrf_zebra_bc %*% bootcoefz_zebra_bc
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_zebra_bc[bt,][cfmsk_zebra_bc][coefmask_zebra_bc]
  predbt <- lrf_zebra_bc %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_zebra_bc_2d <- cbind(predboot_zebra_bc_2d, predbt)
}
pred_2d_ci_zebra_bc <- t(apply(predboot_zebra_bc_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_zebra_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_zebra_bc, pred_2d_ci_zebra_bc))
colnames(data_2d_zebra_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


zebra_2d_bc <- ggplot(data_2d_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_2d_bc_lw <- ggplot(data_2d_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_2d_bc_hi <- ggplot(data_2d_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_zebra_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_bc_link - predtm_zebra_bc, pred_ci_zebra_bc_link - pred_2d_ci_zebra_bc))
colnames(data_not2d_zebra_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


zebra_not2d_bc <- ggplot(data_not2d_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_not2d_bc_lw <- ggplot(data_not2d_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_not2d_bc_hi <- ggplot(data_not2d_zebra_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))


zebra_pred_bc_link <- ggplot(data_zebra_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_bc_lw_link <- ggplot(data_zebra_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_bc_hi_link <- ggplot(data_zebra_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,4), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "ZebraModel32d.png"))
zebra_2d_bc
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel32dLw.png"))
zebra_2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel32dHi.png"))
zebra_2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel3NotPrd.png"))
zebra_not2d_bc
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel3NotLw.png"))
zebra_not2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel3NotHi.png"))
zebra_not2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel3Lnk.png"))
zebra_pred_bc_link
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel3LnkLw.png"))
zebra_pred_bc_lw_link
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel3LnkHi.png"))
zebra_pred_bc_hi_link
dev.off()

### Included in thesis below

## 1d calculations
# Elev, dist track, any road, dist minor, any water dam road
# create prediction vals for variables

dist_med_zebra_bc <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_track + dist_med_anyrd + 
  dist_med_minor + dist_med_anywater + dist_med_dam
index_med_zebra_bc <- which.min(dist_med_zebra_bc)
median_vals_zebra_bc <- sectionz_ndvi[index_med_zebra_bc,]
distMats_pred_zebra_bc <- makeDists(
  cbind(rep(median_vals_zebra_bc$x.pos, 200), rep(median_vals_zebra_bc$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_zebra_bc <- as.data.frame(cbind(pred_zebra_bc, sectionz_ndvi$Elev, sectionz_ndvi$dist_track_clip, 
                                           sectionz_ndvi$any_road_clip,sectionz_ndvi$dist_minor_clip, 
                                           sectionz_ndvi$any_water_clip, sectionz_ndvi$dist_dam_clip))
colnames(rug_vals_zebra_bc) <- c("predictions", "Elev", "dist_track", "any_road", "dist_minor", "any_water", "dist_dam")
rug_vals_zebra_bc$dist_track_km <- rug_vals_zebra_bc$dist_track / 1000
rug_vals_zebra_bc$any_road_km <- rug_vals_zebra_bc$any_road / 1000
rug_vals_zebra_bc$dist_minor_km <- rug_vals_zebra_bc$dist_minor / 1000
rug_vals_zebra_bc$any_water_km <- rug_vals_zebra_bc$any_water / 1000
rug_vals_zebra_bc$dist_dam_km <- rug_vals_zebra_bc$dist_dam / 1000

predict_zebra_bc <- as.data.frame(cbind(rep(median_vals_zebra_bc$x.pos, 200), 
                                          rep(median_vals_zebra_bc$y.pos, 200), 
                                          rep(median_vals_zebra_bc$Elev, 200),
                                          rep(median_vals_zebra_bc$dist_track_clip, 200),
                                          rep(median_vals_zebra_bc$any_road_clip, 200),
                                          rep(median_vals_zebra_bc$dist_minor_clip, 200),
                                          rep(median_vals_zebra_bc$any_water_clip, 200),
                                          rep(median_vals_zebra_bc$dist_dam_clip, 200)))
colnames(predict_zebra_bc) <- c("x.pos", "y.pos", "Elev", "dist_track_clip", "any_road_clip", "dist_minor_clip", 
                                "any_water_clip", "dist_dam_clip")

predict_elev_zebra_bc <- predict_zebra_bc
predict_elev_zebra_bc$Elev <- elev_predvals
pred_elev_zebra_bc <- predict.gamMRSea(newdata=predict_elev_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist)
pred_bootz_zebra_bc_elev <- predict.gamMRSea(newdata=predict_elev_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[bt,])
  pred_bootz_zebra_bc_elev <- cbind(pred_bootz_zebra_bc_elev, pred_boot)
}
pred_ci_elev_zebra_bc <- t(apply(pred_bootz_zebra_bc_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_zebra_bc <- cbind(predict_elev_zebra_bc, pred_elev_zebra_bc, pred_ci_elev_zebra_bc)
colnames(pred_df_elev_zebra_bc) <- c(colnames(predict_elev_zebra_bc), "predictions", "LowerCI", "UpperCI")

predict_track_zebra_bc <- predict_zebra_bc
predict_track_zebra_bc$dist_track_clip <- track_predvals
pred_track_zebra_bc <- predict.gamMRSea(newdata=predict_track_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist)
pred_bootz_zebra_bc_track <- predict.gamMRSea(newdata=predict_track_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[bt,])
  pred_bootz_zebra_bc_track <- cbind(pred_bootz_zebra_bc_track, pred_boot)
}
pred_ci_track_zebra_bc <- t(apply(pred_bootz_zebra_bc_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_zebra_bc <- cbind(predict_track_zebra_bc, pred_track_zebra_bc, pred_ci_track_zebra_bc)
colnames(pred_df_track_zebra_bc) <- c(colnames(predict_track_zebra_bc), "predictions", "LowerCI", "UpperCI")
pred_df_track_zebra_bc$dist_track_km <- pred_df_track_zebra_bc$dist_track_clip / 1000

predict_anyrd_zebra_bc <- predict_zebra_bc
predict_anyrd_zebra_bc$any_road_clip <- anyrd_predvals
pred_anyrd_zebra_bc <- predict.gamMRSea(newdata=predict_anyrd_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist)
pred_bootz_zebra_bc_anyrd <- predict.gamMRSea(newdata=predict_anyrd_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[bt,])
  pred_bootz_zebra_bc_anyrd <- cbind(pred_bootz_zebra_bc_anyrd, pred_boot)
}
pred_ci_anyrd_zebra_bc <- t(apply(pred_bootz_zebra_bc_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_zebra_bc <- cbind(predict_anyrd_zebra_bc, pred_anyrd_zebra_bc, pred_ci_anyrd_zebra_bc)
colnames(pred_df_anyrd_zebra_bc) <- c(colnames(predict_anyrd_zebra_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anyrd_zebra_bc$any_road_km <- pred_df_anyrd_zebra_bc$any_road_clip / 1000

predict_minor_zebra_bc <- predict_zebra_bc
predict_minor_zebra_bc$dist_minor_clip <- minor_predvals
pred_minor_zebra_bc <- predict.gamMRSea(newdata=predict_minor_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist)
pred_bootz_zebra_bc_minor <- predict.gamMRSea(newdata=predict_minor_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_minor_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[bt,])
  pred_bootz_zebra_bc_minor <- cbind(pred_bootz_zebra_bc_minor, pred_boot)
}
pred_ci_minor_zebra_bc <- t(apply(pred_bootz_zebra_bc_minor, 1, quantile, probs=c(0.025, 0.975)))
pred_df_minor_zebra_bc <- cbind(predict_minor_zebra_bc, pred_minor_zebra_bc, pred_ci_minor_zebra_bc)
colnames(pred_df_minor_zebra_bc) <- c(colnames(predict_minor_zebra_bc), "predictions", "LowerCI", "UpperCI")
pred_df_minor_zebra_bc$dist_minor_km <- pred_df_minor_zebra_bc$dist_minor_clip / 1000

predict_anywt_zebra_bc <- predict_zebra_bc
predict_anywt_zebra_bc$any_water_clip <- anywater_predvals
pred_anywt_zebra_bc <- predict.gamMRSea(newdata=predict_anywt_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist)
pred_bootz_zebra_bc_anywt <- predict.gamMRSea(newdata=predict_anywt_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[bt,])
  pred_bootz_zebra_bc_anywt <- cbind(pred_bootz_zebra_bc_anywt, pred_boot)
}
pred_ci_anywt_zebra_bc <- t(apply(pred_bootz_zebra_bc_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_zebra_bc <- cbind(predict_anywt_zebra_bc, pred_anywt_zebra_bc, pred_ci_anywt_zebra_bc)
colnames(pred_df_anywt_zebra_bc) <- c(colnames(predict_anywt_zebra_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anywt_zebra_bc$any_water_km <- pred_df_anywt_zebra_bc$any_water_clip / 1000

predict_dam_zebra_bc <- predict_zebra_bc
predict_dam_zebra_bc$dist_dam_clip <- dam_predvals
pred_dam_zebra_bc <- predict.gamMRSea(newdata=predict_dam_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist)
pred_bootz_zebra_bc_dam <- predict.gamMRSea(newdata=predict_dam_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_dam_zebra_bc, object=zebra_mod_bc$mod_2d$bestModel, g2k=distMats_pred_zebra_bc$dataDist, coeff=rcoefs_zebra_bc[bt,])
  pred_bootz_zebra_bc_dam <- cbind(pred_bootz_zebra_bc_dam, pred_boot)
}
pred_ci_dam_zebra_bc <- t(apply(pred_bootz_zebra_bc_dam, 1, quantile, probs=c(0.025, 0.975)))
pred_df_dam_zebra_bc <- cbind(predict_dam_zebra_bc, pred_dam_zebra_bc, pred_ci_dam_zebra_bc)
colnames(pred_df_dam_zebra_bc) <- c(colnames(predict_dam_zebra_bc), "predictions", "LowerCI", "UpperCI")
pred_df_dam_zebra_bc$dist_dam_km <- pred_df_dam_zebra_bc$dist_dam_clip / 1000

pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_zebra_bc] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")
df_pointz <- cbind(df_pointz, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(df_pointz) <- c("Lat", "Long", "Points", "UTMX", "UTMY", "distX", "distY")

pt_med3 <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of median point used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.05, 0.05)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))
pt_med3 <- ggplot(df_pointz, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(6, 6)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))

plotout_elev_zebra_bc <- ggplot(pred_df_elev_zebra_bc, aes(x=Elev, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Elevation (m)') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_track_zebra_bc <- ggplot(pred_df_track_zebra_bc, aes(x=dist_track_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to track (km)') + geom_ribbon(aes(x=dist_track_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anyrd_zebra_bc <- ggplot(pred_df_anyrd_zebra_bc, aes(x=any_road_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any road (km)') + geom_ribbon(aes(x=any_road_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_minor_zebra_bc <- ggplot(pred_df_minor_zebra_bc, aes(x=dist_minor_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to minor river (km)') + geom_ribbon(aes(x=dist_minor_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anywt_zebra_bc <- ggplot(pred_df_anywt_zebra_bc, aes(x=any_water_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any water (km)') + geom_ribbon(aes(x=any_water_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_dam_zebra_bc <- ggplot(pred_df_dam_zebra_bc, aes(x=dist_dam_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to dam (km)') + geom_ribbon(aes(x=dist_dam_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed3_zebra_bc.png"))
pt_med3
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Elev_zebra.png"))
plotout_elev_zebra_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Track_zebra.png"))
plotout_track_zebra_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Anyrd_zebra.png"))
plotout_anyrd_zebra_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod3minor_zebra.png"))
plotout_minor_zebra_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod3anywater_zebra.png"))
plotout_anywt_zebra_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Dam_zebra.png"))
plotout_dam_zebra_bc
dev.off()


















## not included so graphs not updated


# Fourth model forward selection
var_list_fw <- c("Elev", "dist_track_clip", "any_road_clip", "dist_minor_clip", "any_water_clip", "dist_road_clip")
zebra_mod_fw <- rez_for_var_list(var_list_fw, PA_mod, distMats)

coefz_zebra_fw <- zebra_mod_fw$mod_2d$bestModel$coefficients
covmat_zebra_fw <- summary(zebra_mod_fw$mod_2d$bestModel)$cov.robust
rcoefs_zebra_fw <- rmvnorm(1000, coefz_zebra_fw, sigma=covmat_zebra_fw)

pred_bootz_zebra_fw_all <- predict.gamMRSea(object=zebra_mod_fw$mod_2d$bestModel, coeff=rcoefs_zebra_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=zebra_mod_fw$mod_2d$bestModel, coeff=rcoefs_zebra_fw[bt,])
  pred_bootz_zebra_fw_all <- cbind(pred_bootz_zebra_fw_all, pred_boot)
}
pred_ci_zebra_fw <- t(apply(pred_bootz_zebra_fw_all, 1, quantile, probs=c(0.025, 0.975)))

pred_zebra_fw <- predict.gamMRSea(object=zebra_mod_fw$mod_2d$bestModel)
data_zebra_fw <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_fw, pred_ci_zebra_fw))
colnames(data_zebra_fw) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

graph_pred_zebra_fw <- ggplot(data_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.025,0.05,0.075,0.1,0.15,0.2,1.0), limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_zebra_fw_lw <- ggplot(data_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.025,0.05,0.075,0.1,0.15,0.2,1.0), limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_zebra_fw_hi <- ggplot(data_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.025,0.05,0.075,0.1,0.15,0.2,1.0), limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "ZebraModel4Pred.png"))
graph_pred_zebra_fw
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel4Lw.png"))
graph_pred_zebra_fw_lw
dev.off()

png(filename=paste0(img_save_dir, "ZebraModel4Hi.png"))
graph_pred_zebra_fw_hi
dev.off()


pred_bootz_zebra_fw_all_link <- predict.gamMRSea(object=zebra_mod_fw$mod_2d$bestModel, coeff=rcoefs_zebra_fw[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=zebra_mod_fw$mod_2d$bestModel, coeff=rcoefs_zebra_fw[bt,], type='link')
  pred_bootz_zebra_fw_all_link <- cbind(pred_bootz_zebra_fw_all, pred_boot)
}
pred_ci_zebra_fw_link <- t(apply(pred_bootz_zebra_fw_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_zebra_fw_link <- predict.gamMRSea(object=zebra_mod_fw$mod_2d$bestModel)
data_zebra_fw_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_fw_link, pred_ci_zebra_fw_link))
colnames(data_zebra_fw_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


# filter out any NAs
testcoefs_zebra_fw <- zebra_mod_fw$mod_2d$bestModel$coefficients
cfmsk_zebra_fw <- !is.na(testcoefs_zebra_fw)
tstcfs_zebra_fw <- testcoefs_zebra_fw[cfmsk_zebra_fw]

sp_col_zebra_fw <- zebra_mod_fw$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_zebra_fw <- str_sub(names(tstcfs_zebra_fw),1,strlenclnam)
coefmask_zebra_fw <- coefnamsplt_zebra_fw == nam2d

# create radial gaussian bases
radii_zebra_fw <- sp_col_zebra_fw$radii
radiiInd_zebra_fw <- sp_col_zebra_fw$radiusIndices
aR_zebra_fw <- sp_col_zebra_fw$knotPos
lrf_zebra_fw <- LRF.g(radiiInd_zebra_fw, distMats$dataDist, radii_zebra_fw, aR_zebra_fw)

# combine coefmask and facts
coefz_zebra_fw <- tstcfs_zebra_fw[coefmask_zebra_fw]
# get predicted values on link scale
predtm_zebra_fw <- lrf_zebra_fw %*% coefz_zebra_fw
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_zebra_fw <- rcoefs_zebra_fw[1,][cfmsk_zebra_fw][coefmask_zebra_fw]
predboot_zebra_fw_2d <- lrf_zebra_fw %*% bootcoefz_zebra_fw
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_zebra_fw[bt,][cfmsk_zebra_fw][coefmask_zebra_fw]
  predbt <- lrf_zebra_fw %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_zebra_fw_2d <- cbind(predboot_zebra_fw_2d, predbt)
}
pred_2d_ci_zebra_fw <- t(apply(predboot_zebra_fw_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_zebra_fw <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_zebra_fw, pred_2d_ci_zebra_fw))
colnames(data_2d_zebra_fw) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


zebra_2d_fw <- ggplot(data_2d_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

zebra_2d_fw_lw <- ggplot(data_2d_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

zebra_2d_fw_hi <- ggplot(data_2d_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_zebra_fw <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_zebra_fw_link - predtm_zebra_fw, pred_ci_zebra_fw_link - pred_2d_ci_zebra_fw))
colnames(data_not2d_zebra_fw) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


zebra_not2d_fw <- ggplot(data_not2d_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

zebra_not2d_fw_lw <- ggplot(data_not2d_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

zebra_not2d_fw_hi <- ggplot(data_not2d_zebra_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))


zebra_pred_fw_link <- ggplot(data_zebra_fw_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_fw_lw_link <- ggplot(data_zebra_fw_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

zebra_pred_fw_hi_link <- ggplot(data_zebra_fw_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "ZebraModel42d.png"))
zebra_2d_fw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel42dLw.png"))
zebra_2d_fw_lw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel42dHi.png"))
zebra_2d_fw_hi
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel4NotPrd.png"))
zebra_not2d_fw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel4NotLw.png"))
zebra_not2d_fw_lw
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel4NotHi.png"))
zebra_not2d_fw_hi
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel4Lnk.png"))
zebra_pred_fw_link
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel4LnkLw.png"))
zebra_pred_fw_lw_link
dev.off()
png(filename=paste0(img_save_dir, "ZebraModel4LnkHi.png"))
zebra_pred_fw_hi_link
dev.off()

## 1d calculations
# Elev, dist track, any road, dist minor, any water dam road
# create prediction vals for variables

dist_med_zebra_fw <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_track + dist_med_anyrd + 
  dist_med_minor + dist_med_anywater + dist_med_road
index_med_zebra_fw <- which.min(dist_med_zebra_fw)
median_vals_zebra_fw <- sectionz_ndvi[index_med_zebra_fw,]
distMats_pred_zebra_fw <- makeDists(
  cbind(rep(median_vals_zebra_fw$x.pos, 200), rep(median_vals_zebra_fw$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_zebra_fw <- as.data.frame(cbind(pred_zebra_fw, sectionz_ndvi$Elev, sectionz_ndvi$dist_track_clip, 
                                         sectionz_ndvi$any_road_clip,sectionz_ndvi$dist_minor_clip, 
                                         sectionz_ndvi$any_water_clip, sectionz_ndvi$dist_road_clip))
colnames(rug_vals_zebra_fw) <- c("predictions", "Elev", "dist_track", "any_road", "dist_minor", "any_water", "dist_road")

predict_zebra_fw <- as.data.frame(cbind(rep(median_vals_zebra_fw$x.pos, 200), 
                                        rep(median_vals_zebra_fw$y.pos, 200), 
                                        rep(median_vals_zebra_fw$Elev, 200),
                                        rep(median_vals_zebra_fw$dist_track_clip, 200),
                                        rep(median_vals_zebra_fw$any_road_clip, 200),
                                        rep(median_vals_zebra_fw$dist_minor_clip, 200),
                                        rep(median_vals_zebra_fw$any_water_clip, 200),
                                        rep(median_vals_zebra_fw$dist_road_clip, 200)))
colnames(predict_zebra_fw) <- c("x.pos", "y.pos", "Elev", "dist_track_clip", "any_road_clip", "dist_minor_clip", 
                                "any_water_clip", "dist_road_clip")

predict_elev_zebra_fw <- predict_zebra_fw
predict_elev_zebra_fw$Elev <- elev_predvals
pred_elev_zebra_fw <- predict.gamMRSea(newdata=predict_elev_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist)
pred_bootz_zebra_fw_elev <- predict.gamMRSea(newdata=predict_elev_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[bt,])
  pred_bootz_zebra_fw_elev <- cbind(pred_bootz_zebra_fw_elev, pred_boot)
}
pred_ci_elev_zebra_fw <- t(apply(pred_bootz_zebra_fw_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_zebra_fw <- cbind(predict_elev_zebra_fw, pred_elev_zebra_fw, pred_ci_elev_zebra_fw)
colnames(pred_df_elev_zebra_fw) <- c(colnames(predict_elev_zebra_fw), "predictions", "LowerCI", "UpperCI")

predict_track_zebra_fw <- predict_zebra_fw
predict_track_zebra_fw$dist_track_clip <- track_predvals
pred_track_zebra_fw <- predict.gamMRSea(newdata=predict_track_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist)
pred_bootz_zebra_fw_track <- predict.gamMRSea(newdata=predict_track_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[bt,])
  pred_bootz_zebra_fw_track <- cbind(pred_bootz_zebra_fw_track, pred_boot)
}
pred_ci_track_zebra_fw <- t(apply(pred_bootz_zebra_fw_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_zebra_fw <- cbind(predict_track_zebra_fw, pred_track_zebra_fw, pred_ci_track_zebra_fw)
colnames(pred_df_track_zebra_fw) <- c(colnames(predict_track_zebra_fw), "predictions", "LowerCI", "UpperCI")

predict_anyrd_zebra_fw <- predict_zebra_fw
predict_anyrd_zebra_fw$any_road_clip <- anyrd_predvals
pred_anyrd_zebra_fw <- predict.gamMRSea(newdata=predict_anyrd_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist)
pred_bootz_zebra_fw_anyrd <- predict.gamMRSea(newdata=predict_anyrd_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[bt,])
  pred_bootz_zebra_fw_anyrd <- cbind(pred_bootz_zebra_fw_anyrd, pred_boot)
}
pred_ci_anyrd_zebra_fw <- t(apply(pred_bootz_zebra_fw_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_zebra_fw <- cbind(predict_anyrd_zebra_fw, pred_anyrd_zebra_fw, pred_ci_anyrd_zebra_fw)
colnames(pred_df_anyrd_zebra_fw) <- c(colnames(predict_anyrd_zebra_fw), "predictions", "LowerCI", "UpperCI")

predict_minor_zebra_fw <- predict_zebra_fw
predict_minor_zebra_fw$dist_minor_clip <- minor_predvals
pred_minor_zebra_fw <- predict.gamMRSea(newdata=predict_minor_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist)
pred_bootz_zebra_fw_minor <- predict.gamMRSea(newdata=predict_minor_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_minor_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[bt,])
  pred_bootz_zebra_fw_minor <- cbind(pred_bootz_zebra_fw_minor, pred_boot)
}
pred_ci_minor_zebra_fw <- t(apply(pred_bootz_zebra_fw_minor, 1, quantile, probs=c(0.025, 0.975)))
pred_df_minor_zebra_fw <- cbind(predict_minor_zebra_fw, pred_minor_zebra_fw, pred_ci_minor_zebra_fw)
colnames(pred_df_minor_zebra_fw) <- c(colnames(predict_minor_zebra_fw), "predictions", "LowerCI", "UpperCI")

predict_anywt_zebra_fw <- predict_zebra_fw
predict_anywt_zebra_fw$any_water_clip <- anywater_predvals
pred_anywt_zebra_fw <- predict.gamMRSea(newdata=predict_anywt_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist)
pred_bootz_zebra_fw_anywt <- predict.gamMRSea(newdata=predict_anywt_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[bt,])
  pred_bootz_zebra_fw_anywt <- cbind(pred_bootz_zebra_fw_anywt, pred_boot)
}
pred_ci_anywt_zebra_fw <- t(apply(pred_bootz_zebra_fw_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_zebra_fw <- cbind(predict_anywt_zebra_fw, pred_anywt_zebra_fw, pred_ci_anywt_zebra_fw)
colnames(pred_df_anywt_zebra_fw) <- c(colnames(predict_anywt_zebra_fw), "predictions", "LowerCI", "UpperCI")

predict_road_zebra_fw <- predict_zebra_fw
predict_road_zebra_fw$dist_road_clip <- road_predvals
pred_road_zebra_fw <- predict.gamMRSea(newdata=predict_road_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist)
pred_bootz_zebra_fw_road <- predict.gamMRSea(newdata=predict_road_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_road_zebra_fw, object=zebra_mod_fw$mod_2d$bestModel, g2k=distMats_pred_zebra_fw$dataDist, coeff=rcoefs_zebra_fw[bt,])
  pred_bootz_zebra_fw_road <- cbind(pred_bootz_zebra_fw_road, pred_boot)
}
pred_ci_road_zebra_fw <- t(apply(pred_bootz_zebra_fw_road, 1, quantile, probs=c(0.025, 0.975)))
pred_df_road_zebra_fw <- cbind(predict_road_zebra_fw, pred_road_zebra_fw, pred_ci_road_zebra_fw)
colnames(pred_df_road_zebra_fw) <- c(colnames(predict_road_zebra_fw), "predictions", "LowerCI", "UpperCI")

pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_zebra_fw] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")

pt_med4 <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of median point used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.05, 0.05)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))

plotout_elev_zebra_fw <- ggplot(pred_df_elev_zebra_fw, aes(x=Elev, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Elevation') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_track_zebra_fw <- ggplot(pred_df_track_zebra_fw, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anyrd_zebra_fw <- ggplot(pred_df_anyrd_zebra_fw, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to any road') + geom_ribbon(aes(x=any_road_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_minor_zebra_fw <- ggplot(pred_df_minor_zebra_fw, aes(x=dist_minor_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor river') + geom_ribbon(aes(x=dist_minor_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anywt_zebra_fw <- ggplot(pred_df_anywt_zebra_fw, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to any water') + geom_ribbon(aes(x=any_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_road_zebra_fw <- ggplot(pred_df_road_zebra_fw, aes(x=dist_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to road') + geom_ribbon(aes(x=dist_road_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_zebra_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed4_zebra.png"))
pt_med4
dev.off()
png(filename=paste0(img_save_dir, "PtMod4Elev_zebra.png"))
plotout_elev_zebra_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod4Track_zebra.png"))
plotout_track_zebra_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod4Anyrd_zebra.png"))
plotout_anyrd_zebra_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod4minor_zebra.png"))
plotout_minor_zebra_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod4anywater_zebra.png"))
plotout_anywt_zebra_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod4Road_zebra.png"))
plotout_road_zebra_fw
dev.off()
















