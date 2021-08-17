# This file fits the binomial kudu model and creates output images
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

sectionz_ndvi$response <- sectionz_ndvi$KuduPA
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
kudu_2d_mod <- runSALSA2D(
  PA_mod, 
  salsa2dlist_in, 
  d2k=distMats$dataDist, 
  k2k=distMats$knotDist,
  suppress.printout = TRUE
)
an_res_kudu <- Anova(kudu_2d_mod$bestModel, test='F')
rownames(an_res_kudu) <- c(kudu_2d_mod$bestModel$varshortnames, "2dsmooth", "Residuals")

coefz_kudu_sm <- kudu_2d_mod$bestModel$coefficients
covmat_kudu_sm <- summary(kudu_2d_mod$bestModel)$cov.robust
rcoefs_kudu_sm <- rmvnorm(1000, coefz_kudu_sm, sigma=covmat_kudu_sm)

pred_bootz_kudu_sm_all <- predict.gamMRSea(object=kudu_2d_mod$bestModel, coeff=rcoefs_kudu_sm[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=kudu_2d_mod$bestModel, coeff=rcoefs_kudu_sm[bt,])
  pred_bootz_kudu_sm_all <- cbind(pred_bootz_kudu_sm_all, pred_boot)
}
pred_ci_kudu_sm <- t(apply(pred_bootz_kudu_sm_all, 1, quantile, probs=c(0.025, 0.975)))

pred_kudu_sm <- predict.gamMRSea(object=kudu_2d_mod$bestModel)
data_kudu_sm <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_kudu_sm, pred_ci_kudu_sm))
colnames(data_kudu_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_kudu_sm <- cbind(data_kudu_sm, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_kudu_sm) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")

kudu_pred_sm <- ggplot(data_kudu_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette =  "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
kudu_pred_sm <- ggplot(data_kudu_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_sm_lw <- ggplot(data_kudu_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
kudu_pred_sm_lw <- ggplot(data_kudu_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_sm_hi <- ggplot(data_kudu_sm, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
kudu_pred_sm_hi <- ggplot(data_kudu_sm, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "KuduModel1Pred.png"))
kudu_pred_sm
dev.off()

png(filename=paste0(img_save_dir, "KuduModel1Lw.png"))
kudu_pred_sm_lw
dev.off()

png(filename=paste0(img_save_dir, "KuduModel1Hi.png"))
kudu_pred_sm_hi
dev.off()







# Second model backward selection
var_list_back <- c("NDVI250", "Elev", "dist_fence_clip", "dist_road_clip", "dist_track_clip", "any_road_clip",  
                   "dist_dam_clip", "dist_river_clip", "any_water_clip")
kudu_mod_bc <- rez_for_var_list(var_list_back, PA_mod, distMats)

coefz_kudu_bc <- kudu_mod_bc$mod_2d$bestModel$coefficients
covmat_kudu_bc <- summary(kudu_mod_bc$mod_2d$bestModel)$cov.robust
rcoefs_kudu_bc <- rmvnorm(1000, coefz_kudu_bc, sigma=covmat_kudu_bc)

pred_bootz_kudu_bc_all <- predict.gamMRSea(object=kudu_mod_bc$mod_2d$bestModel, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=kudu_mod_bc$mod_2d$bestModel, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_all <- cbind(pred_bootz_kudu_bc_all, pred_boot)
}
pred_ci_kudu_bc <- t(apply(pred_bootz_kudu_bc_all, 1, quantile, probs=c(0.025, 0.975)))

pred_kudu_bc <- predict.gamMRSea(object=kudu_mod_bc$mod_2d$bestModel)
data_kudu_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_kudu_bc, pred_ci_kudu_bc))
colnames(data_kudu_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")
data_kudu_bc <- cbind(data_kudu_bc, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(data_kudu_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")

graph_pred_kudu_bc <- ggplot(data_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_kudu_bc <- ggplot(data_kudu_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_kudu_bc_lw <- ggplot(data_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette =  "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_kudu_bc_lw <- ggplot(data_kudu_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_kudu_bc_hi <- ggplot(data_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))
graph_pred_kudu_bc_hi <- ggplot(data_kudu_bc, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.1), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "KuduModel2Pred.png"))
graph_pred_kudu_bc
dev.off()

png(filename=paste0(img_save_dir, "KuduModel2Lw.png"))
graph_pred_kudu_bc_lw
dev.off()

png(filename=paste0(img_save_dir, "KuduModel2Hi.png"))
graph_pred_kudu_bc_hi
dev.off()

### not included in thesis

pred_bootz_kudu_bc_all_link <- predict.gamMRSea(object=kudu_mod_bc$mod_2d$bestModel, coeff=rcoefs_kudu_bc[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=kudu_mod_bc$mod_2d$bestModel, coeff=rcoefs_kudu_bc[bt,], type='link')
  pred_bootz_kudu_bc_all_link <- cbind(pred_bootz_kudu_bc_all, pred_boot)
}
pred_ci_kudu_bc_link <- t(apply(pred_bootz_kudu_bc_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_kudu_bc_link <- predict.gamMRSea(object=kudu_mod_bc$mod_2d$bestModel)
data_kudu_bc_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_kudu_bc_link, pred_ci_kudu_bc_link))
colnames(data_kudu_bc_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


# filter out any NAs
testcoefs_kudu_bc <- kudu_mod_bc$mod_2d$bestModel$coefficients
cfmsk_kudu_bc <- !is.na(testcoefs_kudu_bc)
tstcfs_kudu_bc <- testcoefs_kudu_bc[cfmsk_kudu_bc]

sp_col_kudu_bc <- kudu_mod_bc$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_kudu_bc <- str_sub(names(tstcfs_kudu_bc),1,strlenclnam)
coefmask_kudu_bc <- coefnamsplt_kudu_bc == nam2d

# create radial gaussian bases
radii_kudu_bc <- sp_col_kudu_bc$radii
radiiInd_kudu_bc <- sp_col_kudu_bc$radiusIndices
aR_kudu_bc <- sp_col_kudu_bc$knotPos
lrf_kudu_bc <- LRF.g(radiiInd_kudu_bc, distMats$dataDist, radii_kudu_bc, aR_kudu_bc)

# combine coefmask and facts
coefz_kudu_bc <- tstcfs_kudu_bc[coefmask_kudu_bc]
# get predicted values on link scale
predtm_kudu_bc <- lrf_kudu_bc %*% coefz_kudu_bc
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_kudu_bc <- rcoefs_kudu_bc[1,][cfmsk_kudu_bc][coefmask_kudu_bc]
predboot_kudu_bc_2d <- lrf_kudu_bc %*% bootcoefz_kudu_bc
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_kudu_bc[bt,][cfmsk_kudu_bc][coefmask_kudu_bc]
  predbt <- lrf_kudu_bc %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_kudu_bc_2d <- cbind(predboot_kudu_bc_2d, predbt)
}
pred_2d_ci_kudu_bc <- t(apply(predboot_kudu_bc_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_kudu_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_kudu_bc, pred_2d_ci_kudu_bc))
colnames(data_2d_kudu_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


kudu_2d_bc <- ggplot(data_2d_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

kudu_2d_bc_lw <- ggplot(data_2d_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

kudu_2d_bc_hi <- ggplot(data_2d_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_kudu_bc <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_kudu_bc_link - predtm_kudu_bc, pred_ci_kudu_bc_link - pred_2d_ci_kudu_bc))
colnames(data_not2d_kudu_bc) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


kudu_not2d_bc <- ggplot(data_not2d_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

kudu_not2d_bc_lw <- ggplot(data_not2d_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

kudu_not2d_bc_hi <- ggplot(data_not2d_kudu_bc, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))


kudu_pred_bc_link <- ggplot(data_kudu_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_bc_lw_link <- ggplot(data_kudu_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_bc_hi_link <- ggplot(data_kudu_bc_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-10,2), oob=squish, name="") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "KuduModel22d.png"))
  kudu_2d_bc
dev.off()
png(filename=paste0(img_save_dir, "KuduModel22dLw.png"))
  kudu_2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "KuduModel22dHi.png"))
  kudu_2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "KuduModel2NotPrd.png"))
  kudu_not2d_bc
dev.off()
png(filename=paste0(img_save_dir, "KuduModel2NotLw.png"))
  kudu_not2d_bc_lw
dev.off()
png(filename=paste0(img_save_dir, "KuduModel2NotHi.png"))
  kudu_not2d_bc_hi
dev.off()
png(filename=paste0(img_save_dir, "KuduModel2Lnk.png"))
kudu_pred_bc_link
dev.off()
png(filename=paste0(img_save_dir, "KuduModel2LnkLw.png"))
kudu_pred_bc_lw_link
dev.off()
png(filename=paste0(img_save_dir, "KuduModel2LnkHi.png"))
kudu_pred_bc_hi_link
dev.off()

### included in thesis below

## 1d calculations

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
ndvi_predvals <- seq(min(sectionz_ndvi$NDVI250), max(sectionz_ndvi$NDVI250), length.out=200)
fence_predvals <- seq(min(sectionz_ndvi$dist_fence_clip), max(sectionz_ndvi$dist_fence_clip), length.out=200)
track_predvals <- seq(min(sectionz_ndvi$dist_track_clip), max(sectionz_ndvi$dist_track_clip), length.out=200)
anyrd_predvals <- seq(min(sectionz_ndvi$any_road_clip), max(sectionz_ndvi$any_road_clip), length.out=200)
dam_predvals <- seq(min(sectionz_ndvi$dist_dam_clip), max(sectionz_ndvi$dist_dam_clip), length.out=200)
anywater_predvals <- seq(min(sectionz_ndvi$any_water_clip), max(sectionz_ndvi$any_water_clip), length.out=200)
road_predvals <- seq(min(sectionz_ndvi$dist_road_clip), max(sectionz_ndvi$dist_road_clip), length.out=200)
river_predvals <- seq(min(sectionz_ndvi$dist_river_clip), max(sectionz_ndvi$dist_river_clip), length.out=200)


# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
elev_normal <- (sectionz_ndvi$Elev - mean(sectionz_ndvi$Elev)) / sd(sectionz_ndvi$Elev)
ndvi_normal <- (sectionz_ndvi$NDVI250 - mean(sectionz_ndvi$NDVI250)) / sd(sectionz_ndvi$NDVI250)
fence_normal <- (sectionz_ndvi$dist_fence_clip - mean(sectionz_ndvi$dist_fence_clip)) / sd(sectionz_ndvi$dist_fence_clip)
track_normal <- (sectionz_ndvi$dist_track_clip - mean(sectionz_ndvi$dist_track_clip)) / sd(sectionz_ndvi$dist_track_clip)
anyrd_normal <- (sectionz_ndvi$any_road_clip - mean(sectionz_ndvi$any_road_clip)) / sd(sectionz_ndvi$any_road_clip)
dam_normal <- (sectionz_ndvi$dist_dam_clip - mean(sectionz_ndvi$dist_dam_clip)) / sd(sectionz_ndvi$dist_dam_clip)
anywater_normal <- (sectionz_ndvi$any_water_clip - mean(sectionz_ndvi$any_water_clip)) / sd(sectionz_ndvi$any_water_clip)
road_normal <- (sectionz_ndvi$dist_road_clip - mean(sectionz_ndvi$dist_road_clip)) / sd(sectionz_ndvi$dist_road_clip)
river_normal <- (sectionz_ndvi$dist_river_clip - mean(sectionz_ndvi$dist_river_clip)) / sd(sectionz_ndvi$dist_river_clip)

dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_elev <- (elev_normal - median(elev_normal))^2
dist_med_ndvi <- (ndvi_normal - median(ndvi_normal))^2
dist_med_fence <- (fence_normal - median(fence_normal))^2
dist_med_track <- (track_normal - median(track_normal))^2
dist_med_anyrd <- (anyrd_normal - median(anyrd_normal))^2
dist_med_dam <- (dam_normal - median(dam_normal))^2
dist_med_anywater <- (anywater_normal - median(anywater_normal))^2
dist_med_road <- (road_normal - median(road_normal))^2
dist_med_river <- (river_normal - median(river_normal))^2



dist_med_kudu_bc <- dist_med_xpos + dist_med_ypos + dist_med_elev + dist_med_ndvi + dist_med_fence + 
  dist_med_track + dist_med_anyrd + dist_med_dam + dist_med_anywater + dist_med_road + dist_med_river
index_med_kudu_bc <- which.min(dist_med_kudu_bc)
median_vals_kudu_bc <- sectionz_ndvi[index_med_kudu_bc,]
distMats_pred_kudu_bc <- makeDists(
  cbind(rep(median_vals_kudu_bc$x.pos, 200), rep(median_vals_kudu_bc$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_kudu_bc <- as.data.frame(cbind(pred_kudu_bc, sectionz_ndvi$Elev, sectionz_ndvi$NDVI250, 
                                   sectionz_ndvi$dist_fence_clip, sectionz_ndvi$dist_track_clip, 
                                   sectionz_ndvi$any_road_clip, 
                                   sectionz_ndvi$dist_dam_clip, sectionz_ndvi$any_water_clip,
                                   sectionz_ndvi$dist_road_clip, sectionz_ndvi$dist_river_clip))
colnames(rug_vals_kudu_bc) <- c("predictions", "Elev", "NDVI250", "dist_fence", "dist_track",
                               "any_road", "dist_dam", "any_water", "dist_road", "dist_river")
rug_vals_kudu_bc$dist_fence_km <- rug_vals_kudu_bc$dist_fence / 1000
rug_vals_kudu_bc$dist_track_km <- rug_vals_kudu_bc$dist_track / 1000
rug_vals_kudu_bc$any_road_km <- rug_vals_kudu_bc$any_road / 1000
rug_vals_kudu_bc$dist_dam_km <- rug_vals_kudu_bc$dist_dam / 1000
rug_vals_kudu_bc$any_water_km <- rug_vals_kudu_bc$any_water / 1000
rug_vals_kudu_bc$dist_road_km <- rug_vals_kudu_bc$dist_road / 1000
rug_vals_kudu_bc$dist_river_km <- rug_vals_kudu_bc$dist_river / 1000

predict_kudu_bc <- as.data.frame(cbind(rep(median_vals_kudu_bc$x.pos, 200), 
                                          rep(median_vals_kudu_bc$y.pos, 200), 
                                          rep(median_vals_kudu_bc$Elev, 200),
                                          rep(median_vals_kudu_bc$NDVI250, 200),
                                          rep(median_vals_kudu_bc$dist_fence_clip, 200),
                                          rep(median_vals_kudu_bc$dist_track_clip, 200),
                                          rep(median_vals_kudu_bc$any_road_clip, 200),
                                          rep(median_vals_kudu_bc$dist_dam_clip, 200),
                                          rep(median_vals_kudu_bc$any_water_clip, 200),
                                          rep(median_vals_kudu_bc$dist_road_clip, 200),
                                          rep(median_vals_kudu_bc$dist_river_clip, 200)))
colnames(predict_kudu_bc) <- c("x.pos", "y.pos", "Elev", "NDVI250", "dist_fence_clip", 
                                  "dist_track_clip", "any_road_clip",
                                  "dist_dam_clip", "any_water_clip", "dist_road_clip", "dist_river_clip")

predict_elev_kudu_bc <- predict_kudu_bc
predict_elev_kudu_bc$Elev <- elev_predvals
pred_elev_kudu_bc <- predict.gamMRSea(newdata=predict_elev_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_elev <- predict.gamMRSea(newdata=predict_elev_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elev_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_elev <- cbind(pred_bootz_kudu_bc_elev, pred_boot)
}
pred_ci_elev_kudu_bc <- t(apply(pred_bootz_kudu_bc_elev, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elev_kudu_bc <- cbind(predict_elev_kudu_bc, pred_elev_kudu_bc, pred_ci_elev_kudu_bc)
colnames(pred_df_elev_kudu_bc) <- c(colnames(predict_elev_kudu_bc), "predictions", "LowerCI", "UpperCI")

predict_ndvi_kudu_bc <- predict_kudu_bc
predict_ndvi_kudu_bc$NDVI250 <- ndvi_predvals
pred_ndvi_kudu_bc <- predict.gamMRSea(newdata=predict_ndvi_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_ndvi <- predict.gamMRSea(newdata=predict_ndvi_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_ndvi_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_ndvi <- cbind(pred_bootz_kudu_bc_ndvi, pred_boot)
}
pred_ci_ndvi_kudu_bc <- t(apply(pred_bootz_kudu_bc_ndvi, 1, quantile, probs=c(0.025, 0.975)))
pred_df_ndvi_kudu_bc <- cbind(predict_ndvi_kudu_bc, pred_ndvi_kudu_bc, pred_ci_ndvi_kudu_bc)
colnames(pred_df_ndvi_kudu_bc) <- c(colnames(predict_ndvi_kudu_bc), "predictions", "LowerCI", "UpperCI")

predict_fence_kudu_bc <- predict_kudu_bc
predict_fence_kudu_bc$dist_fence_clip <- fence_predvals
pred_fence_kudu_bc <- predict.gamMRSea(newdata=predict_fence_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_fence <- predict.gamMRSea(newdata=predict_fence_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_fence_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_fence <- cbind(pred_bootz_kudu_bc_fence, pred_boot)
}
pred_ci_fence_kudu_bc <- t(apply(pred_bootz_kudu_bc_fence, 1, quantile, probs=c(0.025, 0.975)))
pred_df_fence_kudu_bc <- cbind(predict_fence_kudu_bc, pred_fence_kudu_bc, pred_ci_fence_kudu_bc)
colnames(pred_df_fence_kudu_bc) <- c(colnames(predict_fence_kudu_bc), "predictions", "LowerCI", "UpperCI")
pred_df_fence_kudu_bc$dist_fence_km <- pred_df_fence_kudu_bc$dist_fence_clip / 1000

predict_track_kudu_bc <- predict_kudu_bc
predict_track_kudu_bc$dist_track_clip <- track_predvals
pred_track_kudu_bc <- predict.gamMRSea(newdata=predict_track_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_track <- predict.gamMRSea(newdata=predict_track_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_track <- cbind(pred_bootz_kudu_bc_track, pred_boot)
}
pred_ci_track_kudu_bc <- t(apply(pred_bootz_kudu_bc_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_kudu_bc <- cbind(predict_track_kudu_bc, pred_track_kudu_bc, pred_ci_track_kudu_bc)
colnames(pred_df_track_kudu_bc) <- c(colnames(predict_track_kudu_bc), "predictions", "LowerCI", "UpperCI")
pred_df_track_kudu_bc$dist_track_km <- pred_df_track_kudu_bc$dist_track_clip / 1000

predict_anyrd_kudu_bc <- predict_kudu_bc
predict_anyrd_kudu_bc$any_road_clip <- anyrd_predvals
pred_anyrd_kudu_bc <- predict.gamMRSea(newdata=predict_anyrd_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_anyrd <- predict.gamMRSea(newdata=predict_anyrd_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_anyrd <- cbind(pred_bootz_kudu_bc_anyrd, pred_boot)
}
pred_ci_anyrd_kudu_bc <- t(apply(pred_bootz_kudu_bc_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_kudu_bc <- cbind(predict_anyrd_kudu_bc, pred_anyrd_kudu_bc, pred_ci_anyrd_kudu_bc)
colnames(pred_df_anyrd_kudu_bc) <- c(colnames(predict_anyrd_kudu_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anyrd_kudu_bc$any_road_km <- pred_df_anyrd_kudu_bc$any_road_clip / 1000

predict_dam_kudu_bc <- predict_kudu_bc
predict_dam_kudu_bc$dist_dam_clip <- dam_predvals
pred_dam_kudu_bc <- predict.gamMRSea(newdata=predict_dam_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_dam <- predict.gamMRSea(newdata=predict_dam_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_dam_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_dam <- cbind(pred_bootz_kudu_bc_dam, pred_boot)
}
pred_ci_dam_kudu_bc <- t(apply(pred_bootz_kudu_bc_dam, 1, quantile, probs=c(0.025, 0.975)))
pred_df_dam_kudu_bc <- cbind(predict_dam_kudu_bc, pred_dam_kudu_bc, pred_ci_dam_kudu_bc)
colnames(pred_df_dam_kudu_bc) <- c(colnames(predict_dam_kudu_bc), "predictions", "LowerCI", "UpperCI")
pred_df_dam_kudu_bc$dist_dam_km <- pred_df_dam_kudu_bc$dist_dam_clip / 1000

predict_anywt_kudu_bc <- predict_kudu_bc
predict_anywt_kudu_bc$any_water_clip <- anywater_predvals
pred_anywt_kudu_bc <- predict.gamMRSea(newdata=predict_anywt_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_anywt <- predict.gamMRSea(newdata=predict_anywt_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_anywt <- cbind(pred_bootz_kudu_bc_anywt, pred_boot)
}
pred_ci_anywt_kudu_bc <- t(apply(pred_bootz_kudu_bc_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_kudu_bc <- cbind(predict_anywt_kudu_bc, pred_anywt_kudu_bc, pred_ci_anywt_kudu_bc)
colnames(pred_df_anywt_kudu_bc) <- c(colnames(predict_anywt_kudu_bc), "predictions", "LowerCI", "UpperCI")
pred_df_anywt_kudu_bc$any_water_km <- pred_df_anywt_kudu_bc$any_water_clip / 1000

predict_road_kudu_bc <- predict_kudu_bc
predict_road_kudu_bc$dist_road_clip <- road_predvals
pred_road_kudu_bc <- predict.gamMRSea(newdata=predict_road_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_road <- predict.gamMRSea(newdata=predict_road_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_road_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_road <- cbind(pred_bootz_kudu_bc_road, pred_boot)
}
pred_ci_road_kudu_bc <- t(apply(pred_bootz_kudu_bc_road, 1, quantile, probs=c(0.025, 0.975)))
pred_df_road_kudu_bc <- cbind(predict_road_kudu_bc, pred_road_kudu_bc, pred_ci_road_kudu_bc)
colnames(pred_df_road_kudu_bc) <- c(colnames(predict_road_kudu_bc), "predictions", "LowerCI", "UpperCI")
pred_df_road_kudu_bc$dist_road_km <- pred_df_road_kudu_bc$dist_road_clip / 1000

predict_river_kudu_bc <- predict_kudu_bc
predict_river_kudu_bc$dist_river_clip <- river_predvals
pred_river_kudu_bc <- predict.gamMRSea(newdata=predict_river_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist)
pred_bootz_kudu_bc_river <- predict.gamMRSea(newdata=predict_river_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_river_kudu_bc, object=kudu_mod_bc$mod_2d$bestModel, g2k=distMats_pred_kudu_bc$dataDist, coeff=rcoefs_kudu_bc[bt,])
  pred_bootz_kudu_bc_river <- cbind(pred_bootz_kudu_bc_river, pred_boot)
}
pred_ci_river_kudu_bc <- t(apply(pred_bootz_kudu_bc_river, 1, quantile, probs=c(0.025, 0.975)))
pred_df_river_kudu_bc <- cbind(predict_river_kudu_bc, pred_river_kudu_bc, pred_ci_river_kudu_bc)
colnames(pred_df_river_kudu_bc) <- c(colnames(predict_river_kudu_bc), "predictions", "LowerCI", "UpperCI")
pred_df_river_kudu_bc$dist_river_km <- pred_df_river_kudu_bc$dist_river_clip / 1000

pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_kudu_bc] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")
df_pointz <- cbind(df_pointz, (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170)
colnames(df_pointz) <- c("Lat", "Long", "Points", "UTMX", "UTMY")

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


plotout_elev_kudu_bc <- ggplot(pred_df_elev_kudu_bc, aes(x=Elev, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Elevation (m)') + geom_ribbon(aes(x=Elev, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_ndvi_kudu_bc <- ggplot(pred_df_ndvi_kudu_bc, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_fence_kudu_bc <- ggplot(pred_df_fence_kudu_bc, aes(x=dist_fence_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to fence (km)') + geom_ribbon(aes(x=dist_fence_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_track_kudu_bc <- ggplot(pred_df_track_kudu_bc, aes(x=dist_track_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to track (km)') + geom_ribbon(aes(x=dist_track_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anyrd_kudu_bc <- ggplot(pred_df_anyrd_kudu_bc, aes(x=any_road_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any road (km)') + geom_ribbon(aes(x=any_road_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_dam_kudu_bc <- ggplot(pred_df_dam_kudu_bc, aes(x=dist_dam_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to dam (km)') + geom_ribbon(aes(x=dist_dam_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anywt_kudu_bc <- ggplot(pred_df_anywt_kudu_bc, aes(x=any_water_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to any water (km)') + geom_ribbon(aes(x=any_water_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_road_kudu_bc <- ggplot(pred_df_road_kudu_bc, aes(x=dist_road_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to road (km)') + geom_ribbon(aes(x=dist_road_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_river_kudu_bc <- ggplot(pred_df_river_kudu_bc, aes(x=dist_river_km, y=predictions)) + theme_bw() + ylab('Predicted probability') + xlab('Distance to river (km)') + geom_ribbon(aes(x=dist_river_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_bc, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed2_kudu.png"))
pt_med2
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Elev_kudu.png"))
plotout_elev_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2NDVI_kudu.png"))
plotout_ndvi_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Fence_kudu.png"))
plotout_fence_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Track_kudu.png"))
plotout_track_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Anyrd_kudu.png"))
plotout_anyrd_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2dam_kudu.png"))
plotout_dam_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2anywater_kudu.png"))
plotout_anywt_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2road_kudu.png"))
plotout_road_kudu_bc
dev.off()
png(filename=paste0(img_save_dir, "PtMod2river_kudu.png"))
plotout_river_kudu_bc
dev.off()





## forward selection models not included graphs not regenrated
var_list_fw <- c("NDVI250", "dist_fence_clip", "dist_track_clip", "any_road_clip", "dist_dam_clip",  "any_water_clip",  "dist_road_clip", "dist_river_clip","dist_elec_clip")
kudu_mod_fw <- rez_for_var_list(var_list_fw, PA_mod, distMats)


coefz_kudu_fw <- kudu_mod_fw$mod_2d$bestModel$coefficients
covmat_kudu_fw <- summary(kudu_mod_fw$mod_2d$bestModel)$cov.robust
rcoefs_kudu_fw <- rmvnorm(1000, coefz_kudu_fw, sigma=covmat_kudu_fw)

pred_bootz_kudu_fw_all <- predict.gamMRSea(object=kudu_mod_fw$mod_2d$bestModel, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=kudu_mod_fw$mod_2d$bestModel, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_all <- cbind(pred_bootz_kudu_fw_all, pred_boot)
}
pred_ci_kudu_fw <- t(apply(pred_bootz_kudu_fw_all, 1, quantile, probs=c(0.025, 0.975)))


pred_kudu_fw <- predict.gamMRSea(object=kudu_mod_fw$mod_2d$bestModel)
data_kudu_fw <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_kudu_fw, pred_ci_kudu_fw))
colnames(data_kudu_fw) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

graph_pred_kudu_fw <- ggplot(data_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.025,0.05,0.075,0.1,0.15,0.2,1.0), limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_kudu_fw_lw <- ggplot(data_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.025,0.05,0.075,0.1,0.15,0.2,1.0), limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_pred_kudu_fw_hi <- ggplot(data_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu", values=c(0,0.025,0.05,0.075,0.1,0.15,0.2,1.0), limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "KuduModel3Pred.png"))
graph_pred_kudu_fw
dev.off()

png(filename=paste0(img_save_dir, "KuduModel3Lw.png"))
graph_pred_kudu_fw_lw
dev.off()

png(filename=paste0(img_save_dir, "KuduModel3Hi.png"))
graph_pred_kudu_fw_hi
dev.off()


pred_bootz_kudu_fw_all_link <- predict.gamMRSea(object=kudu_mod_fw$mod_2d$bestModel, coeff=rcoefs_kudu_fw[1,], type='link')
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(object=kudu_mod_fw$mod_2d$bestModel, coeff=rcoefs_kudu_fw[bt,], type='link')
  pred_bootz_kudu_fw_all_link <- cbind(pred_bootz_kudu_fw_all, pred_boot)
}
pred_ci_kudu_fw_link <- t(apply(pred_bootz_kudu_fw_all_link, 1, quantile, probs=c(0.025, 0.975)))


pred_kudu_fw_link <- predict.gamMRSea(object=kudu_mod_fw$mod_2d$bestModel)
data_kudu_fw_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_kudu_fw_link, pred_ci_kudu_fw_link))
colnames(data_kudu_fw_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")




# filter out any NAs
testcoefs_kudu_fw <- kudu_mod_fw$mod_2d$bestModel$coefficients
cfmsk_kudu_fw <- !is.na(testcoefs_kudu_fw)
tstcfs_kudu_fw <- testcoefs_kudu_fw[cfmsk_kudu_fw]

sp_col_kudu_fw <- kudu_mod_fw$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt_kudu_fw <- str_sub(names(tstcfs_kudu_fw),1,strlenclnam)
coefmask_kudu_fw <- coefnamsplt_kudu_fw == nam2d

# create radial gaussian bases
radii_kudu_fw <- sp_col_kudu_fw$radii
radiiInd_kudu_fw <- sp_col_kudu_fw$radiusIndices
aR_kudu_fw <- sp_col_kudu_fw$knotPos
lrf_kudu_fw <- LRF.g(radiiInd_kudu_fw, distMats$dataDist, radii_kudu_fw, aR_kudu_fw)

# combine coefmask and facts
coefz_kudu_fw <- tstcfs_kudu_fw[coefmask_kudu_fw]
# get predicted values on link scale
predtm_kudu_fw <- lrf_kudu_fw %*% coefz_kudu_fw
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

bootcoefz_kudu_fw <- rcoefs_kudu_fw[1,][cfmsk_kudu_fw][coefmask_kudu_fw]
predboot_kudu_fw_2d <- lrf_kudu_fw %*% bootcoefz_kudu_fw
# predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
for (bt in 2:1000){
  bootcf <- rcoefs_kudu_fw[bt,][cfmsk_kudu_fw][coefmask_kudu_fw]
  predbt <- lrf_kudu_fw %*% bootcf
  # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
  predboot_kudu_fw_2d <- cbind(predboot_kudu_fw_2d, predbt)
}
pred_2d_ci_kudu_fw <- t(apply(predboot_kudu_fw_2d, 1, quantile, probs=c(0.025, 0.975)))


data_2d_kudu_fw <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm_kudu_fw, pred_2d_ci_kudu_fw))
colnames(data_2d_kudu_fw) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


kudu_2d_fw <- ggplot(data_2d_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

kudu_2d_fw_lw <- ggplot(data_2d_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

kudu_2d_fw_hi <- ggplot(data_2d_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_kudu_fw <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_kudu_fw_link - predtm_kudu_fw, pred_ci_kudu_fw_link - pred_2d_ci_kudu_fw))
colnames(data_not2d_kudu_fw) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")


kudu_not2d_fw <- ggplot(data_not2d_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

kudu_not2d_fw_lw <- ggplot(data_not2d_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

kudu_not2d_fw_hi <- ggplot(data_not2d_kudu_fw, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))


kudu_pred_fw_link <- ggplot(data_kudu_fw_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_fw_lw_link <- ggplot(data_kudu_fw_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

kudu_pred_fw_hi_link <- ggplot(data_kudu_fw_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "KuduModel32d.png"))
kudu_2d_fw
dev.off()
png(filename=paste0(img_save_dir, "KuduModel32dLw.png"))
kudu_2d_fw_lw
dev.off()
png(filename=paste0(img_save_dir, "KuduModel32dHi.png"))
kudu_2d_fw_hi
dev.off()
png(filename=paste0(img_save_dir, "KuduModel3NotPrd.png"))
kudu_not2d_fw
dev.off()
png(filename=paste0(img_save_dir, "KuduModel3NotLw.png"))
kudu_not2d_fw_lw
dev.off()
png(filename=paste0(img_save_dir, "KuduModel3NotHi.png"))
kudu_not2d_fw_hi
dev.off()
png(filename=paste0(img_save_dir, "KuduModel3Lnk.png"))
kudu_pred_fw_link
dev.off()
png(filename=paste0(img_save_dir, "KuduModel3LnkLw.png"))
kudu_pred_fw_lw_link
dev.off()
png(filename=paste0(img_save_dir, "KuduModel3LnkHi.png"))
kudu_pred_fw_hi_link
dev.off()

## 1d calculations

# create prediction vals for variables
elec_predvals <- seq(min(sectionz_ndvi$dist_elec_clip), max(sectionz_ndvi$dist_elec_clip), length.out=200)
#ndvi, elev, fence, track, any road, wat, dam, any water

elec_normal <- (sectionz_ndvi$dist_elec_clip - mean(sectionz_ndvi$dist_elec_clip)) / sd(sectionz_ndvi$dist_elec_clip)

dist_med_elec <- (elec_normal - median(elec_normal))^2

dist_med_kudu_fw <- dist_med_xpos + dist_med_ypos + dist_med_ndvi + dist_med_fence + 
  dist_med_track + dist_med_anyrd + dist_med_dam + dist_med_anywater + dist_med_road + dist_med_river + dist_med_elec
index_med_kudu_fw <- which.min(dist_med_kudu_fw)
median_vals_kudu_fw <- sectionz_ndvi[index_med_kudu_fw,]
distMats_pred_kudu_fw <- makeDists(
  cbind(rep(median_vals_kudu_fw$x.pos, 200), rep(median_vals_kudu_fw$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_kudu_fw <- as.data.frame(cbind(pred_kudu_fw, sectionz_ndvi$NDVI250, 
                                        sectionz_ndvi$dist_fence_clip, sectionz_ndvi$dist_track_clip, 
                                        sectionz_ndvi$any_road_clip, 
                                        sectionz_ndvi$dist_dam_clip, sectionz_ndvi$any_water_clip,
                                        sectionz_ndvi$dist_road_clip, sectionz_ndvi$dist_river_clip,
                                        sectionz_ndvi$dist_elec_clip))
colnames(rug_vals_kudu_fw) <- c("predictions", "NDVI250", "dist_fence", "dist_track",
                                "any_road", "dist_dam", "any_water", "dist_road", "dist_river", "dist_elec")


predict_kudu_fw <- as.data.frame(cbind(rep(median_vals_kudu_fw$x.pos, 200), 
                                       rep(median_vals_kudu_fw$y.pos, 200), 
                                       rep(median_vals_kudu_fw$NDVI250, 200),
                                       rep(median_vals_kudu_fw$dist_fence_clip, 200),
                                       rep(median_vals_kudu_fw$dist_track_clip, 200),
                                       rep(median_vals_kudu_fw$any_road_clip, 200),
                                       rep(median_vals_kudu_fw$dist_dam_clip, 200),
                                       rep(median_vals_kudu_fw$any_water_clip, 200),
                                       rep(median_vals_kudu_fw$dist_road_clip, 200),
                                       rep(median_vals_kudu_fw$dist_river_clip, 200),
                                       rep(median_vals_kudu_fw$dist_elec_clip, 200)))
colnames(predict_kudu_fw) <- c("x.pos", "y.pos", "NDVI250", "dist_fence_clip", "dist_track_clip", "any_road_clip",
                               "dist_dam_clip", "any_water_clip", "dist_road_clip", "dist_river_clip", "dist_elec_clip")

predict_elec_kudu_fw <- predict_kudu_fw
predict_elec_kudu_fw$dist_elec_clip <- elec_predvals
pred_elec_kudu_fw <- predict.gamMRSea(newdata=predict_elec_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_elec <- predict.gamMRSea(newdata=predict_elec_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_elec_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_elec <- cbind(pred_bootz_kudu_fw_elec, pred_boot)
}
pred_ci_elec_kudu_fw <- t(apply(pred_bootz_kudu_fw_elec, 1, quantile, probs=c(0.025, 0.975)))
pred_df_elec_kudu_fw <- cbind(predict_elec_kudu_fw, pred_elec_kudu_fw, pred_ci_elec_kudu_fw)
colnames(pred_df_elec_kudu_fw) <- c(colnames(predict_elec_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_ndvi_kudu_fw <- predict_kudu_fw
predict_ndvi_kudu_fw$NDVI250 <- ndvi_predvals
pred_ndvi_kudu_fw <- predict.gamMRSea(newdata=predict_ndvi_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_ndvi <- predict.gamMRSea(newdata=predict_ndvi_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_ndvi_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_ndvi <- cbind(pred_bootz_kudu_fw_ndvi, pred_boot)
}
pred_ci_ndvi_kudu_fw <- t(apply(pred_bootz_kudu_fw_ndvi, 1, quantile, probs=c(0.025, 0.975)))
pred_df_ndvi_kudu_fw <- cbind(predict_ndvi_kudu_fw, pred_ndvi_kudu_fw, pred_ci_ndvi_kudu_fw)
colnames(pred_df_ndvi_kudu_fw) <- c(colnames(predict_ndvi_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_fence_kudu_fw <- predict_kudu_fw
predict_fence_kudu_fw$dist_fence_clip <- fence_predvals
pred_fence_kudu_fw <- predict.gamMRSea(newdata=predict_fence_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_fence <- predict.gamMRSea(newdata=predict_fence_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_fence_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_fence <- cbind(pred_bootz_kudu_fw_fence, pred_boot)
}
pred_ci_fence_kudu_fw <- t(apply(pred_bootz_kudu_fw_fence, 1, quantile, probs=c(0.025, 0.975)))
pred_df_fence_kudu_fw <- cbind(predict_fence_kudu_fw, pred_fence_kudu_fw, pred_ci_fence_kudu_fw)
colnames(pred_df_fence_kudu_fw) <- c(colnames(predict_fence_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_track_kudu_fw <- predict_kudu_fw
predict_track_kudu_fw$dist_track_clip <- track_predvals
pred_track_kudu_fw <- predict.gamMRSea(newdata=predict_track_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_track <- predict.gamMRSea(newdata=predict_track_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_track_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_track <- cbind(pred_bootz_kudu_fw_track, pred_boot)
}
pred_ci_track_kudu_fw <- t(apply(pred_bootz_kudu_fw_track, 1, quantile, probs=c(0.025, 0.975)))
pred_df_track_kudu_fw <- cbind(predict_track_kudu_fw, pred_track_kudu_fw, pred_ci_track_kudu_fw)
colnames(pred_df_track_kudu_fw) <- c(colnames(predict_track_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_anyrd_kudu_fw <- predict_kudu_fw
predict_anyrd_kudu_fw$any_road_clip <- anyrd_predvals
pred_anyrd_kudu_fw <- predict.gamMRSea(newdata=predict_anyrd_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_anyrd <- predict.gamMRSea(newdata=predict_anyrd_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anyrd_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_anyrd <- cbind(pred_bootz_kudu_fw_anyrd, pred_boot)
}
pred_ci_anyrd_kudu_fw <- t(apply(pred_bootz_kudu_fw_anyrd, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anyrd_kudu_fw <- cbind(predict_anyrd_kudu_fw, pred_anyrd_kudu_fw, pred_ci_anyrd_kudu_fw)
colnames(pred_df_anyrd_kudu_fw) <- c(colnames(predict_anyrd_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_dam_kudu_fw <- predict_kudu_fw
predict_dam_kudu_fw$dist_dam_clip <- dam_predvals
pred_dam_kudu_fw <- predict.gamMRSea(newdata=predict_dam_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_dam <- predict.gamMRSea(newdata=predict_dam_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_dam_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_dam <- cbind(pred_bootz_kudu_fw_dam, pred_boot)
}
pred_ci_dam_kudu_fw <- t(apply(pred_bootz_kudu_fw_dam, 1, quantile, probs=c(0.025, 0.975)))
pred_df_dam_kudu_fw <- cbind(predict_dam_kudu_fw, pred_dam_kudu_fw, pred_ci_dam_kudu_fw)
colnames(pred_df_dam_kudu_fw) <- c(colnames(predict_dam_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_anywt_kudu_fw <- predict_kudu_fw
predict_anywt_kudu_fw$any_water_clip <- anywater_predvals
pred_anywt_kudu_fw <- predict.gamMRSea(newdata=predict_anywt_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_anywt <- predict.gamMRSea(newdata=predict_anywt_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_anywt_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_anywt <- cbind(pred_bootz_kudu_fw_anywt, pred_boot)
}
pred_ci_anywt_kudu_fw <- t(apply(pred_bootz_kudu_fw_anywt, 1, quantile, probs=c(0.025, 0.975)))
pred_df_anywt_kudu_fw <- cbind(predict_anywt_kudu_fw, pred_anywt_kudu_fw, pred_ci_anywt_kudu_fw)
colnames(pred_df_anywt_kudu_fw) <- c(colnames(predict_anywt_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_road_kudu_fw <- predict_kudu_fw
predict_road_kudu_fw$dist_road_clip <- road_predvals
pred_road_kudu_fw <- predict.gamMRSea(newdata=predict_road_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_road <- predict.gamMRSea(newdata=predict_road_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_road_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_road <- cbind(pred_bootz_kudu_fw_road, pred_boot)
}
pred_ci_road_kudu_fw <- t(apply(pred_bootz_kudu_fw_road, 1, quantile, probs=c(0.025, 0.975)))
pred_df_road_kudu_fw <- cbind(predict_road_kudu_fw, pred_road_kudu_fw, pred_ci_road_kudu_fw)
colnames(pred_df_road_kudu_fw) <- c(colnames(predict_road_kudu_fw), "predictions", "LowerCI", "UpperCI")

predict_river_kudu_fw <- predict_kudu_fw
predict_river_kudu_fw$dist_river_clip <- riverr_predvals
pred_river_kudu_fw <- predict.gamMRSea(newdata=predict_river_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist)
pred_bootz_kudu_fw_river <- predict.gamMRSea(newdata=predict_river_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_river_kudu_fw, object=kudu_mod_fw$mod_2d$bestModel, g2k=distMats_pred_kudu_fw$dataDist, coeff=rcoefs_kudu_fw[bt,])
  pred_bootz_kudu_fw_river <- cbind(pred_bootz_kudu_fw_river, pred_boot)
}
pred_ci_river_kudu_fw <- t(apply(pred_bootz_kudu_fw_river, 1, quantile, probs=c(0.025, 0.975)))
pred_df_river_kudu_fw <- cbind(predict_river_kudu_fw, pred_river_kudu_fw, pred_ci_river_kudu_fw)
colnames(pred_df_river_kudu_fw) <- c(colnames(predict_river_kudu_fw), "predictions", "LowerCI", "UpperCI")


pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_kudu_fw] <- 1

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")

pt_med3 <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of median point used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.05, 0.05)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,1)) +
  theme(legend.position="right", text = element_text(size=20))

plotout_elec_kudu_fw <- ggplot(pred_df_elec_kudu_fw, aes(x=dist_elec_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to electricity line') + geom_ribbon(aes(x=dist_elec_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_ndvi_kudu_fw <- ggplot(pred_df_ndvi_kudu_fw, aes(x=NDVI250, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('NDVI') + geom_ribbon(aes(x=NDVI250, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_fence_kudu_fw <- ggplot(pred_df_fence_kudu_fw, aes(x=dist_fence_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to fence') + geom_ribbon(aes(x=dist_fence_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_track_kudu_fw <- ggplot(pred_df_track_kudu_fw, aes(x=dist_track_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to track') + geom_ribbon(aes(x=dist_track_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anyrd_kudu_fw <- ggplot(pred_df_anyrd_kudu_fw, aes(x=any_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to any road') + geom_ribbon(aes(x=any_road_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_dam_kudu_fw <- ggplot(pred_df_dam_kudu_fw, aes(x=dist_dam_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to dam') + geom_ribbon(aes(x=dist_dam_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_anywt_kudu_fw <- ggplot(pred_df_anywt_kudu_fw, aes(x=any_water_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to any water') + geom_ribbon(aes(x=any_water_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_road_kudu_fw <- ggplot(pred_df_road_kudu_fw, aes(x=dist_road_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to road') + geom_ribbon(aes(x=dist_road_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
plotout_river_kudu_fw <- ggplot(pred_df_river_kudu_fw, aes(x=dist_river_clip, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to river') + geom_ribbon(aes(x=dist_river_clip, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_kudu_fw, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))


png(filename=paste0(img_save_dir, "PtMed3_kudu.png"))
pt_med2
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Elec_kudu.png"))
plotout_elec_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3NDVI_kudu.png"))
plotout_ndvi_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Fence_kudu.png"))
plotout_fence_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Track_kudu.png"))
plotout_track_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3Anyrd_kudu.png"))
plotout_anyrd_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3dam_kudu.png"))
plotout_dam_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3anywater_kudu.png"))
plotout_anywt_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3road_kudu.png"))
plotout_road_kudu_fw
dev.off()
png(filename=paste0(img_save_dir, "PtMod3river_kudu.png"))
plotout_river_kudu_fw
dev.off()




















