# This file fits the binomial any animal model and creates output images
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

cbbPalette_colours <- c(dark_blue="#004949",mid_green="#65aa70",mid_pink="#ff6db6",light_pink="#ffb6db",
                        dark_purple="#490092","#006ddb",mid_blue="#b66dff",light_blue="#afeeee",
                        dark_red="#920000",brown="#924900",orange="#db6d00",green="#24ff24",yellow="#ffff6d", 
                        black="#000000", light_grey="#DEDEDE")
cbb_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (cbbPalette_colours)
  
  cbbPalette_colours[cols]
}
cbb_palettes <- list(three_class = cbb_cols("dark_red", "mid_green", "light_blue"))
cbb_palettes <- list(six_class = cbb_cols("light_grey", "mid_green", "mid_pink", "light_pink", "dark_purple", "mid_blue", "dark_red"))
cbb_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- cbb_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}
scale_color_cbb <- function(palette = "three_class", discrete = TRUE, reverse = FALSE, ...) {
  pal <- cbb_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("cbb_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}
scale_fill_cbb <- function(palette = "three_class", discrete = TRUE, reverse = FALSE, ...) {
  pal <- cbb_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("cbb_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
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

sectionz_ndvi$response <- sectionz_ndvi$AnimalsPA
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

var_list48 <- c("dist_dam_clip", "dist_minor_clip")
if (!exists("PA_48_mod")) {
  PA_48_mod <- rez_for_var_list(var_list48, PA_mod, distMats)
}

coefz48 <- PA_48_mod$mod_2d$bestModel$coefficients
covmat48 <- summary(PA_48_mod$mod_2d$bestModel)$cov.scaled
rcoefs48 <- rmvnorm(1000, coefz48, sigma=covmat48)

if (!exists("pred_bootz48_all")){
  pred_bootz48_all <- predict.gamMRSea(object=PA_48_mod$mod_2d$bestModel, coeff=rcoefs48[1,])
  for (bt in 2:1000){
    pred_boot <- predict.gamMRSea(object=PA_48_mod$mod_2d$bestModel, coeff=rcoefs48[bt,])
    pred_bootz48_all <- cbind(pred_bootz48_all, pred_boot)
  }
  pred_ci48 <- t(apply(pred_bootz48_all, 1, quantile, probs=c(0.025, 0.975)))
}

pred_48 <- predict.gamMRSea(object=PA_48_mod$mod_2d$bestModel)
data_48 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_48, pred_ci48))
colnames(data_48) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

data_48 <- as.data.frame(cbind(data_48, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_48) <- c(colnames(data_48[1:5]), "UTMX", "UTMY", "distX", "distY")

if (!exists("pred_bootz48_all_link")){
  pred_bootz48_all_link <- predict.gamMRSea(object=PA_48_mod$mod_2d$bestModel, coeff=rcoefs48[1,], type="link")
  for (bt in 2:1000){
    pred_boot <- predict.gamMRSea(object=PA_48_mod$mod_2d$bestModel, coeff=rcoefs48[bt,], type="link")
    pred_bootz48_all <- cbind(pred_bootz48_all_link, pred_boot)
  }
  pred_ci48_link <- t(apply(pred_bootz48_all_link, 1, quantile, probs=c(0.025, 0.975)))
}

pred_48_link <- predict.gamMRSea(object=PA_48_mod$mod_2d$bestModel, type="link")
data_48_link <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_48_link, pred_ci48_link))
colnames(data_48_link) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

data_48_link <- as.data.frame(cbind(data_48_link, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_48_link) <- c(colnames(data_48_link[1:5]), "UTMX", "UTMY", "distX", "distY")

graph_pred48 <- ggplot(data_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred48 <- ggplot(data_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred48_lw <- ggplot(data_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n lower limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred48_lw <- ggplot(data_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred48_hi <- ggplot(data_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Predictions from model \n upper limit of confidence interval")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred48_hi <- ggplot(data_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "GraphModel2Pred.png"))
graph_pred48
dev.off()

png(filename=paste0(img_save_dir, "GraphModel2Lw.png"))
graph_pred48_lw
dev.off()

png(filename=paste0(img_save_dir, "GraphModel2Hi.png"))
graph_pred48_hi
dev.off()

# filter out any NAs
testcoefs48 <- PA_48_mod$mod_2d$bestModel$coefficients
cfmsk48 <- !is.na(testcoefs48)
tstcfs48 <- testcoefs48[cfmsk48]

sp_col48 <- PA_48_mod$mod_2d$bestModel$splineParams[[1]]

# get 2d columns
nam2d <- "LRF.g" 
strlenclnam <- str_length(nam2d)
coefnamsplt48 <- str_sub(names(tstcfs48),1,strlenclnam)
coefmask48 <- coefnamsplt48 == nam2d

# create radial gaussian bases
radii48 <- sp_col48$radii
radiiInd48 <- sp_col48$radiusIndices
aR48 <- sp_col48$knotPos
lrf48 <- LRF.g(radiiInd48, distMats$dataDist, radii48, aR48)

# combine coefmask and facts
coefz48 <- tstcfs48[coefmask48]
# get predicted values on link scale
predtm48 <- lrf48 %*% coefz48
# convert to response
# predtm48 <- PA_48_mod$mod_2d$bestModel$family$linkinv(predtm48)

if (!exists("pred_2d_ci48")) {
  bootcoefz48 <- rcoefs48[1,][cfmsk48][coefmask48]
  predboot48_2d <- lrf48 %*% bootcoefz48
  # predboot48_2d <- PA_48_mod$mod_2d$bestModel$family$linkinv(predboot48_2d)
  for (bt in 2:1000){
    bootcf <- rcoefs48[bt,][cfmsk48][coefmask48]
    predbt <- lrf48 %*% bootcf
    # predbt <- PA_48_mod$mod_2d$bestModel$family$linkinv(predbt)
    predboot48_2d <- cbind(predboot48_2d, predbt)
  }
  pred_2d_ci48 <- t(apply(predboot48_2d, 1, quantile, probs=c(0.025, 0.975)))
}

data_2d_48 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, predtm48, pred_2d_ci48))
colnames(data_2d_48) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

data_2d_48 <- as.data.frame(cbind(data_2d_48, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_2d_48) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")

graph_2d_48 <- ggplot(data_2d_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_2d_48 <- ggplot(data_2d_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_2d_48_lw <- ggplot(data_2d_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_2d_48_lw <- ggplot(data_2d_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_2d_48_hi <- ggplot(data_2d_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_2d_48_hi <- ggplot(data_2d_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

data_not2d_48 <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_48_link - predtm48, pred_ci48_link - pred_2d_ci48))
colnames(data_not2d_48) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI")

data_not2d_48 <- as.data.frame(cbind(data_not2d_48, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_not2d_48) <- c("Lat", "Long", "Prediction", "LowerCI", "UpperCI", "UTMX", "UTMY", "distX", "distY")


graph_not2d_48 <- ggplot(data_not2d_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_not2d_48 <- ggplot(data_not2d_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_not2d_48_lw <- ggplot(data_not2d_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_not2d_48_lw <- ggplot(data_not2d_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_not2d_48_hi <- ggplot(data_not2d_48, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from 2d smooth")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_not2d_48_hi <- ggplot(data_not2d_48, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))


graph_pred48_link <- ggplot(data_48_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred48_link <- ggplot(data_48_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Prediction), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred48_lw_link <- ggplot(data_48_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model lower limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred48_lw_link <- ggplot(data_48_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=LowerCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_pred48_hi_link <- ggplot(data_48_link, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("GFRC predictions from model upper limit of confidence interval on link scale")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_pred48_hi_link <- ggplot(data_48_link, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=UpperCI), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(-6,2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "GraphModel22d.png"))
  graph_2d_48
dev.off()
png(filename=paste0(img_save_dir, "GraphModel22dLw.png"))
  graph_2d_48_lw
dev.off()
png(filename=paste0(img_save_dir, "GraphModel22dHi.png"))
  graph_2d_48_hi
dev.off()
png(filename=paste0(img_save_dir, "GraphModel2NotPrd.png"))
  graph_not2d_48
dev.off()
png(filename=paste0(img_save_dir, "GraphModel2NotLw.png"))
  graph_not2d_48_lw
dev.off()
png(filename=paste0(img_save_dir, "GraphModel2NotHi.png"))
  graph_not2d_48_hi
dev.off()
png(filename=paste0(img_save_dir, "GraphModel2Lnk.png"))
graph_pred48_link
dev.off()
png(filename=paste0(img_save_dir, "GraphModel2LnkLw.png"))
graph_pred48_lw_link
dev.off()
png(filename=paste0(img_save_dir, "GraphModel2LnkHi.png"))
graph_pred48_hi_link
dev.off()

## 1d calculations

# create prediction vals for variables
elev_predvals <- seq(min(sectionz_ndvi$Elev), max(sectionz_ndvi$Elev), length.out=200)
ndvi_predvals <- seq(min(sectionz_ndvi$NDVI250), max(sectionz_ndvi$NDVI250), length.out=200)
dam_predvals <- seq(min(sectionz_ndvi$dist_dam_clip), max(sectionz_ndvi$dist_dam_clip), length.out=200)
minor_predvals <- seq(min(sectionz_ndvi$dist_minor_clip), max(sectionz_ndvi$dist_minor_clip), length.out=200)
fence_predvals <- seq(min(sectionz_ndvi$dist_fence_clip), max(sectionz_ndvi$dist_fence_clip), length.out=200)
track_predvals <- seq(min(sectionz_ndvi$dist_track_clip), max(sectionz_ndvi$dist_track_clip), length.out=200)
wat_predvals <- seq(min(sectionz_ndvi$dist_wat_clip), max(sectionz_ndvi$dist_wat_clip), length.out=200)


# normalise data
xpos_normal <- (sectionz_ndvi$x.pos - mean(sectionz_ndvi$x.pos)) / sd(sectionz_ndvi$x.pos)
ypos_normal <- (sectionz_ndvi$y.pos - mean(sectionz_ndvi$y.pos)) / sd(sectionz_ndvi$y.pos)
dam_normal <- (sectionz_ndvi$dist_dam_clip - mean(sectionz_ndvi$dist_dam_clip)) / sd(sectionz_ndvi$dist_dam_clip)
minor_normal <- (sectionz_ndvi$dist_minor_clip - mean(sectionz_ndvi$dist_minor_clip)) / sd(sectionz_ndvi$dist_minor_clip)

dist_med_xpos <- (xpos_normal - median(xpos_normal))^2
dist_med_ypos <- (ypos_normal - median(ypos_normal))^2
dist_med_dam <- (dam_normal - median(dam_normal))^2
dist_med_minor <- (minor_normal - median(minor_normal))^2

dist_med_48 <- dist_med_xpos + dist_med_ypos + dist_med_dam + dist_med_minor
index_med_48 <- which.min(dist_med_48)
median_vals_48 <- sectionz_ndvi[index_med_48,]
distMats_pred_48 <- makeDists(
  cbind(rep(median_vals_48$x.pos, 200), rep(median_vals_48$y.pos, 200)),
  na.omit(knotgrid)
)
rug_vals_48 <- as.data.frame(cbind(pred_48, sectionz_ndvi$dist_dam_clip, sectionz_ndvi$dist_minor_clip))
colnames(rug_vals_48) <- c("predictions", "dist_dam_clip", "dist_minor_clip")

rug_vals_48$dist_minor_km <- rug_vals_48$dist_minor_clip / 1000
rug_vals_48$dist_dam_km <- rug_vals_48$dist_dam_clip / 1000


predict_dam48 <- as.data.frame(cbind(rep(median_vals_48$x.pos, 200), rep(median_vals_48$y.pos, 200), dam_predvals,
                                     rep(median_vals_48$dist_minor_clip, 200)))
colnames(predict_dam48) <- c("x.pos", "y.pos", "dist_dam_clip", "dist_minor_clip")

pred_dam48 <- predict.gamMRSea(newdata=predict_dam48, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pred_48$dataDist)

pred_bootz48_dam <- predict.gamMRSea(newdata=predict_dam48, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pred_48$dataDist, coeff=rcoefs48[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_dam48, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pred_48$dataDist, coeff=rcoefs48[bt,])
  pred_bootz48_dam <- cbind(pred_bootz48_dam, pred_boot)
}
pred_ci_dam48 <- t(apply(pred_bootz48_dam, 1, quantile, probs=c(0.025, 0.975)))

pred_df_dam48 <- cbind(predict_dam48, pred_dam48, pred_ci_dam48)
colnames(pred_df_dam48) <- c(colnames(predict_dam48), "predictions", "LowerCI", "UpperCI")

pred_df_dam48$dist_dam_km <- pred_df_dam48$dist_dam_clip / 1000

predict_minor48 <- as.data.frame(cbind(rep(median_vals_48$x.pos, 200), rep(median_vals_48$y.pos, 200), rep(median_vals_48$dist_dam_clip, 200), minor_predvals))
colnames(predict_minor48) <- c("x.pos", "y.pos", "dist_dam_clip", "dist_minor_clip")

pred_minor48 <- predict.gamMRSea(newdata=predict_minor48, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pred_48$dataDist)

pred_bootz48_minor <- predict.gamMRSea(newdata=predict_minor48, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pred_48$dataDist, coeff=rcoefs48[1,])
for (bt in 2:1000){
  pred_boot <- predict.gamMRSea(newdata=predict_minor48, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pred_48$dataDist, coeff=rcoefs48[bt,])
  pred_bootz48_minor <- cbind(pred_bootz48_minor, pred_boot)
}
pred_ci_minor48 <- t(apply(pred_bootz48_minor, 1, quantile, probs=c(0.025, 0.975)))

pred_df_minor48 <- cbind(predict_minor48, pred_minor48, pred_ci_minor48)
colnames(pred_df_minor48) <- c(colnames(predict_minor48), "predictions", "LowerCI", "UpperCI")

pred_df_minor48$dist_minor_km <- pred_df_minor48$dist_minor_clip / 1000

pointz_vec <- rep(0, nrow(sectionz_ndvi))
pointz_vec[index_med_48] <- 0

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")

df_pointz <- as.data.frame(cbind(df_pointz, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(df_pointz) <- c("Lat", "Long", "Points", "UTMX", "UTMY", "distX", "distY")
df_pointz$PointsFac <- as.factor(df_pointz$Points)

max_fac_func <- function(vector_in){
  max_vector <- max(vector_in)
  string_out <- as.character(max_vector)
  return (string_out)
}

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
  xlim(670, 810) + ylim(-3170, -2970) +
  coord_fixed(ratio = 1) +
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max_fac_func, binwidth = c(3, 3)) +
  scale_fill_cbb(guide = guide_legend(title=""), labels=c("survey area", "median point"), palette = "six_class") +
  theme(legend.position="right", text = element_text(size=20))

plotout_minor48 <- ggplot(pred_df_minor48, aes(x=dist_minor_km, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor river (km)') + geom_ribbon(aes(x=dist_minor_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_48, sides='b', size=0.5, length=unit(0.02,'npc'))  + theme(text = element_text(size=20))

plotout_dam48 <- ggplot(pred_df_dam48, aes(x=dist_dam_km, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to dam (km)') + geom_ribbon(aes(x=dist_dam_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_48, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))

png(filename=paste0(img_save_dir, "PtMed.png"))
pt_med
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Minor.png"))
plotout_minor48
dev.off()
png(filename=paste0(img_save_dir, "PtMod2Dam.png"))
plotout_dam48
dev.off()

var_list_dm <- c("dist_dam_clip")

if (!exists("PA_mod_dam")) {
  PA_mod_dam <- rez_for_var_list(var_list_dm, PA_mod, distMats)
}

pred_d <- predict.gamMRSea(object=PA_mod_dam$mod_2d$bestModel)
data_d <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pred_48, pred_d, pred_48 - pred_d))
colnames(data_d) <- c("Lat", "Long", "PredOrig", "PredDam", "PredOrigMinusDam")

data_d <- as.data.frame(cbind(data_d, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(data_d) <- c("Lat", "Long", "PredOrig", "PredDam", "PredOrigMinusDam", "UTMX", "UTMY", "distX", "distY")

graph_predd <- ggplot(data_d, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") +
  # ggtitle(paste("GFRC predictions from model")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary_hex(aes(z=PredDam), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.15), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_predd <- ggplot(data_d, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=PredDam), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

graph_predm <- ggplot(data_d, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") +
  # ggtitle(paste("GFRC predictions from model")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary_hex(aes(z=PredOrigMinusDam), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.15), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_predm <- ggplot(data_d, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(670, 810) + ylim(-3170, -2970) + 
  # ggtitle(paste("Springbok predictions from multinomial model")) + 
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=PredOrigMinusDam), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.2), oob=squish, name="probability") +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "PredD.png"))
graph_predd
dev.off()
png(filename=paste0(img_save_dir, "PredM.png"))
graph_predm
dev.off()

pointz_vec <- rep(0, nrow(sectionz_ndvi))
pt1 <- c(18, -26.95)
dist_pt1 <- (sectionz_ndvi$Long - pt1[1])^2 + (sectionz_ndvi$Lat - pt1[2])^2
index_pt1 <- which.min(dist_pt1)
pt2 <- c(18.05, -27.4)
dist_pt2 <- (sectionz_ndvi$Long - pt2[1])^2 + (sectionz_ndvi$Lat - pt2[2])^2
index_pt2 <- which.min(dist_pt2)
pt3 <- c(17.99, -27.58)
dist_pt3 <- (sectionz_ndvi$Long - pt3[1])^2 + (sectionz_ndvi$Lat - pt3[2])^2
index_pt3 <- which.min(dist_pt3)
pt4 <- c(17.82, -27.625)
dist_pt4 <- (sectionz_ndvi$Long - pt4[1])^2 + (sectionz_ndvi$Lat - pt4[2])^2
index_pt4 <- which.min(dist_pt4)
pt5 <- c(17.7, -27.55)
dist_pt5 <- (sectionz_ndvi$Long - pt5[1])^2 + (sectionz_ndvi$Lat - pt5[2])^2
index_pt5 <- which.min(dist_pt5)
pt6 <- c(17.6, -27.625)
dist_pt6 <- (sectionz_ndvi$Long - pt6[1])^2 + (sectionz_ndvi$Lat - pt6[2])^2
index_pt6 <- which.min(dist_pt6)

pointz_vec[index_pt1] <- 1
pointz_vec[index_pt2] <- 2
pointz_vec[index_pt3] <- 3
pointz_vec[index_pt4] <- 4
pointz_vec[index_pt5] <- 5
pointz_vec[index_pt6] <- 6

df_pointz <- as.data.frame(cbind(sectionz_ndvi$Lat, sectionz_ndvi$Long, pointz_vec))
colnames(df_pointz) <- c("Lat", "Long", "Points")

df_pointz <- as.data.frame(cbind(df_pointz, (sectionz_ndvi$x.pos / 1000), (sectionz_ndvi$y.pos / 1000), (sectionz_ndvi$x.pos / 1000)-680, (sectionz_ndvi$y.pos / 1000)+3170))
colnames(df_pointz) <- c("Lat", "Long", "Points", "UTMX", "UTMY", "distX", "distY")

pt_locs <- ggplot(df_pointz, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Locations of points used to predict")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "Set1", limits=c(0.1,6), name="") +
  theme(legend.position="right", text = element_text(size=20))
pt_locs <- ggplot(df_pointz, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  xlim(670, 810) + ylim(-3170, -2970) +
  coord_fixed(ratio = 1) +
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=Points), fun=max_fac_func, binwidth = c(3, 3)) +
  scale_fill_cbb(guide = guide_legend(title=""), labels=c("survey area", "point a", "point b", "point c", "point d", "point e", "point f"), palette = "six_class") +
  theme(legend.position="right", text = element_text(size=20))

df_pt_ex <- data_d
df_pt_ex$pointz <- df_pt_ex$PredOrig + df_pointz$Points / 20

graph_pt_ex <- ggplot(df_pt_ex, aes(Long, Lat)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(16.8, 18.1) + ylim(-28.7, -26.8) + xlab("") + ylab("") + 
  # ggtitle(paste("Points over explore")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=pointz), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.15), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))
graph_pt_ex <- ggplot(df_pt_ex, aes(UTMX, UTMY)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlab("UTM coordinate (km)") + ylab("UTM coordinate (km)") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_hex(aes(z=pointz), fun=max, binwidth = c(2, 2)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, limits=c(0,0.15), name = "", oob=squish) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "PtLocs.png"))
pt_locs
dev.off()
png(filename=paste0(img_save_dir, "PtEx.png"))
graph_pt_ex
dev.off()

minor_preds_at_point <- function(point){
  dist_point <- (sectionz_ndvi$Long - point[1])^2 + (sectionz_ndvi$Lat - point[2])^2
  index_point <- which.min(dist_point)
  median_vals_pt <- sectionz_ndvi[index_point,]
  distMats_pt <- makeDists(
    cbind(rep(median_vals_pt$x.pos, 200), rep(median_vals_pt$y.pos, 200)),
    na.omit(knotgrid)
  )
  
  predict_minorpt <- as.data.frame(cbind(rep(median_vals_pt$x.pos, 200), rep(median_vals_pt$y.pos, 200), rep(median_vals_pt$dist_dam_clip, 200), minor_predvals))
  colnames(predict_minorpt) <- c("x.pos", "y.pos", "dist_dam_clip", "dist_minor_clip")
  
  pred_minorpt <- predict.gamMRSea(newdata=predict_minorpt, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pt$dataDist)
  
  pred_bootzpt_minor <- predict.gamMRSea(newdata=predict_minorpt, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pt$dataDist, coeff=rcoefs48[1,])
  for (bt in 2:1000){
    pred_boot <- predict.gamMRSea(newdata=predict_minorpt, object=PA_48_mod$mod_2d$bestModel, g2k=distMats_pt$dataDist, coeff=rcoefs48[bt,])
    pred_bootzpt_minor <- cbind(pred_bootzpt_minor, pred_boot)
  }
  pred_ci_minorpt <- t(apply(pred_bootzpt_minor, 1, quantile, probs=c(0.025, 0.975)))
  
  pred_df_minorpt <- cbind(predict_minorpt, pred_minorpt, pred_ci_minorpt)
  colnames(pred_df_minorpt) <- c(colnames(predict_minorpt), "predictions", "LowerCI", "UpperCI")
  pred_df_minorpt$dist_minor_km <- pred_df_minorpt$dist_minor_clip / 1000
  
  plotout_minorpt <- ggplot(pred_df_minorpt, aes(x=dist_minor_km, y=predictions)) + theme_bw() + ylab('predicted probability') + xlab('Distance to minor river (km)') + geom_ribbon(aes(x=dist_minor_km, ymin=LowerCI, ymax=UpperCI), fill='#d3d3d3') + geom_path(color='red') + geom_rug(data=rug_vals_48, sides='b', size=0.5, length=unit(0.02,'npc')) + theme(text = element_text(size=20))
  
  return(plotout_minorpt)
  
}

plot_pt1 <- minor_preds_at_point(c(18, -26.95))
plot_pt2 <- minor_preds_at_point(c(17.9, -27.4))
plot_pt3 <- minor_preds_at_point(c(17.95, -27.625)) 
plot_pt4 <- minor_preds_at_point(c(17.75, -27.7))
plot_pt5 <- minor_preds_at_point(c(17.2, -27.75))
plot_pt6 <- minor_preds_at_point(c(17.7, -28.5))

png(filename=paste0(img_save_dir, "PlotPt1.png"))
plot_pt1
dev.off()
png(filename=paste0(img_save_dir, "PlotPt2.png"))
plot_pt2
dev.off()
png(filename=paste0(img_save_dir, "PlotPt3.png"))
plot_pt3
dev.off()
png(filename=paste0(img_save_dir, "PlotPt4.png"))
plot_pt4
dev.off()
png(filename=paste0(img_save_dir, "PlotPt5.png"))
plot_pt5
dev.off()
png(filename=paste0(img_save_dir, "PlotPt6.png"))
plot_pt6
dev.off()

cv48 <- cv.gamMRSea(sectionz_ndvi, PA_48_mod$mod_2d$bestModel, K=5)
# 0.01523283 0.01523037

elev_median <- median(sectionz_ndvi$Elev)
ndvi_median <- median(sectionz_ndvi$NDVI250)
dam_median <- median(sectionz_ndvi$dist_dam_clip)
minor_median <- median(sectionz_ndvi$dist_minor_clip)
fence_median <- median(sectionz_ndvi$dist_fence_clip)
track_median <- median(sectionz_ndvi$dist_track_clip)
wat_median <- median(sectionz_ndvi$dist_wat_clip)
xpos_median <- median(sectionz_ndvi$x.pos)
ypos_median <- median(sectionz_ndvi$y.pos)

# calculate distance to medians

dist_med_elev <- (sectionz_ndvi$Elev - elev_median)^2
dist_med_ndvi <- (sectionz_ndvi$NDVI250 - ndvi_median)^2
dist_med_dam <- (sectionz_ndvi$dist_dam_clip - dam_median)^2
dist_med_minor <- (sectionz_ndvi$dist_minor_clip - minor_median)^2
dist_med_fence <- (sectionz_ndvi$dist_fence_clip - fence_median)^2
dist_med_track <- (sectionz_ndvi$dist_track_clip - track_median)^2
dist_med_wat <- (sectionz_ndvi$dist_wat_clip - wat_median)^2
dist_med_xpos <- (sectionz_ndvi$x.pos - xpos_median)^2
dist_med_ypos <- (sectionz_ndvi$y.pos - ypos_median)^2

