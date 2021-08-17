# This file collates mean error results writen out from CRESS, MGCV and SVM models into tables for comparison
# set the directory locations in lines 4-9 then run whole file.

# base directory to save all output
base_dir_on_comp <- "C:/Users/christina/"
base_dir_sim <- base_dir_on_comp + "multinomial_simulations/"
# subdirectories to save csv output and generated images
file_save_dir <- paste0(base_dir_sim, "csv_files/")
img_save_dir <- paste0(base_dir_sim, "image_files/")


library(tidyverse)
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


cbbPalette <- c("#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")


simp_mgcv <- read.csv(paste0(file_save_dir,"mgcv_simp_dat2.csv"), sep=",", header=TRUE)
simp_svm <- read.csv(paste0(file_save_dir,"svm_simp_dat2.csv"), sep=",", header=TRUE)
simp_cress <- read.csv(paste0(file_save_dir,"cress_simp_dat2_new_cv.csv"), sep=",", header=TRUE)
island_mgcv <- read.csv(paste0(file_save_dir,"mgcv_island_dat2_new.csv"), sep=",", header=TRUE)
island_svm <- read.csv(paste0(file_save_dir,"svm_island_dat2_new.csv"), sep=",", header=TRUE)
island_cress <- read.csv(paste0(file_save_dir,"cress_island_dat2_new_cv.csv"), sep=",", header=TRUE)
comp_mgcv <- read.csv(paste0(file_save_dir,"mgcv_comp_dat2.csv"), sep=",", header=TRUE)
comp_svm <- read.csv(paste0(file_save_dir,"svm_comp_dat2.csv"), sep=",", header=TRUE)
comp_cress <- read.csv(paste0(file_save_dir,"cress_comp_dat2_new_cv.csv"), sep=",", header=TRUE)

rss_cr_si <- read.csv(paste0(img_save_dir, "rss_cress_simp.csv"))
rss_mg_si <- read.csv(paste0(img_save_dir, "rss_mgcv_simp.csv"))
rss_sv_si <- read.csv(paste0(img_save_dir, "rss_svm_simp.csv"))
rss_cr_co <- read.csv(paste0(img_save_dir, "rss_cress_comp.csv"))
rss_mg_co <- read.csv(paste0(img_save_dir, "rss_mgcv_comp.csv"))
rss_sv_co <- read.csv(paste0(img_save_dir, "rss_svm_comp.csv"))
rss_cr_is <- read.csv(paste0(img_save_dir, "rss_cress_island.csv"))
rss_mg_is <- read.csv(paste0(img_save_dir, "rss_mgcv_island.csv"))
rss_sv_is <- read.csv(paste0(img_save_dir, "rss_svm_island.csv"))


get_tab_rss <- function(out_store, dataset, model, time_mat){
  repz <- 100
  timecol_lw <- paste0("time_", model, "_", dataset, "_lw")
  timecol_hi <- paste0("time_", model, "_", dataset, "_hi")
  times100_lw <- time_mat[seq(1, 3*repz, by=3), timecol_lw] 
  times500_lw <- time_mat[seq(2, 3*repz, by=3), timecol_lw] 
  times5000_lw <- time_mat[seq(3, 3*repz, by=3), timecol_lw] 
  times100_hi <- time_mat[seq(1, 3*repz, by=3), timecol_hi] 
  times500_hi <- time_mat[seq(2, 3*repz, by=3), timecol_hi] 
  times5000_hi <- time_mat[seq(3, 3*repz, by=3), timecol_hi] 
  dat100_lw <- out_store[,"lw_1"]
  dat500_lw <- out_store[,"lw_2"]
  dat5000_lw <- out_store[,"lw_3"]
  dat100_hi <- out_store[,"hi_1"]
  dat500_hi <- out_store[,"hi2"]
  dat5000_hi <- out_store[,"hi3"]
  col1 <- c(rep("low", 3), rep("high", 3))
  col2 <- c(100,500,5000,100,500,5000)
  col5 <- round(c(mean(dat100_lw), mean(dat500_lw), mean(dat5000_lw), mean(dat100_hi), mean(dat500_hi), mean(dat5000_hi)), 1)
  col6 <- round(c(sd(dat100_lw), sd(dat500_lw), sd(dat5000_lw), sd(dat100_hi), sd(dat500_hi), sd(dat5000_hi)), 1)
  col7 <- round(c(mean(times100_lw), mean(times500_lw), mean(times5000_lw), mean(times100_hi), mean(times500_hi), mean(times5000_hi)), 1)
  col8 <- rep(dataset, 6)
  col9 <- rep(model, 6)
  col10 <- round(c(median(dat100_lw), median(dat500_lw), median(dat5000_lw), median(dat100_hi), median(dat500_hi), median(dat5000_hi)), 1)
  col11 <- round(c(IQR(dat100_lw), IQR(dat500_lw), IQR(dat5000_lw), IQR(dat100_hi), IQR(dat500_hi), IQR(dat5000_hi)), 1)
  
  
  colz <- c("noise_level", "no_points", "mean", "sd",  "time", "dataset", "model", "median", "IQR")
  
  tab_all <- as.data.frame(cbind(col1,col2,col5,col6,col7,col8,col9,col10,col11))
  colnames(tab_all) <- colz
  
  return(tab_all)
}

tab_out <- get_tab_rss(rss_mg_si, "simp", "mgcv", simp_mgcv)
tab_out <- rbind(tab_out, get_tab_rss(rss_sv_si, "simp", "svm", simp_svm))
tab_out <- rbind(tab_out, get_tab_rss(rss_cr_si, "simp", "cress", simp_cress))
tab_out <- rbind(tab_out, get_tab_rss(rss_mg_co, "comp", "mgcv", comp_mgcv))
tab_out <- rbind(tab_out, get_tab_rss(rss_sv_co, "comp", "svm", comp_svm))
tab_out <- rbind(tab_out, get_tab_rss(rss_cr_co, "comp", "cress", comp_cress))
tab_out <- rbind(tab_out, get_tab_rss(rss_mg_is, "island", "mgcv", island_mgcv))
tab_out <- rbind(tab_out, get_tab_rss(rss_sv_is, "island", "svm", island_svm))
tab_out <- rbind(tab_out, get_tab_rss(rss_cr_is, "island", "cress", island_cress))

write.csv(tab_out, paste0(img_save_dir, "rss_comparison_table_new.csv"), row.names=F)

get_tab_disagg_rss <- function(out_store, dataset, model, time_mat){
  repz <- 100
  timecol_lw <- paste0("time_", model, "_", dataset, "_lw")
  timecol_hi <- paste0("time_", model, "_", dataset, "_hi")
  times100_lw <- time_mat[seq(1, 3*repz, by=3), timecol_lw] 
  times500_lw <- time_mat[seq(2, 3*repz, by=3), timecol_lw] 
  times5000_lw <- time_mat[seq(3, 3*repz, by=3), timecol_lw] 
  times100_hi <- time_mat[seq(1, 3*repz, by=3), timecol_hi] 
  times500_hi <- time_mat[seq(2, 3*repz, by=3), timecol_hi] 
  times5000_hi <- time_mat[seq(3, 3*repz, by=3), timecol_hi] 
  dat100_lw <- out_store[,"lw_1"]
  dat500_lw <- out_store[,"lw_2"]
  dat5000_lw <- out_store[,"lw_3"]
  dat100_hi <- out_store[,"hi_1"]
  dat500_hi <- out_store[,"hi2"]
  dat5000_hi <- out_store[,"hi3"]
  rss_out <- c(dat100_lw, dat500_lw, dat5000_lw, dat100_hi, dat500_hi, dat5000_hi)
  tim_out <- c(times100_lw, times500_lw, times5000_lw, times100_hi, times500_hi, times5000_hi)
  colz <- c("rss", "time")
  
  tab_all <- as.data.frame(cbind(rss_out, tim_out))
  colnames(tab_all) <- colz
  tab_all$dataset <- rep(dataset, repz*6)
  tab_all$model <- rep(model, repz*6)
  tab_all$noise <- c(rep("lw", repz*3), rep("hi", repz*3))
  tab_all$noPoints <- c(rep(100, repz), rep(500, repz), rep(5000, repz), rep(100, repz), rep(500, repz), rep(5000, repz))
  
  return(tab_all)
}

tab_out_disagg <- get_tab_disagg_rss(rss_mg_si, "simp", "mgcv", simp_mgcv)
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_sv_si, "simp", "svm", simp_svm))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_cr_si, "simp", "cress", simp_cress))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_mg_co, "comp", "mgcv", comp_mgcv))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_sv_co, "comp", "svm", comp_svm))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_cr_co, "comp", "cress", comp_cress))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_mg_is, "island", "mgcv", island_mgcv))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_sv_is, "island", "svm", island_svm))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg_rss(rss_cr_is, "island", "cress", island_cress))

write.csv(tab_out_disagg, paste0(img_save_dir, "rss_disagg_comparison_table_new.csv"), row.names=F)

get_sig_vecs_rss <- function(tab_out_disagg, datset, thresh){
  
  npts <- sort(unique(tab_out_disagg$noPoints))
  nz_lev <- unique(tab_out_disagg$noise)
  mgcv_vec <- rep(" ", 18)
  svm_vec <- rep(" ", 18)
  cress_vec <- rep(" ", 18)
  col <- 1
  
  for (nz in 1:length(nz_lev)){
    noize = nz_lev[nz]
    filter_tab <- tab_out_disagg[tab_out_disagg$dataset==datset&tab_out_disagg$noise==noize,]
    for (np in 1:length(npts)){
      pts_tab <- filter_tab[filter_tab$noPoints==npts[np],]
      vec_rw <- (col - 1) * length(npts) * length(nz_lev) + (nz - 1) * length(npts) + np
      wtest <- wilcox.test(pts_tab[pts_tab$model=='mgcv', col], 
                           pts_tab[pts_tab$model=='svm', col], 
                           paired=TRUE)
      if (wtest$p.value < thresh){
        svm_vec[vec_rw] <- "*"
        mgcv_vec[(length(npts)*length(nz_lev)) + vec_rw] <- "*"
      } 
      wtest <- wilcox.test(pts_tab[pts_tab$model=='mgcv', col], 
                           pts_tab[pts_tab$model=='cress', col], 
                           paired=TRUE)
      if (wtest$p.value < thresh){
        cress_vec[vec_rw] <- "*"
        mgcv_vec[(length(npts)*length(nz_lev)*2) + vec_rw] <- "*"
      } 
      wtest <- wilcox.test(pts_tab[pts_tab$model=='svm', col], 
                           pts_tab[pts_tab$model=='cress', col], 
                           paired=TRUE)
      if (wtest$p.value < thresh){
        cress_vec[(length(npts)*length(nz_lev)) + vec_rw] <- "*"
        svm_vec[(length(npts)*length(nz_lev)*2) + vec_rw] <- "*"
      } 
    }
  }
  
  sig_vec_out <- cbind(mgcv_vec, svm_vec, cress_vec)
  colnames(sig_vec_out) <- c("Diff2mgcv", "Diff2svm", "Diff2cress")

  return(sig_vec_out)
}

sigtest_mat <- get_sig_vecs_rss(tab_out_disagg, "simp", 0.05)
sigtest_mat <- rbind(sigtest_mat, get_sig_vecs_rss(tab_out_disagg, "comp", 0.05))
sigtest_mat <- rbind(sigtest_mat, get_sig_vecs_rss(tab_out_disagg, "island", 0.05))


tab_out <- cbind(tab_out, sigtest_mat)

write.csv(tab_out, paste0(img_save_dir, "rss_comparison_table_new.csv"), row.names=F)


wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)



wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$rss[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)

filter_valz <- cbind(c(rep("simp", 6), rep("comp", 6), rep("island", 6)),rep(c(rep("lw", 3), rep("hi", 3)), 3),rep(rep(c(100, 500, 5000), 2), 3))

f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


scale1 <- scale_x_continuous(limits=c(0,2500), breaks=seq(0,2500, by=500)) 
scale2 <- scale_x_continuous(limits=c(0,2500), breaks=seq(0,2500, by=250)) 
scale3 <- scale_x_continuous(limits=c(0,2500), breaks=seq(0,2500, by=100)) 
scale4 <- scale_x_continuous(limits=c(0,2500), breaks=seq(0,2500, by=50)) 
scale5 <- scale_x_continuous(limits=c(0,2500), breaks=seq(0,2500, by=25)) 
scale6 <- scale_x_continuous(limits=c(0,2500), breaks=seq(0,2500, by=10)) 
theme_thesis <- theme(text=element_text(size=20), panel.grid.major.x=element_line(color='gray'), panel.background=element_rect(fill="white", color="black"),  panel.ontop=F)
xl <- xlab("Mean error")
yl <- ylab("Modelling approach")
statsum <- stat_summary(fun.data = f, geom="boxplot", outlier.shape=NA, fill="gray", alpha=1)
yscale <- scale_y_discrete(labels = c('SVM','MGCV','CRESS'))

for (ii in 1:nrow(filter_valz)){
  dst <- filter_valz[ii, 1]
  nz <- filter_valz[ii, 2]
  np <- filter_valz[ii, 3]
  filenameout <- paste0("errorBox", dst, nz, np,".png")
  qntm <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="mgcv") %>% select(rss), use.names=F, recursive=T), probs=0.9)
  qnts <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="cress") %>% select(rss), use.names=F, recursive=T), probs=0.9)
  qntc <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="svm") %>% select(rss), use.names=F, recursive=T), probs=0.9)
  maxval <- max(qntm, qnts, qntc)
  if (maxval > 1500){
    scaleout <- scale1
  } else if (maxval > 1000){
    scaleout <- scale2
  } else if (maxval > 500) {
    scaleout <- scale3
  } else if (maxval > 250) {
    scaleout <- scale4
  } else if (maxval > 100) {
    scaleout <- scale5
  } else {
    scaleout <- scale6
  }
  print(filenameout)
  print(maxval)
  print(max(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np) %>% select(rss)))
  plt <- ggplot(data=tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np)%>% mutate(model=factor(model, levels=c("svm", "mgcv", "cress"))), aes(rss, model)) +
    statsum + yscale + scaleout + theme_thesis + xl + yl + coord_cartesian(xlim=c(0,maxval))
  png(filename=paste0(img_save_dir, filenameout))
  print(plt)
  dev.off()
}

# same scale for comparisons

for (ii in 1:nrow(filter_valz)){
  dst <- filter_valz[ii, 1]
  nz <- filter_valz[ii, 2]
  np <- filter_valz[ii, 3]
  filenameout <- paste0("errorBox", dst, nz, np,"_samescale.png")
    plt <- ggplot(data=tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np)%>% mutate(model=factor(model, levels=c("svm", "mgcv", "cress"))), aes(rss, model)) +
    statsum + yscale + scale2 + theme_thesis + xl + yl + coord_cartesian(xlim=c(0,1750))
  png(filename=paste0(img_save_dir, filenameout))
  print(plt)
  dev.off()
}


#facet_grid(cols=vars(noPoints), scales = "free") +
