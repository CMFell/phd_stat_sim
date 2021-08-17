# This file collates accuracy results writen out from CRESS, MGCV and SVM models into tables for comaprison
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

sm_valz <- read.csv(paste0(file_save_dir,"mgcv_simp_valz2.csv"), sep=",", header=TRUE)
ss_valz <- read.csv(paste0(file_save_dir,"svm_simp_valz2.csv"), sep=",", header=TRUE)
sc_valz <- read.csv(paste0(file_save_dir,"cress_simp_valz2_new_cv.csv"), sep=",", header=TRUE)
cm_valz <- read.csv(paste0(file_save_dir,"mgcv_comp_valz2.csv"), sep=",", header=TRUE)
cs_valz <- read.csv(paste0(file_save_dir,"svm_comp_valz2.csv"), sep=",", header=TRUE)
cc_valz <- read.csv(paste0(file_save_dir,"cress_comp_valz2_new_cv.csv"), sep=",", header=TRUE)
im_valz <- read.csv(paste0(file_save_dir,"mgcv_island_valz2_new.csv"), sep=",", header=TRUE)
is_valz <- read.csv(paste0(file_save_dir,"svm_island_valz2_new.csv"), sep=",", header=TRUE)
ic_valz <- read.csv(paste0(file_save_dir,"cress_island_valz2_new_cv.csv"), sep=",", header=TRUE)

get_tab_data <- function(out_store, dataset, model, valz){
  repz <- as.numeric(as.character(valz[1,1]))
  dat50 <- out_store[seq(from=1, by=3, length.out=repz),]
  dat500 <- out_store[seq(from=2, by=3, length.out=repz),]
  dat5000 <- out_store[seq(from=3, by=3, length.out=repz),]
  acc_tru_naml <- paste0("acc_truth_", model, "_all_", dataset, "_lw")
  acc_tru_namh <- paste0("acc_truth_", model, "_all_", dataset, "_hi")
  acc_dat_naml <- paste0("acc_data_", model, "_all_", dataset, "_lw")
  acc_dat_namh <- paste0("acc_data_", model, "_all_", dataset, "_hi")
  acc_fit_naml <- paste0("acc_truth_", model, "_fit_", dataset, "_lw")
  acc_fit_namh <- paste0("acc_truth_", model, "_fit_", dataset, "_hi")
  acc_fitdat_naml <- paste0("acc_data_", model, "_fit_", dataset, "_lw")
  acc_fitdat_namh <- paste0("acc_data_", model, "_fit_", dataset, "_hi")
  tim_naml <- paste0("time_", model, "_", dataset, "_lw")
  tim_namh <- paste0("time_", model, "_", dataset, "_hi")
  col1 <- rep(c(rep("low", 3), rep("high", 3)), 4)
  col2 <- rep(c(100,500,5000,100,500,5000), 4)
  col3 <- c(rep("all", 12), rep("fit", 12))
  col4 <- rep(c(rep("truth", 6), rep("data", 6)), 2)
  col5 <- round(c(mean(dat50[,acc_tru_naml]), mean(dat500[,acc_tru_naml]), mean(dat5000[,acc_tru_naml]), mean(dat50[,acc_tru_namh]), mean(dat500[,acc_tru_namh]), mean(dat5000[,acc_tru_namh]),mean(dat50[,acc_dat_naml]), mean(dat500[,acc_dat_naml]), mean(dat5000[,acc_dat_naml]), mean(dat50[,acc_dat_namh]), mean(dat500[,acc_dat_namh]), mean(dat5000[,acc_dat_namh]), mean(dat50[,acc_fit_naml]), mean(dat500[,acc_fit_naml]), mean(dat5000[,acc_fit_naml]), mean(dat50[,acc_fit_namh]), mean(dat500[,acc_fit_namh]), mean(dat5000[,acc_fit_namh]), mean(dat50[,acc_fitdat_naml]), mean(dat500[,acc_fitdat_naml]), mean(dat5000[,acc_fitdat_naml]), mean(dat50[,acc_fitdat_namh]), mean(dat500[,acc_fitdat_namh]), mean(dat5000[,acc_fitdat_namh])), 3)
  col6 <- round(c(sd(dat50[,acc_tru_naml]), sd(dat500[,acc_tru_naml]), sd(dat5000[,acc_tru_naml]), sd(dat50[,acc_tru_namh]), sd(dat500[,acc_tru_namh]), sd(dat5000[,acc_tru_namh]), sd(dat50[,acc_dat_naml]), sd(dat500[,acc_dat_naml]), sd(dat5000[,acc_dat_naml]), sd(dat50[,acc_dat_namh]), sd(dat500[,acc_dat_namh]), sd(dat5000[,acc_dat_namh]), sd(dat50[,acc_fit_naml]), sd(dat500[,acc_fit_naml]), sd(dat5000[,acc_fit_naml]), sd(dat50[,acc_fit_namh]), sd(dat500[,acc_fit_namh]), sd(dat5000[,acc_fit_namh]), sd(dat50[,acc_fitdat_naml]), sd(dat500[,acc_fitdat_naml]), sd(dat5000[,acc_fitdat_naml]), sd(dat50[,acc_fitdat_namh]), sd(dat500[,acc_fitdat_namh]), sd(dat5000[,acc_fitdat_namh])), 3)
  col7 <- rep(round(c(mean(dat50[,tim_naml]), mean(dat500[,tim_naml]), mean(dat5000[,tim_naml]), mean(dat50[,tim_namh]), mean(dat500[,tim_namh]), mean(dat5000[,tim_namh])), 1),4)
  col8 <- rep(dataset, 24)
  col9 <- rep(model, 24)
  col10 <- round(c(median(dat50[,acc_tru_naml]), median(dat500[,acc_tru_naml]), median(dat5000[,acc_tru_naml]), median(dat50[,acc_tru_namh]), median(dat500[,acc_tru_namh]), median(dat5000[,acc_tru_namh]),median(dat50[,acc_dat_naml]), median(dat500[,acc_dat_naml]), median(dat5000[,acc_dat_naml]), median(dat50[,acc_dat_namh]), median(dat500[,acc_dat_namh]), median(dat5000[,acc_dat_namh]), median(dat50[,acc_fit_naml]), median(dat500[,acc_fit_naml]), median(dat5000[,acc_fit_naml]), median(dat50[,acc_fit_namh]), median(dat500[,acc_fit_namh]), median(dat5000[,acc_fit_namh]), median(dat50[,acc_fitdat_naml]), median(dat500[,acc_fitdat_naml]), median(dat5000[,acc_fitdat_naml]), median(dat50[,acc_fitdat_namh]), median(dat500[,acc_fitdat_namh]), median(dat5000[,acc_fitdat_namh])), 3)
  col11 <- round(c(IQR(dat50[,acc_tru_naml]), IQR(dat500[,acc_tru_naml]), IQR(dat5000[,acc_tru_naml]), IQR(dat50[,acc_tru_namh]), IQR(dat500[,acc_tru_namh]), IQR(dat5000[,acc_tru_namh]),IQR(dat50[,acc_dat_naml]), IQR(dat500[,acc_dat_naml]), IQR(dat5000[,acc_dat_naml]), IQR(dat50[,acc_dat_namh]), IQR(dat500[,acc_dat_namh]), IQR(dat5000[,acc_dat_namh]), IQR(dat50[,acc_fit_naml]), IQR(dat500[,acc_fit_naml]), IQR(dat5000[,acc_fit_naml]), IQR(dat50[,acc_fit_namh]), IQR(dat500[,acc_fit_namh]), IQR(dat5000[,acc_fit_namh]), IQR(dat50[,acc_fitdat_naml]), IQR(dat500[,acc_fitdat_naml]), IQR(dat5000[,acc_fitdat_naml]), IQR(dat50[,acc_fitdat_namh]), IQR(dat500[,acc_fitdat_namh]), IQR(dat5000[,acc_fitdat_namh])), 3)
  
  colz <- c("noise_level", "no_points", "all_fitted", "data_truth", "mean", "sd",  "time", "dataset", "model", "median", "IQR")
  
  tab_all <- as.data.frame(cbind(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11))
  colnames(tab_all) <- colz
  
  return(tab_all)
}

tab_out <- get_tab_data(simp_mgcv, "simp", "mgcv", sm_valz)
tab_out <- rbind(tab_out, get_tab_data(simp_svm, "simp", "svm", ss_valz))
tab_out <- rbind(tab_out, get_tab_data(simp_cress, "simp", "cress", sc_valz))
tab_out <- rbind(tab_out, get_tab_data(comp_mgcv, "comp", "mgcv", cm_valz))
tab_out <- rbind(tab_out, get_tab_data(comp_svm, "comp", "svm", cs_valz))
tab_out <- rbind(tab_out, get_tab_data(comp_cress, "comp", "cress", cc_valz))
tab_out <- rbind(tab_out, get_tab_data(island_mgcv, "island", "mgcv", im_valz))
tab_out <- rbind(tab_out, get_tab_data(island_svm, "island", "svm", is_valz))
tab_out <- rbind(tab_out, get_tab_data(island_cress, "island", "cress", ic_valz))

write.csv(tab_out, paste0(img_save_dir, "simple_comparison_table_new.csv"), row.names=F)



get_tab_disagg <- function(out_store, dataset, model, valz){
  repz <- as.numeric(as.character(valz[1,1]))
  dat50 <- out_store[seq(from=1, by=3, length.out=repz),]
  dat500 <- out_store[seq(from=2, by=3, length.out=repz),]
  dat5000 <- out_store[seq(from=3, by=3, length.out=repz),]
  acc_tru_naml <- paste0("acc_truth_", model, "_all_", dataset, "_lw")
  acc_tru_namh <- paste0("acc_truth_", model, "_all_", dataset, "_hi")
  acc_dat_naml <- paste0("acc_data_", model, "_all_", dataset, "_lw")
  acc_dat_namh <- paste0("acc_data_", model, "_all_", dataset, "_hi")
  acc_fit_naml <- paste0("acc_truth_", model, "_fit_", dataset, "_lw")
  acc_fit_namh <- paste0("acc_truth_", model, "_fit_", dataset, "_hi")
  acc_fitdat_naml <- paste0("acc_data_", model, "_fit_", dataset, "_lw")
  acc_fitdat_namh <- paste0("acc_data_", model, "_fit_", dataset, "_hi")
  tim_naml <- paste0("time_", model, "_", dataset, "_lw")
  tim_namh <- paste0("time_", model, "_", dataset, "_hi")
  acc_tru_all <- c(dat50[,acc_tru_naml], dat500[,acc_tru_naml], dat5000[,acc_tru_naml], dat50[,acc_tru_namh], dat500[,acc_tru_namh], dat5000[,acc_tru_namh])
  acc_dat_all <- c(dat50[,acc_dat_naml], dat500[,acc_dat_naml], dat5000[,acc_dat_naml], dat50[,acc_dat_namh], dat500[,acc_dat_namh], dat5000[,acc_dat_namh])
  acc_tru_fit <- c(dat50[,acc_fit_naml], dat500[,acc_fit_naml], dat5000[,acc_fit_naml], dat50[,acc_fit_namh], dat500[,acc_fit_namh], dat5000[,acc_fit_namh])
  acc_dat_fit <- c(dat50[,acc_fitdat_naml], dat500[,acc_fitdat_naml], dat5000[,acc_fitdat_naml], dat50[,acc_fitdat_namh], dat500[,acc_fitdat_namh], dat5000[,acc_fitdat_namh])
  tim_out <- c(dat50[,tim_naml],dat500[,tim_naml], dat5000[,tim_naml], dat50[,tim_namh], dat500[,tim_namh], dat5000[,tim_namh])
  
  colz <- c("accuracy_to_data_all", "accuracy_to_truth_all", "accuracy_to_data_fitted", "accuracy_to_truth_fitted", "time")
  
  tab_all <- as.data.frame(cbind(acc_dat_all, acc_tru_all, acc_dat_fit, acc_tru_fit, tim_out))
  colnames(tab_all) <- colz
  tab_all$dataset <- rep(dataset, length(acc_tru_all))
  tab_all$model <- rep(model, length(acc_tru_all))
  tab_all$noise <- c(rep("lw", repz*3), rep("hi", repz*3))
  tab_all$noPoints <- c(rep(100, repz), rep(500, repz), rep(5000, repz), rep(100, repz), rep(500, repz), rep(5000, repz))
  
  return(tab_all)
}


tab_out_disagg <- get_tab_disagg(simp_mgcv, "simp", "mgcv", sm_valz)
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(simp_svm, "simp", "svm", ss_valz))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(simp_cress, "simp", "cress", sc_valz))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(comp_mgcv, "comp", "mgcv", cm_valz))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(comp_svm, "comp", "svm", cs_valz))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(comp_cress, "comp", "cress", cc_valz))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(island_mgcv, "island", "mgcv", im_valz))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(island_svm, "island", "svm", is_valz))
tab_out_disagg <- rbind(tab_out_disagg, get_tab_disagg(island_cress, "island", "cress", ic_valz))

write.csv(tab_out_disagg, paste0(img_save_dir, "disagg_comparison_table_new.csv"), row.names=F)

get_sig_vecs <- function(tab_out_disagg, datset, thresh){
  
  npts <- sort(unique(tab_out_disagg$noPoints))
  nz_lev <- unique(tab_out_disagg$noise)
  mgcv_vec <- rep(" ", 72)
  svm_vec <- rep(" ", 72)
  cress_vec <- rep(" ", 72)
  for (col in 1:4) {
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
          mgcv_vec[(length(npts)*length(nz_lev)*4) + vec_rw] <- "*"
        } 
        wtest <- wilcox.test(pts_tab[pts_tab$model=='mgcv', col], 
                             pts_tab[pts_tab$model=='cress', col], 
                             paired=TRUE)
        if (wtest$p.value < thresh){
          cress_vec[vec_rw] <- "*"
          mgcv_vec[(length(npts)*length(nz_lev)*4*2) + vec_rw] <- "*"
        } 
        wtest <- wilcox.test(pts_tab[pts_tab$model=='svm', col], 
                             pts_tab[pts_tab$model=='cress', col], 
                             paired=TRUE)
        if (wtest$p.value < thresh){
          cress_vec[(length(npts)*length(nz_lev)*4) + vec_rw] <- "*"
          svm_vec[(length(npts)*length(nz_lev)*4*2) + vec_rw] <- "*"
        } 
      }
    }
  }
  sig_vec_out <- cbind(mgcv_vec, svm_vec, cress_vec)
  colnames(sig_vec_out) <- c("Diff2mgcv", "Diff2svm", "Diff2cress")

  return(sig_vec_out)
}

sigtest_mat <- get_sig_vecs(tab_out_disagg, "simp", 0.05)
sigtest_mat <- rbind(sigtest_mat, get_sig_vecs(tab_out_disagg, "comp", 0.05))
sigtest_mat <- rbind(sigtest_mat, get_sig_vecs(tab_out_disagg, "island", 0.05))


tab_out <- cbind(tab_out, sigtest_mat)

write.csv(tab_out, paste0(img_save_dir, "simple_comparison_table_new.csv"), row.names=F)


wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='simp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='simp'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='comp'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='comp'], 
            paired=TRUE)



wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='lw'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)


wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='mgcv'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)
wilcox.test(tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='cress'&tab_out_disagg$dataset=='island'], 
            tab_out_disagg$accuracy_to_truth_all[tab_out_disagg$noise=='hi'&tab_out_disagg$model=='svm'&tab_out_disagg$dataset=='island'], 
            paired=TRUE)

filter_valz <- cbind(c(rep("simp", 6), rep("comp", 6), rep("island", 6)),rep(c(rep("lw", 3), rep("hi", 3)), 3),rep(rep(c(100, 500, 5000), 2), 3))

f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


scale1 <- scale_x_continuous(limits=c(0,1), breaks=seq(0,1, by=0.1)) 
scale2 <- scale_x_continuous(limits=c(0,1), breaks=seq(0,1, by=0.05)) 
scale3 <- scale_x_continuous(limits=c(0,1), breaks=seq(0,1, by=0.025)) 
scale4 <- scale_x_continuous(limits=c(0,1), breaks=seq(0,1, by=0.01)) 
scale5 <- scale_x_continuous(limits=c(0,1), breaks=seq(0,1, by=0.005)) 

theme_thesis <- theme(text=element_text(size=20), panel.grid.major.x=element_line(color='gray'), panel.background=element_rect(fill="white", color="black"),  panel.ontop=F)
xl <- xlab("Accuracy")
yl <- ylab("Modelling approach")
statsum <- stat_summary(fun.data = f, geom="boxplot", outlier.shape=NA, fill="gray", alpha=1)
yscale <- scale_y_discrete(labels = c('SVM','MGCV','CRESS'))

for (ii in 1:nrow(filter_valz)){
  dst <- filter_valz[ii, 1]
  nz <- filter_valz[ii, 2]
  np <- filter_valz[ii, 3]
  filenameout <- paste0("accuracyBox", dst, nz, np,".png")
  qntm <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="mgcv") %>% select(accuracy_to_truth_all), use.names=F, recursive=T), probs=0.9)
  qnts <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="cress") %>% select(accuracy_to_truth_all), use.names=F, recursive=T), probs=0.9)
  qntc <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="svm") %>% select(accuracy_to_truth_all), use.names=F, recursive=T), probs=0.9)
  maxval <- max(qntm, qnts, qntc)
  qntm <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="mgcv") %>% select(accuracy_to_truth_all), use.names=F, recursive=T), probs=0.1)
  qnts <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="cress") %>% select(accuracy_to_truth_all), use.names=F, recursive=T), probs=0.1)
  qntc <- quantile(c(tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np, model=="svm") %>% select(accuracy_to_truth_all), use.names=F, recursive=T), probs=0.1)
  minval <- min(qntm, qnts, qntc)
  diff <- maxval - minval
  if (diff > 0.4){
    scaleout <- scale1
  } else if (diff > 0.2){
    scaleout <- scale2
  } else if (diff > 0.1) {
    scaleout <- scale3
  } else if (diff > 0.05) {
    scaleout <- scale4
  } else {
    scaleout <- scale5
  }
  print(filenameout)
  print(maxval)
  print(minval)
  plt <- ggplot(data=tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np) %>% mutate(model=factor(model, levels=c("svm", "mgcv", "cress"))), aes(accuracy_to_truth_all, model)) +
    statsum + yscale + scaleout + theme_thesis + xl + yl + coord_cartesian(xlim=c(minval,maxval))
  png(filename=paste0(img_save_dir, filenameout))
  print(plt)
  dev.off()
}

# Plot all on same scale

for (ii in 1:nrow(filter_valz)){
  dst <- filter_valz[ii, 1]
  nz <- filter_valz[ii, 2]
  np <- filter_valz[ii, 3]
  filenameout <- paste0("accuracyBox", dst, nz, np,"_samescale.png")
    print(filenameout)
  plt <- ggplot(data=tab_out_disagg %>% filter(dataset == dst, noise == nz, noPoints==np)%>% mutate(model=factor(model, levels=c("svm", "mgcv", "cress"))), aes(accuracy_to_truth_all, model)) +
    statsum + yscale + scale1 + theme_thesis + xl + yl 
  png(filename=paste0(img_save_dir, filenameout))
  print(plt)
  dev.off()
}
