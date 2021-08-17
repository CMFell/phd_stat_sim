# This file collates results writen out from CRESS, MGCV and SVM models on the simple data
# It creates images included in thesis for example showing average errors across area
# set the directory locations in lines 5-13 then run whole file.

# base directory to save all output
base_dir_on_comp <- "C:/Users/christina/"
base_dir_sim <- base_dir_on_comp + "multinomial_simulations/"
# subdirectories to save csv output and generated images
file_save_dir <- paste0(base_dir_sim, "csv_files/")
img_save_dir <- paste0(base_dir_sim, "image_files/")
# code currently in multinomial branch of MRSea so need to download package from github to base dir on comp 
# then switch to multinomial branch
devtools::load_all(paste0(base_dir_on_comp, "MRsea"))

library(tidyverse)
library(viridis)
library(hexbin)
library(car)
library(mvtnorm)
library(patchwork)

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output <- sub("1\\*", "", output) # Removes 1*10 but leaves anything else raised to 10 
  parse(text=output)
}



cbbPalette <- c("#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")

file_in_dir <- file_save_dir

grid <- read.csv(paste0(file_save_dir, "location_grid.csv"))
low_s2n_simp <- read.csv(paste0(file_save_dir, "low_s2n_simp.csv"))
high_s2n_simp <- read.csv(paste0(file_save_dir, "high_s2n_simp.csv"))

### CRESS

trup_col_namll <- paste0("truep_simp_lw_1")
trup_col_naml2 <- paste0("truep_simp_lw_2")
trup_col_naml3 <- paste0("truep_simp_lw_3")
trup_col_namhl <- paste0("truep_simp_hi_1")
trup_col_namh2 <- paste0("truep_simp_hi_2")
trup_col_namh3 <- paste0("truep_simp_hi_3")
predp_col_namll <- paste0("predp_cress_simp_lw_1")
predp_col_naml2 <- paste0("predp_cress_simp_lw_2")
predp_col_naml3 <- paste0("predp_cress_simp_lw_3")
predp_col_namhl <- paste0("predp_cress_simp_hi_1")
predp_col_namh2 <- paste0("predp_cress_simp_hi_2")
predp_col_namh3 <- paste0("predp_cress_simp_hi_3")

output_matrix <- read.csv(paste0(file_in_dir, "cress_simp_mat2_new_cv.csv"), sep=",", header=TRUE)

output_cr_si_lw <- as.data.frame(matrix(rep(0, 10000*12), ncol=12))
colnames(output_cr_si_lw) <- c("pA_lw_rez_1", "pA_lw_rez_2", "pA_lw_rez_3",
                               "pB_lw_rez_1", "pB_lw_rez_2", "pB_lw_rez_3",
                               "pC_lw_rez_1", "pC_lw_rez_2", "pC_lw_rez_3",
                               "cl_lw_rez_1", "cl_lw_rez_2", "cl_lw_rez_3")

output_cr_si_hi <- as.data.frame(matrix(rep(0, 10000*12), ncol=12))
colnames(output_cr_si_hi) <- c("pA_hi_rez_1", "pA_hi_rez_2", "pA_hi_rez_3",
                               "pB_hi_rez_1", "pB_hi_rez_2", "pB_hi_rez_3",
                               "pC_hi_rez_1", "pC_hi_rez_2", "pC_hi_rez_3",
                               "cl_hi_rez_1", "cl_hi_rez_2", "cl_hi_rez_3")

repz <- 100

rss_cr_si <- as.data.frame(matrix(rep(0, repz*7), ncol=7))
colnames(rss_cr_si) <- c("sim", "lw_1", "lw_2", "lw_3", "hi_1", "hi2", "hi3")
rss_cr_si$sim <- seq(1,100)

for (rp in 1:repz) {
  print(rp)
  ex_1 <- output_matrix[output_matrix$sim_no == paste0(rp, "_1"),]
  ex_2 <- output_matrix[output_matrix$sim_no == paste0(rp, "_2"),]
  ex_3 <- output_matrix[output_matrix$sim_no == paste0(rp, "_3"),]
  pA_lw_rez_1 <- low_s2n_simp[,trup_col_namll] - ex_1[,predp_col_namll]
  pA_lw_rez_2 <- low_s2n_simp[,trup_col_namll] - ex_2[,predp_col_namll]
  pA_lw_rez_3 <- low_s2n_simp[,trup_col_namll] - ex_3[,predp_col_namll]
  pB_lw_rez_1 <- low_s2n_simp[,trup_col_naml2] - ex_1[,predp_col_naml2]
  pB_lw_rez_2 <- low_s2n_simp[,trup_col_naml2] - ex_2[,predp_col_naml2]
  pB_lw_rez_3 <- low_s2n_simp[,trup_col_naml2] - ex_3[,predp_col_naml2]
  pC_lw_rez_1 <- low_s2n_simp[,trup_col_naml3] - ex_1[,predp_col_naml3]
  pC_lw_rez_2 <- low_s2n_simp[,trup_col_naml3] - ex_2[,predp_col_naml3]
  pC_lw_rez_3 <- low_s2n_simp[,trup_col_naml3] - ex_3[,predp_col_naml3]
  seA_lw_rez_1 <- abs(pA_lw_rez_1)
  seA_lw_rez_2 <- abs(pA_lw_rez_2)
  seA_lw_rez_3 <- abs(pA_lw_rez_3)
  seB_lw_rez_1 <- abs(pB_lw_rez_1)
  seB_lw_rez_2 <- abs(pB_lw_rez_2)
  seB_lw_rez_3 <- abs(pB_lw_rez_3)
  seC_lw_rez_1 <- abs(pC_lw_rez_1)
  seC_lw_rez_2 <- abs(pC_lw_rez_2)
  seC_lw_rez_3 <- abs(pC_lw_rez_3)
  cl_lw_rez_1 <- ex_1[,"resid_true_cress_simp_lw"]
  cl_lw_rez_2 <- ex_2[,"resid_true_cress_simp_lw"]
  cl_lw_rez_3 <- ex_3[,"resid_true_cress_simp_lw"]
  #all_lw <- cbind(pA_lw_rez_1, pA_lw_rez_2, pA_lw_rez_3,
  #                pB_lw_rez_1, pB_lw_rez_2, pB_lw_rez_3,
  #                pC_lw_rez_1, pC_lw_rez_2, pC_lw_rez_3,
  #                cl_lw_rez_1, cl_lw_rez_2, cl_lw_rez_3)
  all_lw <- cbind(seA_lw_rez_1, seA_lw_rez_2, seA_lw_rez_3,
                  seB_lw_rez_1, seB_lw_rez_2, seB_lw_rez_3,
                  seC_lw_rez_1, seC_lw_rez_2, seC_lw_rez_3,
                  cl_lw_rez_1, cl_lw_rez_2, cl_lw_rez_3)
  output_cr_si_lw <- output_cr_si_lw + all_lw
  rss_lw_1 <- sum(pA_lw_rez_1^2) + sum(pB_lw_rez_1^2) + sum(pC_lw_rez_1^2)
  rss_lw_2 <- sum(pA_lw_rez_2^2) + sum(pB_lw_rez_2^2) + sum(pC_lw_rez_2^2) 
  rss_lw_3 <- sum(pA_lw_rez_3^2) + sum(pB_lw_rez_3^2) + sum(pC_lw_rez_3^2) 
  rss_cr_si[rp, 2:4] <- c(rss_lw_1, rss_lw_2, rss_lw_3)
  
  pA_hi_rez_1 <- high_s2n_simp[,trup_col_namhl] - ex_1[,predp_col_namhl]
  pA_hi_rez_2 <- high_s2n_simp[,trup_col_namhl] - ex_2[,predp_col_namhl]
  pA_hi_rez_3 <- high_s2n_simp[,trup_col_namhl] - ex_3[,predp_col_namhl]
  pB_hi_rez_1 <- high_s2n_simp[,trup_col_namh2] - ex_1[,predp_col_namh2]
  pB_hi_rez_2 <- high_s2n_simp[,trup_col_namh2] - ex_2[,predp_col_namh2]
  pB_hi_rez_3 <- high_s2n_simp[,trup_col_namh2] - ex_3[,predp_col_namh2]
  pC_hi_rez_1 <- high_s2n_simp[,trup_col_namh3] - ex_1[,predp_col_namh3]
  pC_hi_rez_2 <- high_s2n_simp[,trup_col_namh3] - ex_2[,predp_col_namh3]
  pC_hi_rez_3 <- high_s2n_simp[,trup_col_namh3] - ex_3[,predp_col_namh3]
  seA_hi_rez_1 <- abs(pA_hi_rez_1)
  seA_hi_rez_2 <- abs(pA_hi_rez_2)
  seA_hi_rez_3 <- abs(pA_hi_rez_3)
  seB_hi_rez_1 <- abs(pB_hi_rez_1)
  seB_hi_rez_2 <- abs(pB_hi_rez_2)
  seB_hi_rez_3 <- abs(pB_hi_rez_3)
  seC_hi_rez_1 <- abs(pC_hi_rez_1)
  seC_hi_rez_2 <- abs(pC_hi_rez_2)
  seC_hi_rez_3 <- abs(pC_hi_rez_3)
  cl_hi_rez_1 <- ex_1[,"resid_true_cress_simp_hi"]
  cl_hi_rez_2 <- ex_2[,"resid_true_cress_simp_hi"]
  cl_hi_rez_3 <- ex_3[,"resid_true_cress_simp_hi"]
  #all_hi <- cbind(pA_hi_rez_1, pA_hi_rez_2, pA_hi_rez_3,
  #                pB_hi_rez_1, pB_hi_rez_2, pB_hi_rez_3,
  #                pC_hi_rez_1, pC_hi_rez_2, pC_hi_rez_3,
  #                cl_hi_rez_1, cl_hi_rez_2, cl_hi_rez_3)
  all_hi <- cbind(seA_hi_rez_1, seA_hi_rez_2, seA_hi_rez_3,
                  seB_hi_rez_1, seB_hi_rez_2, seB_hi_rez_3,
                  seC_hi_rez_1, seC_hi_rez_2, seC_hi_rez_3,
                  cl_hi_rez_1, cl_hi_rez_2, cl_hi_rez_3)
  output_cr_si_hi <- output_cr_si_hi + all_hi
  
  rss_hi_1 <- sum(pA_hi_rez_1^2) + sum(pB_hi_rez_1^2) + sum(pC_hi_rez_1^2)
  rss_hi_2 <- sum(pA_hi_rez_2^2) + sum(pB_hi_rez_2^2) + sum(pC_hi_rez_2^2) 
  rss_hi_3 <- sum(pA_hi_rez_3^2) + sum(pB_hi_rez_3^2) + sum(pC_hi_rez_3^2) 
  rss_cr_si[rp, 5:7] <- c(rss_hi_1, rss_hi_2, rss_hi_3)
  
}

output_cr_si_lw <- as.data.frame(output_cr_si_lw / repz)
output_cr_si_hi <- as.data.frame(output_cr_si_hi / repz)
output_cr_si_lw$Model <- "cress"
output_cr_si_hi$Model <- "cress"
output_cr_si_lw <- cbind.data.frame(output_cr_si_lw, grid)
output_cr_si_hi <- cbind.data.frame(output_cr_si_hi, grid)

### MGCV

predp_col_namll <- paste0("predp_mgcv_simp_lw_1")
predp_col_naml2 <- paste0("predp_mgcv_simp_lw_2")
predp_col_naml3 <- paste0("predp_mgcv_simp_lw_3")
predp_col_namhl <- paste0("predp_mgcv_simp_hi_1")
predp_col_namh2 <- paste0("predp_mgcv_simp_hi_2")
predp_col_namh3 <- paste0("predp_mgcv_simp_hi_3")

output_matrix <- read.csv(paste0(file_save_dir, "mgcv_simp_mat2_new.csv"), sep=",", header=TRUE)

output_mg_si_lw <- as.data.frame(matrix(rep(0, 10000*12), ncol=12))
colnames(output_mg_si_lw) <- c("pA_lw_rez_1", "pA_lw_rez_2", "pA_lw_rez_3",
                               "pB_lw_rez_1", "pB_lw_rez_2", "pB_lw_rez_3",
                               "pC_lw_rez_1", "pC_lw_rez_2", "pC_lw_rez_3",
                               "cl_lw_rez_1", "cl_lw_rez_2", "cl_lw_rez_3")

output_mg_si_hi <- as.data.frame(matrix(rep(0, 10000*12), ncol=12))
colnames(output_mg_si_hi) <- c("pA_hi_rez_1", "pA_hi_rez_2", "pA_hi_rez_3",
                               "pB_hi_rez_1", "pB_hi_rez_2", "pB_hi_rez_3",
                               "pC_hi_rez_1", "pC_hi_rez_2", "pC_hi_rez_3",
                               "cl_hi_rez_1", "cl_hi_rez_2", "cl_hi_rez_3")

rss_mg_si <- as.data.frame(matrix(rep(0, repz*7), ncol=7))
colnames(rss_mg_si) <- c("sim", "lw_1", "lw_2", "lw_3", "hi_1", "hi2", "hi3")
rss_mg_si$sim <- seq(1,100)

for (rp in 1:repz) {
  ex_1 <- output_matrix[output_matrix$sim_no == paste0(rp, "_1"),]
  ex_2 <- output_matrix[output_matrix$sim_no == paste0(rp, "_2"),]
  ex_3 <- output_matrix[output_matrix$sim_no == paste0(rp, "_3"),]
  pA_lw_rez_1 <- low_s2n_simp[,trup_col_namll] - ex_1[,predp_col_namll]
  pA_lw_rez_2 <- low_s2n_simp[,trup_col_namll] - ex_2[,predp_col_namll]
  pA_lw_rez_3 <- low_s2n_simp[,trup_col_namll] - ex_3[,predp_col_namll]
  pB_lw_rez_1 <- low_s2n_simp[,trup_col_naml2] - ex_1[,predp_col_naml2]
  pB_lw_rez_2 <- low_s2n_simp[,trup_col_naml2] - ex_2[,predp_col_naml2]
  pB_lw_rez_3 <- low_s2n_simp[,trup_col_naml2] - ex_3[,predp_col_naml2]
  pC_lw_rez_1 <- low_s2n_simp[,trup_col_naml3] - ex_1[,predp_col_naml3]
  pC_lw_rez_2 <- low_s2n_simp[,trup_col_naml3] - ex_2[,predp_col_naml3]
  pC_lw_rez_3 <- low_s2n_simp[,trup_col_naml3] - ex_3[,predp_col_naml3]
  seA_lw_rez_1 <- abs(pA_lw_rez_1)
  seA_lw_rez_2 <- abs(pA_lw_rez_2)
  seA_lw_rez_3 <- abs(pA_lw_rez_3)
  seB_lw_rez_1 <- abs(pB_lw_rez_1)
  seB_lw_rez_2 <- abs(pB_lw_rez_2)
  seB_lw_rez_3 <- abs(pB_lw_rez_3)
  seC_lw_rez_1 <- abs(pC_lw_rez_1)
  seC_lw_rez_2 <- abs(pC_lw_rez_2)
  seC_lw_rez_3 <- abs(pC_lw_rez_3)
  cl_lw_rez_1 <- ex_1[,"resid_true_mgcv_simp_lw"]
  cl_lw_rez_2 <- ex_2[,"resid_true_mgcv_simp_lw"]
  cl_lw_rez_3 <- ex_3[,"resid_true_mgcv_simp_lw"]
  #all_lw <- cbind(pA_lw_rez_1, pA_lw_rez_2, pA_lw_rez_3,
  #                pB_lw_rez_1, pB_lw_rez_2, pB_lw_rez_3,
  #                pC_lw_rez_1, pC_lw_rez_2, pC_lw_rez_3,
  #                cl_lw_rez_1, cl_lw_rez_2, cl_lw_rez_3)
  all_lw <- cbind(seA_lw_rez_1, seA_lw_rez_2, seA_lw_rez_3,
                  seB_lw_rez_1, seB_lw_rez_2, seB_lw_rez_3,
                  seC_lw_rez_1, seC_lw_rez_2, seC_lw_rez_3,
                  cl_lw_rez_1, cl_lw_rez_2, cl_lw_rez_3)
  output_mg_si_lw <- output_mg_si_lw + all_lw
  rss_lw_1 <- sum(pA_lw_rez_1^2) + sum(pB_lw_rez_1^2) + sum(pC_lw_rez_1^2)
  rss_lw_2 <- sum(pA_lw_rez_2^2) + sum(pB_lw_rez_2^2) + sum(pC_lw_rez_2^2) 
  rss_lw_3 <- sum(pA_lw_rez_3^2) + sum(pB_lw_rez_3^2) + sum(pC_lw_rez_3^2) 
  rss_mg_si[rp, 2:4] <- c(rss_lw_1, rss_lw_2, rss_lw_3)
  
  pA_hi_rez_1 <- high_s2n_simp[,trup_col_namhl] - ex_1[,predp_col_namhl]
  pA_hi_rez_2 <- high_s2n_simp[,trup_col_namhl] - ex_2[,predp_col_namhl]
  pA_hi_rez_3 <- high_s2n_simp[,trup_col_namhl] - ex_3[,predp_col_namhl]
  pB_hi_rez_1 <- high_s2n_simp[,trup_col_namh2] - ex_1[,predp_col_namh2]
  pB_hi_rez_2 <- high_s2n_simp[,trup_col_namh2] - ex_2[,predp_col_namh2]
  pB_hi_rez_3 <- high_s2n_simp[,trup_col_namh2] - ex_3[,predp_col_namh2]
  pC_hi_rez_1 <- high_s2n_simp[,trup_col_namh3] - ex_1[,predp_col_namh3]
  pC_hi_rez_2 <- high_s2n_simp[,trup_col_namh3] - ex_2[,predp_col_namh3]
  pC_hi_rez_3 <- high_s2n_simp[,trup_col_namh3] - ex_3[,predp_col_namh3]
  seA_hi_rez_1 <- abs(pA_hi_rez_1)
  seA_hi_rez_2 <- abs(pA_hi_rez_2)
  seA_hi_rez_3 <- abs(pA_hi_rez_3)
  seB_hi_rez_1 <- abs(pB_hi_rez_1)
  seB_hi_rez_2 <- abs(pB_hi_rez_2)
  seB_hi_rez_3 <- abs(pB_hi_rez_3)
  seC_hi_rez_1 <- abs(pC_hi_rez_1)
  seC_hi_rez_2 <- abs(pC_hi_rez_2)
  seC_hi_rez_3 <- abs(pC_hi_rez_3)
  cl_hi_rez_1 <- ex_1[,"resid_true_mgcv_simp_hi"]
  cl_hi_rez_2 <- ex_2[,"resid_true_mgcv_simp_hi"]
  cl_hi_rez_3 <- ex_3[,"resid_true_mgcv_simp_hi"]
  #all_hi <- cbind(pA_hi_rez_1, pA_hi_rez_2, pA_hi_rez_3,
  #                pB_hi_rez_1, pB_hi_rez_2, pB_hi_rez_3,
  #                pC_hi_rez_1, pC_hi_rez_2, pC_hi_rez_3,
  #                cl_hi_rez_1, cl_hi_rez_2, cl_hi_rez_3)
  all_hi <- cbind(seA_hi_rez_1, seA_hi_rez_2, seA_hi_rez_3,
                  seB_hi_rez_1, seB_hi_rez_2, seB_hi_rez_3,
                  seC_hi_rez_1, seC_hi_rez_2, seC_hi_rez_3,
                  cl_hi_rez_1, cl_hi_rez_2, cl_hi_rez_3)
  output_mg_si_hi <- output_mg_si_hi + all_hi
  
  rss_hi_1 <- sum(pA_hi_rez_1^2) + sum(pB_hi_rez_1^2) + sum(pC_hi_rez_1^2)
  rss_hi_2 <- sum(pA_hi_rez_2^2) + sum(pB_hi_rez_2^2) + sum(pC_hi_rez_2^2) 
  rss_hi_3 <- sum(pA_hi_rez_3^2) + sum(pB_hi_rez_3^2) + sum(pC_hi_rez_3^2) 
  rss_mg_si[rp, 5:7] <- c(rss_hi_1, rss_hi_2, rss_hi_3)
  
}

output_mg_si_lw <- as.data.frame(output_mg_si_lw / repz)
output_mg_si_hi <- as.data.frame(output_mg_si_hi / repz)
output_mg_si_lw$Model <- "mgcv"
output_mg_si_hi$Model <- "mgcv"
output_mg_si_lw <- cbind.data.frame(output_mg_si_lw, grid)
output_mg_si_hi <- cbind.data.frame(output_mg_si_hi, grid)

### SVM

predp_col_namll <- paste0("predp_svm_simp_lw_1")
predp_col_naml2 <- paste0("predp_svm_simp_lw_2")
predp_col_naml3 <- paste0("predp_svm_simp_lw_3")
predp_col_namhl <- paste0("predp_svm_simp_hi_1")
predp_col_namh2 <- paste0("predp_svm_simp_hi_2")
predp_col_namh3 <- paste0("predp_svm_simp_hi_3")

output_matrix <- read.csv(paste0(file_save_dir, "svm_simp_mat2_new.csv"), sep=",", header=TRUE)

output_sv_si_lw <- as.data.frame(matrix(rep(0, 10000*12), ncol=12))
colnames(output_sv_si_lw) <- c("pA_lw_rez_1", "pA_lw_rez_2", "pA_lw_rez_3",
                               "pB_lw_rez_1", "pB_lw_rez_2", "pB_lw_rez_3",
                               "pC_lw_rez_1", "pC_lw_rez_2", "pC_lw_rez_3",
                               "cl_lw_rez_1", "cl_lw_rez_2", "cl_lw_rez_3")

output_sv_si_hi <- as.data.frame(matrix(rep(0, 10000*12), ncol=12))
colnames(output_sv_si_hi) <- c("pA_hi_rez_1", "pA_hi_rez_2", "pA_hi_rez_3",
                               "pB_hi_rez_1", "pB_hi_rez_2", "pB_hi_rez_3",
                               "pC_hi_rez_1", "pC_hi_rez_2", "pC_hi_rez_3",
                               "cl_hi_rez_1", "cl_hi_rez_2", "cl_hi_rez_3")

rss_sv_si <- as.data.frame(matrix(rep(0, repz*7), ncol=7))
colnames(rss_sv_si) <- c("sim", "lw_1", "lw_2", "lw_3", "hi_1", "hi2", "hi3")
rss_sv_si$sim <- seq(1,100)

for (rp in 1:repz) {
  ex_1 <- output_matrix[output_matrix$sim_no == paste0(rp, "_1"),]
  ex_2 <- output_matrix[output_matrix$sim_no == paste0(rp, "_2"),]
  ex_3 <- output_matrix[output_matrix$sim_no == paste0(rp, "_3"),]
  pA_lw_rez_1 <- low_s2n_simp[,trup_col_namll] - ex_1[,predp_col_namll]
  pA_lw_rez_2 <- low_s2n_simp[,trup_col_namll] - ex_2[,predp_col_namll]
  pA_lw_rez_3 <- low_s2n_simp[,trup_col_namll] - ex_3[,predp_col_namll]
  pB_lw_rez_1 <- low_s2n_simp[,trup_col_naml2] - ex_1[,predp_col_naml2]
  pB_lw_rez_2 <- low_s2n_simp[,trup_col_naml2] - ex_2[,predp_col_naml2]
  pB_lw_rez_3 <- low_s2n_simp[,trup_col_naml2] - ex_3[,predp_col_naml2]
  pC_lw_rez_1 <- low_s2n_simp[,trup_col_naml3] - ex_1[,predp_col_naml3]
  pC_lw_rez_2 <- low_s2n_simp[,trup_col_naml3] - ex_2[,predp_col_naml3]
  pC_lw_rez_3 <- low_s2n_simp[,trup_col_naml3] - ex_3[,predp_col_naml3]
  seA_lw_rez_1 <- abs(pA_lw_rez_1)
  seA_lw_rez_2 <- abs(pA_lw_rez_2)
  seA_lw_rez_3 <- abs(pA_lw_rez_3)
  seB_lw_rez_1 <- abs(pB_lw_rez_1)
  seB_lw_rez_2 <- abs(pB_lw_rez_2)
  seB_lw_rez_3 <- abs(pB_lw_rez_3)
  seC_lw_rez_1 <- abs(pC_lw_rez_1)
  seC_lw_rez_2 <- abs(pC_lw_rez_2)
  seC_lw_rez_3 <- abs(pC_lw_rez_3)
  cl_lw_rez_1 <- ex_1[,"resid_true_svm_simp_lw"]
  cl_lw_rez_2 <- ex_2[,"resid_true_svm_simp_lw"]
  cl_lw_rez_3 <- ex_3[,"resid_true_svm_simp_lw"]
  #all_lw <- cbind(pA_lw_rez_1, pA_lw_rez_2, pA_lw_rez_3,
  #                pB_lw_rez_1, pB_lw_rez_2, pB_lw_rez_3,
  #                pC_lw_rez_1, pC_lw_rez_2, pC_lw_rez_3,
  #                cl_lw_rez_1, cl_lw_rez_2, cl_lw_rez_3)
  all_lw <- cbind(seA_lw_rez_1, seA_lw_rez_2, seA_lw_rez_3,
                  seB_lw_rez_1, seB_lw_rez_2, seB_lw_rez_3,
                  seC_lw_rez_1, seC_lw_rez_2, seC_lw_rez_3,
                  cl_lw_rez_1, cl_lw_rez_2, cl_lw_rez_3)
  output_sv_si_lw <- output_sv_si_lw + all_lw
  
  rss_lw_1 <- sum(pA_lw_rez_1^2) + sum(pB_lw_rez_1^2) + sum(pC_lw_rez_1^2)
  rss_lw_2 <- sum(pA_lw_rez_2^2) + sum(pB_lw_rez_2^2) + sum(pC_lw_rez_2^2) 
  rss_lw_3 <- sum(pA_lw_rez_3^2) + sum(pB_lw_rez_3^2) + sum(pC_lw_rez_3^2) 
  rss_sv_si[rp, 2:4] <- c(rss_lw_1, rss_lw_2, rss_lw_3)
  
  pA_hi_rez_1 <- high_s2n_simp[,trup_col_namhl] - ex_1[,predp_col_namhl]
  pA_hi_rez_2 <- high_s2n_simp[,trup_col_namhl] - ex_2[,predp_col_namhl]
  pA_hi_rez_3 <- high_s2n_simp[,trup_col_namhl] - ex_3[,predp_col_namhl]
  pB_hi_rez_1 <- high_s2n_simp[,trup_col_namh2] - ex_1[,predp_col_namh2]
  pB_hi_rez_2 <- high_s2n_simp[,trup_col_namh2] - ex_2[,predp_col_namh2]
  pB_hi_rez_3 <- high_s2n_simp[,trup_col_namh2] - ex_3[,predp_col_namh2]
  pC_hi_rez_1 <- high_s2n_simp[,trup_col_namh3] - ex_1[,predp_col_namh3]
  pC_hi_rez_2 <- high_s2n_simp[,trup_col_namh3] - ex_2[,predp_col_namh3]
  pC_hi_rez_3 <- high_s2n_simp[,trup_col_namh3] - ex_3[,predp_col_namh3]
  seA_hi_rez_1 <- abs(pA_hi_rez_1)
  seA_hi_rez_2 <- abs(pA_hi_rez_2)
  seA_hi_rez_3 <- abs(pA_hi_rez_3)
  seB_hi_rez_1 <- abs(pB_hi_rez_1)
  seB_hi_rez_2 <- abs(pB_hi_rez_2)
  seB_hi_rez_3 <- abs(pB_hi_rez_3)
  seC_hi_rez_1 <- abs(pC_hi_rez_1)
  seC_hi_rez_2 <- abs(pC_hi_rez_2)
  seC_hi_rez_3 <- abs(pC_hi_rez_3)
  cl_hi_rez_1 <- ex_1[,"resid_true_svm_simp_hi"]
  cl_hi_rez_2 <- ex_2[,"resid_true_svm_simp_hi"]
  cl_hi_rez_3 <- ex_3[,"resid_true_svm_simp_hi"]
  #all_hi <- cbind(pA_hi_rez_1, pA_hi_rez_2, pA_hi_rez_3,
  #                pB_hi_rez_1, pB_hi_rez_2, pB_hi_rez_3,
  #                pC_hi_rez_1, pC_hi_rez_2, pC_hi_rez_3,
  #                cl_hi_rez_1, cl_hi_rez_2, cl_hi_rez_3)
  all_hi <- cbind(seA_hi_rez_1, seA_hi_rez_2, seA_hi_rez_3,
                  seB_hi_rez_1, seB_hi_rez_2, seB_hi_rez_3,
                  seC_hi_rez_1, seC_hi_rez_2, seC_hi_rez_3,
                  cl_hi_rez_1, cl_hi_rez_2, cl_hi_rez_3)
  output_sv_si_hi <- output_sv_si_hi + all_hi
  
  rss_hi_1 <- sum(pA_hi_rez_1^2) + sum(pB_hi_rez_1^2) + sum(pC_hi_rez_1^2)
  rss_hi_2 <- sum(pA_hi_rez_2^2) + sum(pB_hi_rez_2^2) + sum(pC_hi_rez_2^2) 
  rss_hi_3 <- sum(pA_hi_rez_3^2) + sum(pB_hi_rez_3^2) + sum(pC_hi_rez_3^2) 
  rss_sv_si[rp, 5:7] <- c(rss_hi_1, rss_hi_2, rss_hi_3)
  
}

output_sv_si_lw <- as.data.frame(output_sv_si_lw / repz)
output_sv_si_hi <- as.data.frame(output_sv_si_hi / repz)
output_sv_si_lw$Model <- "svm"
output_sv_si_hi$Model <- "svm"
output_sv_si_lw <- cbind.data.frame(output_sv_si_lw, grid)
output_sv_si_hi <- cbind.data.frame(output_sv_si_hi, grid)

write.csv(rss_cr_si, paste0(img_save_dir, "rss_cress_simp.csv"))
write.csv(rss_mg_si, paste0(img_save_dir, "rss_mgcv_simp.csv"))
write.csv(rss_sv_si, paste0(img_save_dir, "rss_svm_simp.csv"))

library(scales)

graph_cress_A_simp_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_simp_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_simp_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_simp_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_simp_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_simp_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_simp_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_simp_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_simp_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_simp_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_simp_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_simp_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_simp_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_simp_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_simp_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_simp_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_simp_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_simp_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_simp_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_simp_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_simp_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_simp_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_simp_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_simp_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_simp_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_simp_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_simp_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))



### High noise


graph_cress_A_simp_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_simp_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_simp_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_simp_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_simp_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_simp_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_simp_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_simp_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_simp_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_simp_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_simp_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_simp_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_simp_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_simp_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_simp_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_simp_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_simp_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_simp_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))


graph_svm_A_simp_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_simp_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_simp_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_simp_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_simp_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_simp_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_simp_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_simp_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_simp_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))





### Classes

graph_cress_cl_simp_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_simp_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_simp_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_simp_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_simp_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_simp_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_simp_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_simp_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_simp_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_simp_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_simp_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_simp_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_simp_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_simp_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_simp_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_simp_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_simp_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_simp_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))



### RSS

output_cr_si_lw$lw_rss_1 <- (output_cr_si_lw$pA_lw_rez_1 + output_cr_si_lw$pB_lw_rez_1 + output_cr_si_lw$pC_lw_rez_1) / 3
output_cr_si_lw$lw_rss_2 <- (output_cr_si_lw$pA_lw_rez_2 + output_cr_si_lw$pB_lw_rez_2 + output_cr_si_lw$pC_lw_rez_2) / 3
output_cr_si_lw$lw_rss_3 <- (output_cr_si_lw$pA_lw_rez_3 + output_cr_si_lw$pB_lw_rez_3 + output_cr_si_lw$pC_lw_rez_3) / 3
output_cr_si_hi$hi_rss_1 <- (output_cr_si_hi$pA_hi_rez_1 + output_cr_si_hi$pB_hi_rez_1 + output_cr_si_hi$pC_hi_rez_1) / 3
output_cr_si_hi$hi_rss_2 <- (output_cr_si_hi$pA_hi_rez_2 + output_cr_si_hi$pB_hi_rez_2 + output_cr_si_hi$pC_hi_rez_2) / 3
output_cr_si_hi$hi_rss_3 <- (output_cr_si_hi$pA_hi_rez_3 + output_cr_si_hi$pB_hi_rez_3 + output_cr_si_hi$pC_hi_rez_3) / 3
output_mg_si_lw$lw_rss_1 <- (output_mg_si_lw$pA_lw_rez_1 + output_mg_si_lw$pB_lw_rez_1 + output_mg_si_lw$pC_lw_rez_1) / 3
output_mg_si_lw$lw_rss_2 <- (output_mg_si_lw$pA_lw_rez_2 + output_mg_si_lw$pB_lw_rez_2 + output_mg_si_lw$pC_lw_rez_2) / 3
output_mg_si_lw$lw_rss_3 <- (output_mg_si_lw$pA_lw_rez_3 + output_mg_si_lw$pB_lw_rez_3 + output_mg_si_lw$pC_lw_rez_3) / 3
output_mg_si_hi$hi_rss_1 <- (output_mg_si_hi$pA_hi_rez_1 + output_mg_si_hi$pB_hi_rez_1 + output_mg_si_hi$pC_hi_rez_1) / 3
output_mg_si_hi$hi_rss_2 <- (output_mg_si_hi$pA_hi_rez_2 + output_mg_si_hi$pB_hi_rez_2 + output_mg_si_hi$pC_hi_rez_2) / 3
output_mg_si_hi$hi_rss_3 <- (output_mg_si_hi$pA_hi_rez_3 + output_mg_si_hi$pB_hi_rez_3 + output_mg_si_hi$pC_hi_rez_3) / 3
output_sv_si_lw$lw_rss_1 <- (output_sv_si_lw$pA_lw_rez_1 + output_sv_si_lw$pB_lw_rez_1 + output_sv_si_lw$pC_lw_rez_1) / 3
output_sv_si_lw$lw_rss_2 <- (output_sv_si_lw$pA_lw_rez_2 + output_sv_si_lw$pB_lw_rez_2 + output_sv_si_lw$pC_lw_rez_2) / 3
output_sv_si_lw$lw_rss_3 <- (output_sv_si_lw$pA_lw_rez_3 + output_sv_si_lw$pB_lw_rez_3 + output_sv_si_lw$pC_lw_rez_3) / 3
output_sv_si_hi$hi_rss_1 <- (output_sv_si_hi$pA_hi_rez_1 + output_sv_si_hi$pB_hi_rez_1 + output_sv_si_hi$pC_hi_rez_1) / 3
output_sv_si_hi$hi_rss_2 <- (output_sv_si_hi$pA_hi_rez_2 + output_sv_si_hi$pB_hi_rez_2 + output_sv_si_hi$pC_hi_rez_2) / 3
output_sv_si_hi$hi_rss_3 <- (output_sv_si_hi$pA_hi_rez_3 + output_sv_si_hi$pB_hi_rez_3 + output_sv_si_hi$pC_hi_rez_3) / 3

graph_cress_rss_simp_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_simp_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_simp_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_simp_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_simp_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_simp_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_simp_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_simp_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_simp_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_simp_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_simp_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_simp_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_simp_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_simp_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_simp_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_simp_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_simp_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_simp_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))


png(filename=paste0(img_save_dir, "GraphCrASimpLw1.png"))
graph_cress_A_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBSimpLw1.png"))
graph_cress_B_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCSimpLw1.png"))
graph_cress_C_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrASimpLw2.png"))
graph_cress_A_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBSimpLw2.png"))
graph_cress_B_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCSimpLw2.png"))
graph_cress_C_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrASimpLw3.png"))
graph_cress_A_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBSimpLw3.png"))
graph_cress_B_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCSimpLw3.png"))
graph_cress_C_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrASimpHi1.png"))
graph_cress_A_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBSimpHi1.png"))
graph_cress_B_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCSimpHi1.png"))
graph_cress_C_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrASimpHi2.png"))
graph_cress_A_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBSimpHi2.png"))
graph_cress_B_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCSimpHi2.png"))
graph_cress_C_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrASimpHi3.png"))
graph_cress_A_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBSimpHi3.png"))
graph_cress_B_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCSimpHi3.png"))
graph_cress_C_simp_hi3
dev.off()

png(filename=paste0(img_save_dir, "GraphMgASimpLw1.png"))
graph_mgcv_A_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBSimpLw1.png"))
graph_mgcv_B_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCSimpLw1.png"))
graph_mgcv_C_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgASimpLw2.png"))
graph_mgcv_A_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBSimpLw2.png"))
graph_mgcv_B_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCSimpLw2.png"))
graph_mgcv_C_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgASimpLw3.png"))
graph_mgcv_A_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBSimpLw3.png"))
graph_mgcv_B_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCSimpLw3.png"))
graph_mgcv_C_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgASimpHi1.png"))
graph_mgcv_A_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBSimpHi1.png"))
graph_mgcv_B_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCSimpHi1.png"))
graph_mgcv_C_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgASimpHi2.png"))
graph_mgcv_A_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBSimpHi2.png"))
graph_mgcv_B_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCSimpHi2.png"))
graph_mgcv_C_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgASimpHi3.png"))
graph_mgcv_A_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBSimpHi3.png"))
graph_mgcv_B_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCSimpHi3.png"))
graph_mgcv_C_simp_hi3
dev.off()

# combine into one
simp_lw1 <- (graph_cress_A_simp_lw1 | graph_cress_B_simp_lw1 | 
  graph_cress_C_simp_lw1) / (graph_mgcv_A_simp_lw1 |
  graph_mgcv_B_simp_lw1 | graph_mgcv_C_simp_lw3) /
  (graph_svm_A_simp_lw1 | graph_svm_B_simp_lw1 |
  graph_svm_C_simp_lw1) + plot_layout(nrow = 3, byrow = TRUE)
png(filename=paste0(img_save_dir, "GraphSimpLw1.png"))
simp_lw1
dev.off()

wrap_plots(graph_cress_A_simp_lw1, graph_cress_B_simp_lw1, 
  graph_cress_C_simp_lw1, graph_mgcv_A_simp_lw1,
  graph_mgcv_B_simp_lw1, graph_mgcv_C_simp_lw3,
  graph_svm_A_simp_lw1, graph_svm_B_simp_lw1, 
  graph_svm_C_simp_lw1)

png(filename=paste0(img_save_dir, "GraphSvASimpLw1.png"))
graph_svm_A_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBSimpLw1.png"))
graph_svm_B_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCSimpLw1.png"))
graph_svm_C_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvASimpLw2.png"))
graph_svm_A_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBSimpLw2.png"))
graph_svm_B_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCSimpLw2.png"))
graph_svm_C_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvASimpLw3.png"))
graph_svm_A_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBSimpLw3.png"))
graph_svm_B_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCSimpLw3.png"))
graph_svm_C_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvASimpHi1.png"))
graph_svm_A_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBSimpHi1.png"))
graph_svm_B_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCSimpHi1.png"))
graph_svm_C_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvASimpHi2.png"))
graph_svm_A_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBSimpHi2.png"))
graph_svm_B_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCSimpHi2.png"))
graph_svm_C_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvASimpHi3.png"))
graph_svm_A_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBSimpHi3.png"))
graph_svm_B_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCSimpHi3.png"))
graph_svm_C_simp_hi3
dev.off()








png(filename=paste0(img_save_dir, "GraphCrClSimpLw1.png"))
graph_cress_cl_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClSimpLw2.png"))
graph_cress_cl_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClSimpLw3.png"))
graph_cress_cl_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClSimpLw1.png"))
graph_mgcv_cl_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClSimpLw2.png"))
graph_mgcv_cl_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClSimpLw3.png"))
graph_mgcv_cl_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClSimpLw1.png"))
graph_svm_cl_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClSimpLw2.png"))
graph_svm_cl_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClSimpLw3.png"))
graph_svm_cl_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClSimpHi1.png"))
graph_cress_cl_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClSimpHi2.png"))
graph_cress_cl_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClSimpHi3.png"))
graph_cress_cl_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClSimpHi1.png"))
graph_mgcv_cl_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClSimpHi2.png"))
graph_mgcv_cl_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClSimpHi3.png"))
graph_mgcv_cl_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClSimpHi1.png"))
graph_svm_cl_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClSimpHi2.png"))
graph_svm_cl_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClSimpHi3.png"))
graph_svm_cl_simp_hi3
dev.off()



png(filename=paste0(img_save_dir, "GraphCrRssSimpLw1.png"))
graph_cress_rss_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssSimpLw2.png"))
graph_cress_rss_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssSimpLw3.png"))
graph_cress_rss_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssSimpLw1.png"))
graph_mgcv_rss_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssSimpLw2.png"))
graph_mgcv_rss_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssSimpLw3.png"))
graph_mgcv_rss_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssSimpLw1.png"))
graph_svm_rss_simp_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssSimpLw2.png"))
graph_svm_rss_simp_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssSimpLw3.png"))
graph_svm_rss_simp_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssSimpHi1.png"))
graph_cress_rss_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssSimpHi2.png"))
graph_cress_rss_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssSimpHi3.png"))
graph_cress_rss_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssSimpHi1.png"))
graph_mgcv_rss_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssSimpHi2.png"))
graph_mgcv_rss_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssSimpHi3.png"))
graph_mgcv_rss_simp_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssSimpHi1.png"))
graph_svm_rss_simp_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssSimpHi2.png"))
graph_svm_rss_simp_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssSimpHi3.png"))
graph_svm_rss_simp_hi3
dev.off()







# Worst and Best results for CRESS

tab_out_disagg <- read.csv(paste0(img_save_dir, "disagg_comparison_table_new.csv"))

tab_mins <- tab_out_disagg %>% group_by(model, dataset, noise, noPoints) %>% summarise(which.min(accuracy_to_data_all))
tab_maxs <- tab_out_disagg %>% group_by(model, dataset, noise, noPoints) %>% summarise(which.max(accuracy_to_data_all))

colnames(tab_mins) <- c("model", "dataset", "noise", "noPoints", "minrow")
colnames(tab_maxs) <- c("model", "dataset", "noise", "noPoints", "maxrow")

dataset_rows1 <- read.csv(paste0(file_save_dir,"dataset_rows1.csv"))
dataset_rows2 <- read.csv(paste0(file_save_dir,"dataset_rows2.csv"))
dataset_rows3 <- read.csv(paste0(file_save_dir,"dataset_rows3.csv"))
dataset_rows1i <- read.csv(paste0(file_save_dir,"dataset_rows1_island.csv"))
dataset_rows2i <- read.csv(paste0(file_save_dir,"dataset_rows2_island.csv"))
dataset_rows3i <- read.csv(paste0(file_save_dir,"dataset_rows3_island.csv"))

deets_simp <- read.csv(paste0(file_in_dir,"cress_simp_deets_list2_new.csv"))
deets_comp <- read.csv(paste0(file_save_dir,"cress_comp_deets_list2_new.csv"))
deets_island <- read.csv(paste0(file_save_dir,"cress_island_deets_list2_new.csv"))

outmat_nam_save <- ""

for (rw in 1:nrow(tab_mins)){
  print(rw)
  
  row_mn <- tab_mins[rw,]
  row_mx <- tab_maxs[rw,]
  
  if (row_mn$dataset == 'island') {
    grid <- read.csv(paste0(file_save_dir,"location_grid_island.csv"))
  } else {
    grid <- read.csv(paste0(file_save_dir,"location_grid.csv"))
  }
  
  if (row_mn$model == 'cress') {
    file_dir <- file_in_dir
  } else {
    file_dir <- file_save_dir
  }

  outmat_name <- paste0(file_dir, row_mn$model, "_", row_mn$dataset, "_mat2_new.csv")
  if (outmat_name != outmat_nam_save){
    output_matrix <- read.csv(outmat_name, sep=",", header=TRUE)
    outmat_nam_save <- outmat_name
  }
  data_out <- read.csv(paste0(file_save_dir, row_mn$dataset, "_input_data.csv"), sep=",", header=TRUE)

  if (row_mn$model == 'cress') {
    listout_name <- paste0(file_in_dir, "cress_", row_mn$dataset, "_deets_list2_new.csv")
    listout <- read.csv(listout_name, sep=",", header=TRUE)

    matout_name <- paste0(file_in_dir, "cress_", row_mn$dataset, "_deets_mat2_new.csv")
    matout <- read.csv(matout_name, sep=",", header=TRUE)

    coefz_name <- paste0(file_in_dir, "cress_", row_mn$dataset, "_coefs_new.csv")
    coefzout <- read.csv(coefz_name, sep=",", header=TRUE)
    coefzout <- as.matrix(coefzout)
    coefzout <- matrix(coefzout, ncol=2, byrow=T)
    out_deets <- list()
    totknots <- 0
    cfz <- 0
    
    for (md in 1:nrow(listout)) {
      nam <- as.character(listout$namez[md])
      nknotz <- listout$nknt_out[md]
      knotz <- matout[(totknots+1):(totknots+nknotz), 1:2]
      knotz <- as.matrix(knotz)
      intz <- coefzout[(cfz+1),]
      intz <- as.matrix(intz)
      cfz <- cfz+1
      coefz <- coefzout[(cfz+1):(cfz+nknotz),]
      coefz <- matrix(coefz, nrow=nknotz, byrow=TRUE)
      cfz <- cfz + nknotz
      rownames(knotz) <- NULL
      colnames(knotz) <- NULL
      radz <- matout[(totknots+1):(totknots+nknotz), 3]
      totknots <- totknots + nknotz
      list_out <- list("knots"=knotz, "no_knots"=nknotz, "radii"=radz, "intercept"=intz, "coeffs"=coefz)
      out_deets[[nam]] <- list_out
    }
  }

  ex_1 <- output_matrix[output_matrix$sim_no == paste0(row_mn$minrow, "_1"),]
  ex_2 <- output_matrix[output_matrix$sim_no == paste0(row_mn$minrow, "_2"),]
  ex_3 <- output_matrix[output_matrix$sim_no == paste0(row_mn$minrow, "_3"),]
  
  predp_col_nam1 <- paste0("predp_", row_mn$model, "_", row_mn$dataset, "_", row_mn$noise, "_1")
  predp_col_nam2 <- paste0("predp_", row_mn$model, "_", row_mn$dataset, "_", row_mn$noise, "_2")
  predp_col_nam3 <- paste0("predp_", row_mn$model, "_", row_mn$dataset, "_", row_mn$noise, "_3")
  trup_col_nam1 <- paste0("truep_", row_mn$dataset, "_", row_mn$noise, "_1")
  trup_col_nam2 <- paste0("truep_", row_mn$dataset, "_", row_mn$noise, "_2")
  trup_col_nam3 <- paste0("truep_", row_mn$dataset, "_", row_mn$noise, "_3")
  tru_col_nam_cl <- paste0("truthcl_", row_mn$dataset, "_", row_mn$noise)
  bias_1 <- ex_1[,predp_col_nam1] - data_out[,trup_col_nam1]
  bias_2 <- ex_2[,predp_col_nam2] - data_out[,trup_col_nam2]
  bias_3 <- ex_3[,predp_col_nam3] - data_out[,trup_col_nam3]
  
  bias_comb <- cbind.data.frame(grid$xx, grid$yy, bias_1, bias_2, bias_3, 
    ex_1[,predp_col_nam1], ex_2[,predp_col_nam2], ex_3[,predp_col_nam3], 
    data_out[,trup_col_nam1], data_out[,trup_col_nam2], data_out[,trup_col_nam3])
  colnames(bias_comb) <- c("xx", "yy", "bias_1", "bias_2", "bias_3", 
    "pred_1", "pred_2", "pred_3", "true_1", "true_2", "true_3")
  
  ex_1_mx <- output_matrix[output_matrix$sim_no == paste0(row_mx$maxrow, "_1"),]
  ex_2_mx <- output_matrix[output_matrix$sim_no == paste0(row_mx$maxrow, "_2"),]
  ex_3_mx <- output_matrix[output_matrix$sim_no == paste0(row_mx$maxrow, "_3"),]
  
  bias_1_mx <- ex_1_mx[,predp_col_nam1] - data_out[,trup_col_nam1]
  bias_2_mx <- ex_2_mx[,predp_col_nam2] - data_out[,trup_col_nam2]
  bias_3_mx <- ex_3_mx[,predp_col_nam3] - data_out[,trup_col_nam3]
  
  bias_comb_mx <- cbind.data.frame(grid$xx, grid$yy, bias_1_mx, bias_2_mx, bias_3_mx, 
    ex_1_mx[,predp_col_nam1], ex_2_mx[,predp_col_nam2], ex_3_mx[,predp_col_nam3], 
    data_out[,trup_col_nam1], data_out[,trup_col_nam2], data_out[,trup_col_nam3])
  colnames(bias_comb_mx) <- c("xx", "yy", "bias_1", "bias_2", "bias_3", 
                           "pred_1", "pred_2", "pred_3", "true_1", "true_2", "true_3")
  
  if (row_mn$noise == "lw") {
    noize <- "low"
  } else {
    noize <- "high"
  }

  if (row_mn$model == 'cress') {
    knotz <- out_deets[[paste0("cress", noize, "1_", row_mn$minrow)]]$knots
    radii <- out_deets[[paste0("cress", noize, "1_", row_mn$minrow)]]$radii
    coeff <- out_deets[[paste0("cress", noize, "1_", row_mn$minrow)]]$coeffs
    knt_comb <- cbind(knotz, radii)
    colnames(knt_comb) <- c("xx", "yy", "radii")
    knt_comb <- as.data.frame(knt_comb)
    knt_comb$rad <- 1 / knt_comb$radii
    lab_pos <- pmin(knotz, 1)
    lab_pos <- pmax(lab_pos, 0)
    lab_pos <- cbind(lab_pos, round(coeff, 2))
    lab_pos <- data.frame(lab_pos)
    colnames(lab_pos) <- c("xx", "yy", "AA", "BB")
    
    knotz_mx <- out_deets[[paste0("cress", noize, "1_", row_mx$maxrow)]]$knots
    radii_mx <- out_deets[[paste0("cress", noize, "1_", row_mx$maxrow)]]$radii
    coeff_mx <- out_deets[[paste0("cress", noize, "1_", row_mx$maxrow)]]$coeffs
    knt_comb_mx <- cbind(knotz_mx, radii_mx)
    colnames(knt_comb_mx) <- c("xx", "yy", "radii")
    knt_comb_mx <- as.data.frame(knt_comb_mx)
    knt_comb_mx$rad <- 1 / knt_comb_mx$radii
    lab_pos_mx <- pmin(knotz_mx, 1)
    lab_pos_mx <- pmax(lab_pos_mx, 0)
    lab_pos_mx <- cbind(lab_pos_mx, round(coeff_mx, 2))
    lab_pos_mx <- data.frame(lab_pos_mx)
    colnames(lab_pos_mx) <- c("xx", "yy", "AA", "BB")
  }

  if (row_mn$noPoints == 50){
    if (row_mn$dataset == "island"){
      dat_rws <- unlist(dataset_rows1i[row_mn$minrow,])
    } else {
      dat_rws <- unlist(dataset_rows1[row_mn$minrow,])
    }
    dat_pts <- grid[dat_rws,]
    dat_pts$clazz <- data_out[dat_rws,tru_col_nam_cl]
  } else if (row_mn$noPoints == 500){
    if (row_mn$dataset == "island"){
      dat_rws <- unlist(dataset_rows2i[row_mn$minrow,])
    } else {
      dat_rws <- unlist(dataset_rows2[row_mn$minrow,])
    }
    dat_pts <- grid[dat_rws,]
    dat_pts$clazz <- data_out[dat_rws,tru_col_nam_cl]
  } else {
    if (row_mn$dataset == "island"){
      dat_rws <- unlist(dataset_rows3i[row_mn$minrow,])
    } else {
      dat_rws <- unlist(dataset_rows3[row_mn$minrow,])
    }
    dat_pts <- grid[dat_rws,]
    dat_pts$clazz <- data_out[dat_rws,tru_col_nam_cl]
  }

  if (row_mx$noPoints == 50){
    if (row_mx$dataset == "island"){
      dat_rws_mx <- unlist(dataset_rows1i[row_mn$maxrow,])
    } else {
      dat_rws_mx <- unlist(dataset_rows1[row_mx$maxrow,])
    }
    dat_pts_mx <- grid[dat_rws_mx,]
    dat_pts_mx$clazz <- data_out[dat_rws_mx,tru_col_nam_cl]
  } else if (row_mx$noPoints == 500){
    if (row_mx$dataset == "island"){
      dat_rws_mx <- unlist(dataset_rows2i[row_mn$maxrow,])
    } else {
      dat_rws_mx <- unlist(dataset_rows2[row_mx$maxrow,])
    }
    dat_pts_mx <- grid[dat_rws_mx,]
    dat_pts_mx$clazz <- data_out[dat_rws_mx,tru_col_nam_cl]
  } else {
    if (row_mx$dataset == "island"){
      dat_rws_mx <- unlist(dataset_rows3i[row_mn$maxrow,])
    } else {
      dat_rws_mx <- unlist(dataset_rows3[row_mx$maxrow,])
    }
    dat_pts_mx <- grid[dat_rws_mx,]
    dat_pts_mx$clazz <- data_out[dat_rws_mx,tru_col_nam_cl]
  }

  pred_pts <- ggplot(dat_pts, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
    # ggtitle(paste("GFRC Distance to any water")) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_point(aes(colour=clazz), shape=15, size=2) +
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.position="right")
  pred_pts_mx <- ggplot(dat_pts_mx, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
    # ggtitle(paste("GFRC Distance to any water")) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_point(aes(colour=clazz), shape=15, size=2) +
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.position="right")
  
  true_A <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Truth for class A ", row_mn$noise, " noise case, ", row_mn$model)) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=true_1), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  true_B <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Truth for class B ", row_mn$noise, " noise case, ", row_mn$model)) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=true_2), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  true_C <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Truth for class C ", row_mn$noise, " noise case, ", row_mn$model)) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=true_3), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  
  graph_A <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class A ", row_mn$noise, " noise case, ", row_mn$model, ", ", row_mn$noPoints, " datapoints, simulation no ", row_mn$minrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=pred_1), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb)) {
      graph_A <- graph_A + annotate("path",x=knt_comb$xx[rr]+knt_comb$rad[rr]*cos(seq(0,2*pi,length.out=100)),
                y=knt_comb$yy[rr]+knt_comb$rad[rr]*sin(seq(0,2*pi,length.out=100))) +
        geom_text(aes(x=xx, y=yy, label=AA), data=lab_pos)
    }
  }
  graph_B <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class B ", row_mn$noise, " noise case, ", row_mn$model, ", ", row_mn$noPoints, " datapoints, simulation no ", row_mn$minrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=pred_2), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb)) {
      graph_B <- graph_B + annotate("path",x=knt_comb$xx[rr]+knt_comb$rad[rr]*cos(seq(0,2*pi,length.out=100)),
                                    y=knt_comb$yy[rr]+knt_comb$rad[rr]*sin(seq(0,2*pi,length.out=100))) +
        geom_text(aes(x=xx, y=yy, label=BB), data=lab_pos)
    }
  }
  graph_C <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class C ", row_mn$noise, " noise case, ", row_mn$model, ", ", row_mn$noPoints, " datapoints, simulation no ", row_mn$minrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=pred_3), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb)) {
      graph_C <- graph_C + annotate("path",x=knt_comb$xx[rr]+knt_comb$rad[rr]*cos(seq(0,2*pi,length.out=100)),
                                    y=knt_comb$yy[rr]+knt_comb$rad[rr]*sin(seq(0,2*pi,length.out=100)))
    }
  }
  graph_biasA <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class A ", row_mn$noise, " noise case, ", row_mn$model, ", ", row_mn$noPoints, " datapoints, simulation no ", row_mn$minrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=bias_1), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb)) {
      graph_biasA <- graph_biasA + annotate("path",x=knt_comb$xx[rr]+knt_comb$rad[rr]*cos(seq(0,2*pi,length.out=100)),
                                    y=knt_comb$yy[rr]+knt_comb$rad[rr]*sin(seq(0,2*pi,length.out=100)), colour='white') 
    }
  }
  
  graph_biasB <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class B ", row_mn$noise, " noise case, ", row_mn$model, ", ", row_mn$noPoints, " datapoints, simulation no ", row_mn$minrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=bias_2), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb)) {
      graph_biasB <- graph_biasB + annotate("path",x=knt_comb$xx[rr]+knt_comb$rad[rr]*cos(seq(0,2*pi,length.out=100)),
                                            y=knt_comb$yy[rr]+knt_comb$rad[rr]*sin(seq(0,2*pi,length.out=100)), colour='white') 
    }
  }
  graph_biasC <- ggplot(bias_comb, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class C ", row_mn$noise, " noise case, ", row_mn$model, ", ", row_mn$noPoints, " datapoints, simulation no ", row_mn$minrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=bias_3), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb)) {
      graph_biasC <- graph_biasC + annotate("path",x=knt_comb$xx[rr]+knt_comb$rad[rr]*cos(seq(0,2*pi,length.out=100)),
                                            y=knt_comb$yy[rr]+knt_comb$rad[rr]*sin(seq(0,2*pi,length.out=100)), colour='white') 
    }
  }
  
  graph_A_mx <- ggplot(bias_comb_mx, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class A ", row_mx$noise, " noise case, ", row_mx$model, ", ", row_mx$noPoints, " datapoints, simulation no ", row_mx$maxrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=pred_1), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb_mx)) {
      graph_A_mx <- graph_A_mx + annotate("path", 
        x=knt_comb_mx$xx[rr]+knt_comb_mx$rad[rr]*cos(seq(0,2*pi,length.out=100)),
        y=knt_comb_mx$yy[rr]+knt_comb_mx$rad[rr]*sin(seq(0,2*pi,length.out=100))) +
        geom_text(aes(x=xx, y=yy, label=AA), data=lab_pos_mx)
    }
  }
  graph_B_mx <- ggplot(bias_comb_mx, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class B ", row_mx$noise, " noise case, ", row_mx$model, ", ", row_mx$noPoints, " datapoints, simulation no ", row_mx$maxrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=pred_2), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb_mx)) {
      graph_B_mx <- graph_B_mx + annotate("path",
        x=knt_comb_mx$xx[rr]+knt_comb_mx$rad[rr]*cos(seq(0,2*pi,length.out=100)),
        y=knt_comb_mx$yy[rr]+knt_comb_mx$rad[rr]*sin(seq(0,2*pi,length.out=100))) +
        geom_text(aes(x=xx, y=yy, label=BB), data=lab_pos_mx)
    }
  }
  graph_C_mx <- ggplot(bias_comb_mx, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class C ", row_mx$noise, " noise case, ", row_mx$model, ", ", row_mx$noPoints, " datapoints, simulation no ", row_mx$maxrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=pred_3), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb_mx)) {
      graph_C_mx <- graph_C_mx + annotate("path",
        x=knt_comb_mx$xx[rr]+knt_comb_mx$rad[rr]*cos(seq(0,2*pi,length.out=100)),
        y=knt_comb_mx$yy[rr]+knt_comb_mx$rad[rr]*sin(seq(0,2*pi,length.out=100)))
    }
  }
  graph_biasA_mx <- ggplot(bias_comb_mx, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class A ", row_mx$noise, " noise case, ", row_mx$model, ", ", row_mx$noPoints, " datapoints, simulation no ", row_mx$maxrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=bias_1), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb_mx)) {
      graph_biasA_mx <- graph_biasA_mx + annotate("path",
        x=knt_comb_mx$xx[rr]+knt_comb_mx$rad[rr]*cos(seq(0,2*pi,length.out=100)),
        y=knt_comb_mx$yy[rr]+knt_comb_mx$rad[rr]*sin(seq(0,2*pi,length.out=100)), colour='white') 
    }
  }
  graph_biasB_mx <- ggplot(bias_comb_mx, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class B ", row_mx$noise, " noise case, ", row_mx$model, ", ", row_mx$noPoints, " datapoints, simulation no ", row_mx$maxrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=bias_2), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb_mx)) {
      graph_biasB_mx <- graph_biasB_mx + annotate("path",
        x=knt_comb_mx$xx[rr]+knt_comb_mx$rad[rr]*cos(seq(0,2*pi,length.out=100)),
        y=knt_comb_mx$yy[rr]+knt_comb_mx$rad[rr]*sin(seq(0,2*pi,length.out=100)), colour='white') 
    }
  }
  graph_biasC_mx <- ggplot(bias_comb_mx, aes(xx, yy)) +
    theme_bw() + theme(legend.key=element_blank()) +
    coord_fixed(ratio = 1) +
    xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
    ggtitle(paste0("Bias for class C ", row_mx$noise, " noise case, ", row_mx$model, ", ", row_mx$noPoints, " datapoints, simulation no ", row_mx$maxrow)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    stat_summary_2d(aes(z=bias_3), fun=max, binwidth = c(0.02, 0.02)) +
    scale_fill_distiller(palette = "RdYlBu", direction=-1, name = "") +
    theme(legend.position="right")
  if (row_mn$model == 'cress') {
    for (rr in 1:nrow(knt_comb_mx)) {
      graph_biasC_mx <- graph_biasC_mx + annotate("path",
        x=knt_comb_mx$xx[rr]+knt_comb_mx$rad[rr]*cos(seq(0,2*pi,length.out=100)),
        y=knt_comb_mx$yy[rr]+knt_comb_mx$rad[rr]*sin(seq(0,2*pi,length.out=100)), colour='white') 
    }
  }
  
  pts_out_name <- paste0("Pt9Cl", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, ".png")
  pts_out_name_mx <- paste0("Pt9Cl", row_mx$model, row_mx$dataset, row_mx$noise, row_mx$noPoints, "_", row_mx$maxrow, "_best.png")
  tru_out_nam1 <- paste0("Pt9TruP", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "A.png")
  tru_out_nam2 <- paste0("Pt9TruP", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "B.png")
  tru_out_nam3 <- paste0("Pt9TruP", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "C.png")
  pred_out_nam1 <- paste0("Pt9Pred", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "A.png")
  pred_out_nam2 <- paste0("Pt9Pred", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "B.png")
  pred_out_nam3 <- paste0("Pt9Pred", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "C.png")
  bias_out_nam1 <- paste0("Pt9Bias", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "A.png")
  bias_out_nam2 <- paste0("Pt9Bias", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "B.png")
  bias_out_nam3 <- paste0("Pt9Bias", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mn$minrow, "C.png")
  pred_out_nam1_mx <- paste0("Pt9Pred", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mx$maxrow, "A_best.png")
  pred_out_nam2_mx <- paste0("Pt9Pred", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mx$maxrow, "B_best.png")
  pred_out_nam3_mx <- paste0("Pt9Pred", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mx$maxrow, "C_best.png")
  bias_out_nam1_mx <- paste0("Pt9Bias", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mx$maxrow, "A_best.png")
  bias_out_nam2_mx <- paste0("Pt9Bias", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mx$maxrow, "B_best.png")
  bias_out_nam3_mx <- paste0("Pt9Bias", row_mn$model, row_mn$dataset, row_mn$noise, row_mn$noPoints, "_", row_mx$maxrow, "C_best.png")
  
  png(filename=paste0(img_save_dir, pts_out_name))
  print(pred_pts)
  dev.off()
  png(filename=paste0(img_save_dir, tru_out_nam1))
  print(true_A)
  dev.off()
  png(filename=paste0(img_save_dir, tru_out_nam2))
  print(true_B)
  dev.off()
  png(filename=paste0(img_save_dir, tru_out_nam3))
  print(true_C)
  dev.off()
  png(filename=paste0(img_save_dir, pred_out_nam1))
  print(graph_A)
  dev.off()
  png(filename=paste0(img_save_dir, pred_out_nam2))
  print(graph_B)
  dev.off()
  png(filename=paste0(img_save_dir, pred_out_nam3))
  print(graph_C)
  dev.off()
  png(filename=paste0(img_save_dir, bias_out_nam1))
  print(graph_biasA)
  dev.off()
  png(filename=paste0(img_save_dir, bias_out_nam2))
  print(graph_biasB)
  dev.off()
  png(filename=paste0(img_save_dir, bias_out_nam3))
  print(graph_biasC)
  dev.off()
  
  png(filename=paste0(img_save_dir, pts_out_name_mx))
  print(pred_pts_mx)
  dev.off()
  png(filename=paste0(img_save_dir, pred_out_nam1_mx))
  print(graph_A_mx)
  dev.off()
  png(filename=paste0(img_save_dir, pred_out_nam2_mx))
  print(graph_B_mx)
  dev.off()
  png(filename=paste0(img_save_dir, pred_out_nam3_mx))
  print(graph_C_mx)
  dev.off()
  png(filename=paste0(img_save_dir, bias_out_nam1_mx))
  print(graph_biasA_mx)
  dev.off()
  png(filename=paste0(img_save_dir, bias_out_nam2_mx))
  print(graph_biasB_mx)
  dev.off()
  png(filename=paste0(img_save_dir, bias_out_nam3_mx))
  print(graph_biasC_mx)
  dev.off()
}


