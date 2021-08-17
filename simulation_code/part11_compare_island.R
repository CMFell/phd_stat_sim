# This file collates results writen out from CRESS, MGCV and SVM models on the island data
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


cbbPalette <- c("#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")

grid_i <- read.csv(paste0(file_save_dir, "location_grid_island.csv"))
low_s2n_island <- read.csv(paste0(file_save_dir, "low_s2n_island.csv"))
high_s2n_island <- read.csv(paste0(file_save_dir, "high_s2n_island.csv"))

### CRESS

trup_col_namll <- paste0("truep_island_lw_1")
trup_col_naml2 <- paste0("truep_island_lw_2")
trup_col_naml3 <- paste0("truep_island_lw_3")
trup_col_namhl <- paste0("truep_island_hi_1")
trup_col_namh2 <- paste0("truep_island_hi_2")
trup_col_namh3 <- paste0("truep_island_hi_3")
predp_col_namll <- paste0("predp_cress_island_lw_1")
predp_col_naml2 <- paste0("predp_cress_island_lw_2")
predp_col_naml3 <- paste0("predp_cress_island_lw_3")
predp_col_namhl <- paste0("predp_cress_island_hi_1")
predp_col_namh2 <- paste0("predp_cress_island_hi_2")
predp_col_namh3 <- paste0("predp_cress_island_hi_3")

output_matrix <- read.csv(paste0(file_save_dir, "cress_island_mat2_new_cv.csv"), sep=",", header=TRUE)

output_cr_si_lw <- as.data.frame(matrix(rep(0, nrow(low_s2n_island)*12), ncol=12))
colnames(output_cr_si_lw) <- c("pA_lw_rez_1", "pA_lw_rez_2", "pA_lw_rez_3",
                               "pB_lw_rez_1", "pB_lw_rez_2", "pB_lw_rez_3",
                               "pC_lw_rez_1", "pC_lw_rez_2", "pC_lw_rez_3",
                               "cl_lw_rez_1", "cl_lw_rez_2", "cl_lw_rez_3")

output_cr_si_hi <- as.data.frame(matrix(rep(0, nrow(high_s2n_island)*12), ncol=12))
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
  pA_lw_rez_1 <- low_s2n_island[,trup_col_namll] - ex_1[,predp_col_namll]
  pA_lw_rez_2 <- low_s2n_island[,trup_col_namll] - ex_2[,predp_col_namll]
  pA_lw_rez_3 <- low_s2n_island[,trup_col_namll] - ex_3[,predp_col_namll]
  pB_lw_rez_1 <- low_s2n_island[,trup_col_naml2] - ex_1[,predp_col_naml2]
  pB_lw_rez_2 <- low_s2n_island[,trup_col_naml2] - ex_2[,predp_col_naml2]
  pB_lw_rez_3 <- low_s2n_island[,trup_col_naml2] - ex_3[,predp_col_naml2]
  pC_lw_rez_1 <- low_s2n_island[,trup_col_naml3] - ex_1[,predp_col_naml3]
  pC_lw_rez_2 <- low_s2n_island[,trup_col_naml3] - ex_2[,predp_col_naml3]
  pC_lw_rez_3 <- low_s2n_island[,trup_col_naml3] - ex_3[,predp_col_naml3]
  seA_lw_rez_1 <- abs(pA_lw_rez_1)
  seA_lw_rez_2 <- abs(pA_lw_rez_2)
  seA_lw_rez_3 <- abs(pA_lw_rez_3)
  seB_lw_rez_1 <- abs(pB_lw_rez_1)
  seB_lw_rez_2 <- abs(pB_lw_rez_2)
  seB_lw_rez_3 <- abs(pB_lw_rez_3)
  seC_lw_rez_1 <- abs(pC_lw_rez_1)
  seC_lw_rez_2 <- abs(pC_lw_rez_2)
  seC_lw_rez_3 <- abs(pC_lw_rez_3)
  cl_lw_rez_1 <- ex_1[,"resid_true_cress_island_lw"]
  cl_lw_rez_2 <- ex_2[,"resid_true_cress_island_lw"]
  cl_lw_rez_3 <- ex_3[,"resid_true_cress_island_lw"]
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
  
  pA_hi_rez_1 <- high_s2n_island[,trup_col_namhl] - ex_1[,predp_col_namhl]
  pA_hi_rez_2 <- high_s2n_island[,trup_col_namhl] - ex_2[,predp_col_namhl]
  pA_hi_rez_3 <- high_s2n_island[,trup_col_namhl] - ex_3[,predp_col_namhl]
  pB_hi_rez_1 <- high_s2n_island[,trup_col_namh2] - ex_1[,predp_col_namh2]
  pB_hi_rez_2 <- high_s2n_island[,trup_col_namh2] - ex_2[,predp_col_namh2]
  pB_hi_rez_3 <- high_s2n_island[,trup_col_namh2] - ex_3[,predp_col_namh2]
  pC_hi_rez_1 <- high_s2n_island[,trup_col_namh3] - ex_1[,predp_col_namh3]
  pC_hi_rez_2 <- high_s2n_island[,trup_col_namh3] - ex_2[,predp_col_namh3]
  pC_hi_rez_3 <- high_s2n_island[,trup_col_namh3] - ex_3[,predp_col_namh3]
  seA_hi_rez_1 <- abs(pA_hi_rez_1)
  seA_hi_rez_2 <- abs(pA_hi_rez_2)
  seA_hi_rez_3 <- abs(pA_hi_rez_3)
  seB_hi_rez_1 <- abs(pB_hi_rez_1)
  seB_hi_rez_2 <- abs(pB_hi_rez_2)
  seB_hi_rez_3 <- abs(pB_hi_rez_3)
  seC_hi_rez_1 <- abs(pC_hi_rez_1)
  seC_hi_rez_2 <- abs(pC_hi_rez_2)
  seC_hi_rez_3 <- abs(pC_hi_rez_3)
  cl_hi_rez_1 <- ex_1[,"resid_true_cress_island_hi"]
  cl_hi_rez_2 <- ex_2[,"resid_true_cress_island_hi"]
  cl_hi_rez_3 <- ex_3[,"resid_true_cress_island_hi"]
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
output_cr_si_lw <- cbind.data.frame(output_cr_si_lw, grid_i)
output_cr_si_hi <- cbind.data.frame(output_cr_si_hi, grid_i)

### MGCV

predp_col_namll <- paste0("predp_mgcv_island_lw_1")
predp_col_naml2 <- paste0("predp_mgcv_island_lw_2")
predp_col_naml3 <- paste0("predp_mgcv_island_lw_3")
predp_col_namhl <- paste0("predp_mgcv_island_hi_1")
predp_col_namh2 <- paste0("predp_mgcv_island_hi_2")
predp_col_namh3 <- paste0("predp_mgcv_island_hi_3")

output_matrix <- read.csv(paste0(file_save_dir, "mgcv_island_mat2_new.csv"), sep=",", header=TRUE)

output_mg_si_lw <- as.data.frame(matrix(rep(0, nrow(low_s2n_island)*12), ncol=12))
colnames(output_mg_si_lw) <- c("pA_lw_rez_1", "pA_lw_rez_2", "pA_lw_rez_3",
                               "pB_lw_rez_1", "pB_lw_rez_2", "pB_lw_rez_3",
                               "pC_lw_rez_1", "pC_lw_rez_2", "pC_lw_rez_3",
                               "cl_lw_rez_1", "cl_lw_rez_2", "cl_lw_rez_3")

output_mg_si_hi <- as.data.frame(matrix(rep(0, nrow(high_s2n_island)*12), ncol=12))
colnames(output_mg_si_hi) <- c("pA_hi_rez_1", "pA_hi_rez_2", "pA_hi_rez_3",
                               "pB_hi_rez_1", "pB_hi_rez_2", "pB_hi_rez_3",
                               "pC_hi_rez_1", "pC_hi_rez_2", "pC_hi_rez_3",
                               "cl_hi_rez_1", "cl_hi_rez_2", "cl_hi_rez_3")

rss_mg_si <- as.data.frame(matrix(rep(0, repz*7), ncol=7))
colnames(rss_mg_si) <- c("sim", "lw_1", "lw_2", "lw_3", "hi_1", "hi2", "hi3")
rss_mg_si$sim <- seq(1,100)

for (rp in 1:repz) {
  print(rp)
  ex_1 <- output_matrix[output_matrix$sim_no == paste0(rp, "_1"),]
  ex_2 <- output_matrix[output_matrix$sim_no == paste0(rp, "_2"),]
  ex_3 <- output_matrix[output_matrix$sim_no == paste0(rp, "_3"),]
  pA_lw_rez_1 <- low_s2n_island[,trup_col_namll] - ex_1[,predp_col_namll]
  pA_lw_rez_2 <- low_s2n_island[,trup_col_namll] - ex_2[,predp_col_namll]
  pA_lw_rez_3 <- low_s2n_island[,trup_col_namll] - ex_3[,predp_col_namll]
  pB_lw_rez_1 <- low_s2n_island[,trup_col_naml2] - ex_1[,predp_col_naml2]
  pB_lw_rez_2 <- low_s2n_island[,trup_col_naml2] - ex_2[,predp_col_naml2]
  pB_lw_rez_3 <- low_s2n_island[,trup_col_naml2] - ex_3[,predp_col_naml2]
  pC_lw_rez_1 <- low_s2n_island[,trup_col_naml3] - ex_1[,predp_col_naml3]
  pC_lw_rez_2 <- low_s2n_island[,trup_col_naml3] - ex_2[,predp_col_naml3]
  pC_lw_rez_3 <- low_s2n_island[,trup_col_naml3] - ex_3[,predp_col_naml3]
  seA_lw_rez_1 <- abs(pA_lw_rez_1)
  seA_lw_rez_2 <- abs(pA_lw_rez_2)
  seA_lw_rez_3 <- abs(pA_lw_rez_3)
  seB_lw_rez_1 <- abs(pB_lw_rez_1)
  seB_lw_rez_2 <- abs(pB_lw_rez_2)
  seB_lw_rez_3 <- abs(pB_lw_rez_3)
  seC_lw_rez_1 <- abs(pC_lw_rez_1)
  seC_lw_rez_2 <- abs(pC_lw_rez_2)
  seC_lw_rez_3 <- abs(pC_lw_rez_3)
  cl_lw_rez_1 <- ex_1[,"resid_true_mgcv_island_lw"]
  cl_lw_rez_2 <- ex_2[,"resid_true_mgcv_island_lw"]
  cl_lw_rez_3 <- ex_3[,"resid_true_mgcv_island_lw"]
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
  
  pA_hi_rez_1 <- high_s2n_island[,trup_col_namhl] - ex_1[,predp_col_namhl]
  pA_hi_rez_2 <- high_s2n_island[,trup_col_namhl] - ex_2[,predp_col_namhl]
  pA_hi_rez_3 <- high_s2n_island[,trup_col_namhl] - ex_3[,predp_col_namhl]
  pB_hi_rez_1 <- high_s2n_island[,trup_col_namh2] - ex_1[,predp_col_namh2]
  pB_hi_rez_2 <- high_s2n_island[,trup_col_namh2] - ex_2[,predp_col_namh2]
  pB_hi_rez_3 <- high_s2n_island[,trup_col_namh2] - ex_3[,predp_col_namh2]
  pC_hi_rez_1 <- high_s2n_island[,trup_col_namh3] - ex_1[,predp_col_namh3]
  pC_hi_rez_2 <- high_s2n_island[,trup_col_namh3] - ex_2[,predp_col_namh3]
  pC_hi_rez_3 <- high_s2n_island[,trup_col_namh3] - ex_3[,predp_col_namh3]
  seA_hi_rez_1 <- abs(pA_hi_rez_1)
  seA_hi_rez_2 <- abs(pA_hi_rez_2)
  seA_hi_rez_3 <- abs(pA_hi_rez_3)
  seB_hi_rez_1 <- abs(pB_hi_rez_1)
  seB_hi_rez_2 <- abs(pB_hi_rez_2)
  seB_hi_rez_3 <- abs(pB_hi_rez_3)
  seC_hi_rez_1 <- abs(pC_hi_rez_1)
  seC_hi_rez_2 <- abs(pC_hi_rez_2)
  seC_hi_rez_3 <- abs(pC_hi_rez_3)
  cl_hi_rez_1 <- ex_1[,"resid_true_mgcv_island_hi"]
  cl_hi_rez_2 <- ex_2[,"resid_true_mgcv_island_hi"]
  cl_hi_rez_3 <- ex_3[,"resid_true_mgcv_island_hi"]
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
output_mg_si_lw <- cbind.data.frame(output_mg_si_lw, grid_i)
output_mg_si_hi <- cbind.data.frame(output_mg_si_hi, grid_i)

### SVM

predp_col_namll <- paste0("predp_svm_island_lw_1")
predp_col_naml2 <- paste0("predp_svm_island_lw_2")
predp_col_naml3 <- paste0("predp_svm_island_lw_3")
predp_col_namhl <- paste0("predp_svm_island_hi_1")
predp_col_namh2 <- paste0("predp_svm_island_hi_2")
predp_col_namh3 <- paste0("predp_svm_island_hi_3")

output_matrix <- read.csv(paste0(file_save_dir, "svm_island_mat2_new.csv"), sep=",", header=TRUE)

output_sv_si_lw <- as.data.frame(matrix(rep(0, nrow(low_s2n_island)*12), ncol=12))
colnames(output_sv_si_lw) <- c("pA_lw_rez_1", "pA_lw_rez_2", "pA_lw_rez_3",
                               "pB_lw_rez_1", "pB_lw_rez_2", "pB_lw_rez_3",
                               "pC_lw_rez_1", "pC_lw_rez_2", "pC_lw_rez_3",
                               "cl_lw_rez_1", "cl_lw_rez_2", "cl_lw_rez_3")

output_sv_si_hi <- as.data.frame(matrix(rep(0, nrow(high_s2n_island)*12), ncol=12))
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
  pA_lw_rez_1 <- low_s2n_island[,trup_col_namll] - ex_1[,predp_col_namll]
  pA_lw_rez_2 <- low_s2n_island[,trup_col_namll] - ex_2[,predp_col_namll]
  pA_lw_rez_3 <- low_s2n_island[,trup_col_namll] - ex_3[,predp_col_namll]
  pB_lw_rez_1 <- low_s2n_island[,trup_col_naml2] - ex_1[,predp_col_naml2]
  pB_lw_rez_2 <- low_s2n_island[,trup_col_naml2] - ex_2[,predp_col_naml2]
  pB_lw_rez_3 <- low_s2n_island[,trup_col_naml2] - ex_3[,predp_col_naml2]
  pC_lw_rez_1 <- low_s2n_island[,trup_col_naml3] - ex_1[,predp_col_naml3]
  pC_lw_rez_2 <- low_s2n_island[,trup_col_naml3] - ex_2[,predp_col_naml3]
  pC_lw_rez_3 <- low_s2n_island[,trup_col_naml3] - ex_3[,predp_col_naml3]
  seA_lw_rez_1 <- abs(pA_lw_rez_1)
  seA_lw_rez_2 <- abs(pA_lw_rez_2)
  seA_lw_rez_3 <- abs(pA_lw_rez_3)
  seB_lw_rez_1 <- abs(pB_lw_rez_1)
  seB_lw_rez_2 <- abs(pB_lw_rez_2)
  seB_lw_rez_3 <- abs(pB_lw_rez_3)
  seC_lw_rez_1 <- abs(pC_lw_rez_1)
  seC_lw_rez_2 <- abs(pC_lw_rez_2)
  seC_lw_rez_3 <- abs(pC_lw_rez_3)
  cl_lw_rez_1 <- ex_1[,"resid_true_svm_island_lw"]
  cl_lw_rez_2 <- ex_2[,"resid_true_svm_island_lw"]
  cl_lw_rez_3 <- ex_3[,"resid_true_svm_island_lw"]
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
  
  pA_hi_rez_1 <- high_s2n_island[,trup_col_namhl] - ex_1[,predp_col_namhl]
  pA_hi_rez_2 <- high_s2n_island[,trup_col_namhl] - ex_2[,predp_col_namhl]
  pA_hi_rez_3 <- high_s2n_island[,trup_col_namhl] - ex_3[,predp_col_namhl]
  pB_hi_rez_1 <- high_s2n_island[,trup_col_namh2] - ex_1[,predp_col_namh2]
  pB_hi_rez_2 <- high_s2n_island[,trup_col_namh2] - ex_2[,predp_col_namh2]
  pB_hi_rez_3 <- high_s2n_island[,trup_col_namh2] - ex_3[,predp_col_namh2]
  pC_hi_rez_1 <- high_s2n_island[,trup_col_namh3] - ex_1[,predp_col_namh3]
  pC_hi_rez_2 <- high_s2n_island[,trup_col_namh3] - ex_2[,predp_col_namh3]
  pC_hi_rez_3 <- high_s2n_island[,trup_col_namh3] - ex_3[,predp_col_namh3]
  seA_hi_rez_1 <- abs(pA_hi_rez_1)
  seA_hi_rez_2 <- abs(pA_hi_rez_2)
  seA_hi_rez_3 <- abs(pA_hi_rez_3)
  seB_hi_rez_1 <- abs(pB_hi_rez_1)
  seB_hi_rez_2 <- abs(pB_hi_rez_2)
  seB_hi_rez_3 <- abs(pB_hi_rez_3)
  seC_hi_rez_1 <- abs(pC_hi_rez_1)
  seC_hi_rez_2 <- abs(pC_hi_rez_2)
  seC_hi_rez_3 <- abs(pC_hi_rez_3)
  cl_hi_rez_1 <- ex_1[,"resid_true_svm_island_hi"]
  cl_hi_rez_2 <- ex_2[,"resid_true_svm_island_hi"]
  cl_hi_rez_3 <- ex_3[,"resid_true_svm_island_hi"]
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
output_sv_si_lw <- cbind.data.frame(output_sv_si_lw, grid_i)
output_sv_si_hi <- cbind.data.frame(output_sv_si_hi, grid_i)


write.csv(rss_cr_si, paste0(img_save_dir, "rss_cress_island.csv"))
write.csv(rss_mg_si, paste0(img_save_dir, "rss_mgcv_island.csv"))
write.csv(rss_sv_si, paste0(img_save_dir, "rss_svm_island.csv"))

library(scales)


graph_cress_A_island_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_island_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_island_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_island_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_island_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_island_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_island_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_island_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_island_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_island_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_island_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_island_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_island_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_island_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_island_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_island_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_island_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_island_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_island_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_island_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_island_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_island_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_island_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_island_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_island_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_island_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_island_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))



### High noise


graph_cress_A_island_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_island_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_island_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_island_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_island_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_island_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_A_island_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_B_island_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_C_island_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_island_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_island_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_island_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_island_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_island_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_island_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_A_island_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_B_island_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_C_island_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))


graph_svm_A_island_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_island_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_island_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_island_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_island_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_island_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_A_island_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pA_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_B_island_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pB_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_C_island_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=pC_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n square \n error") +
  theme(legend.position="right", text = element_text(size=20))





### Classes

graph_cress_cl_island_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_island_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_island_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_island_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_island_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_island_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_island_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_island_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_island_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_lw_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_island_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_island_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_cl_island_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_island_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_island_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_cl_island_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_island_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_island_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=cl_hi_rez_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=-1, limits=c(0,1), oob=squish, name=" proportion \n of correct \n class \n predictions") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_cl_island_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
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

graph_cress_rss_island_lw1 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_island_lw2 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_island_lw3 <- ggplot(output_cr_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_island_lw1 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_island_lw2 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_island_lw3 <- ggplot(output_mg_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_island_lw1 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_island_lw2 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_island_lw3 <- ggplot(output_sv_si_lw, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=lw_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.25), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_island_hi1 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_island_hi2 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_cress_rss_island_hi3 <- ggplot(output_cr_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_island_hi1 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_island_hi2 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_mgcv_rss_island_hi3 <- ggplot(output_mg_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_island_hi1 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_1), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_island_hi2 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_2), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))

graph_svm_rss_island_hi3 <- ggplot(output_sv_si_hi, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=hi_rss_3), fun=max, binwidth = c(0.01, 0.01)) +
  scale_fill_distiller(palette = "YlGnBu", direction=1, limits=c(0,0.1), oob=squish, name=" probability \n mean \n error") +
  theme(legend.position="right", text = element_text(size=20))






png(filename=paste0(img_save_dir, "GraphCrAIslandLw1.png"))
graph_cress_A_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBIslandLw1.png"))
graph_cress_B_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCIslandLw1.png"))
graph_cress_C_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrAIslandLw2.png"))
graph_cress_A_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBIslandLw2.png"))
graph_cress_B_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCIslandLw2.png"))
graph_cress_C_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrAIslandLw3.png"))
graph_cress_A_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBIslandLw3.png"))
graph_cress_B_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCIslandLw3.png"))
graph_cress_C_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrAIslandHi1.png"))
graph_cress_A_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBIslandHi1.png"))
graph_cress_B_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCIslandHi1.png"))
graph_cress_C_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrAIslandHi2.png"))
graph_cress_A_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBIslandHi2.png"))
graph_cress_B_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCIslandHi2.png"))
graph_cress_C_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrAIslandHi3.png"))
graph_cress_A_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrBIslandHi3.png"))
graph_cress_B_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrCIslandHi3.png"))
graph_cress_C_island_hi3
dev.off()

png(filename=paste0(img_save_dir, "GraphMgAIslandLw1.png"))
graph_mgcv_A_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBIslandLw1.png"))
graph_mgcv_B_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCIslandLw1.png"))
graph_mgcv_C_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgAIslandLw2.png"))
graph_mgcv_A_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBIslandLw2.png"))
graph_mgcv_B_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCIslandLw2.png"))
graph_mgcv_C_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgAIslandLw3.png"))
graph_mgcv_A_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBIslandLw3.png"))
graph_mgcv_B_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCIslandLw3.png"))
graph_mgcv_C_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgAIslandHi1.png"))
graph_mgcv_A_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBIslandHi1.png"))
graph_mgcv_B_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCIslandHi1.png"))
graph_mgcv_C_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgAIslandHi2.png"))
graph_mgcv_A_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBIslandHi2.png"))
graph_mgcv_B_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCIslandHi2.png"))
graph_mgcv_C_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgAIslandHi3.png"))
graph_mgcv_A_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgBIslandHi3.png"))
graph_mgcv_B_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgCIslandHi3.png"))
graph_mgcv_C_island_hi3
dev.off()


png(filename=paste0(img_save_dir, "GraphSvAIslandLw1.png"))
graph_svm_A_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBIslandLw1.png"))
graph_svm_B_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCIslandLw1.png"))
graph_svm_C_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvAIslandLw2.png"))
graph_svm_A_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBIslandLw2.png"))
graph_svm_B_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCIslandLw2.png"))
graph_svm_C_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvAIslandLw3.png"))
graph_svm_A_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBIslandLw3.png"))
graph_svm_B_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCIslandLw3.png"))
graph_svm_C_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvAIslandHi1.png"))
graph_svm_A_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBIslandHi1.png"))
graph_svm_B_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCIslandHi1.png"))
graph_svm_C_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvAIslandHi2.png"))
graph_svm_A_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBIslandHi2.png"))
graph_svm_B_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCIslandHi2.png"))
graph_svm_C_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvAIslandHi3.png"))
graph_svm_A_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvBIslandHi3.png"))
graph_svm_B_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvCIslandHi3.png"))
graph_svm_C_island_hi3
dev.off()








png(filename=paste0(img_save_dir, "GraphCrClIslandLw1.png"))
graph_cress_cl_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClIslandLw2.png"))
graph_cress_cl_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClIslandLw3.png"))
graph_cress_cl_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClIslandLw1.png"))
graph_mgcv_cl_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClIslandLw2.png"))
graph_mgcv_cl_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClIslandLw3.png"))
graph_mgcv_cl_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClIslandLw1.png"))
graph_svm_cl_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClIslandLw2.png"))
graph_svm_cl_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClIslandLw3.png"))
graph_svm_cl_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClIslandHi1.png"))
graph_cress_cl_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClIslandHi2.png"))
graph_cress_cl_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrClIslandHi3.png"))
graph_cress_cl_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClIslandHi1.png"))
graph_mgcv_cl_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClIslandHi2.png"))
graph_mgcv_cl_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgClIslandHi3.png"))
graph_mgcv_cl_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClIslandHi1.png"))
graph_svm_cl_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClIslandHi2.png"))
graph_svm_cl_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvClIslandHi3.png"))
graph_svm_cl_island_hi3
dev.off()



png(filename=paste0(img_save_dir, "GraphCrRssIslandLw1.png"))
graph_cress_rss_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssIslandLw2.png"))
graph_cress_rss_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssIslandLw3.png"))
graph_cress_rss_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssIslandLw1.png"))
graph_mgcv_rss_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssIslandLw2.png"))
graph_mgcv_rss_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssIslandLw3.png"))
graph_mgcv_rss_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssIslandLw1.png"))
graph_svm_rss_island_lw1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssIslandLw2.png"))
graph_svm_rss_island_lw2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssIslandLw3.png"))
graph_svm_rss_island_lw3
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssIslandHi1.png"))
graph_cress_rss_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssIslandHi2.png"))
graph_cress_rss_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphCrRssIslandHi3.png"))
graph_cress_rss_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssIslandHi1.png"))
graph_mgcv_rss_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssIslandHi2.png"))
graph_mgcv_rss_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphMgRssIslandHi3.png"))
graph_mgcv_rss_island_hi3
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssIslandHi1.png"))
graph_svm_rss_island_hi1
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssIslandHi2.png"))
graph_svm_rss_island_hi2
dev.off()
png(filename=paste0(img_save_dir, "GraphSvRssIslandHi3.png"))
graph_svm_rss_island_hi3
dev.off()


