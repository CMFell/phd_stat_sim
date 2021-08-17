# This file fits SVM based models to the simulated data
# This is in 3 sections for simple, complex and island separated by lines of ########
# WARNING - try running each section separately each may take several hours to run
# set the directory locations in lines 6-13 then run whole file.

# base directory to save all output
base_dir_on_comp <- "C:/Users/christina/"
base_dir_sim <- base_dir_on_comp + "multinomial_simulations/"
# subdirectories to save csv output and generated images
file_save_dir <- paste0(base_dir_sim, "csv_files/")
img_save_dir <- paste0(base_dir_sim, "image_files/")

library(fields)
library(VGAM)
library(sf)
library(Rfast)
library(knitr)
library(tidyverse)
library(mgcv)
library(kableExtra)
library(e1071)

str.colnames <- function(basen, nof){
  names.out <- c()
  for (st in 1:nof) {
    names.out <- c(names.out, paste0(basen, "_", st))
  }
  return(names.out)
}

acc_func <- function(pred, truth){
  init_table = table(pred, truth)
  dim_init_tab = dim(init_table)
  if (dim_init_tab[1]==dim_init_tab[2]){
    acc_out <- sum(diag(init_table)) / length(truth)
  } else {
    colz = colnames(init_table)
    rowz = rownames(init_table)
    if (dim_init_tab[2] > dim_init_tab[1]){
      for (dd in 1:dim_init_tab[2]) {
        if (!(colz[dd] %in% rowz)) {
          rws_before = init_table[0:(dd-1), ]
          if (dd > nrow(init_table)) {
            new_rw = rep(0, dim_init_tab[2])
            init_table <- rbind(rws_before, new_rw)
            rownames(init_table)[dd] <- colz[dd]
          } else {
            rws_after = init_table[dd:nrow(init_table),]
            new_rw = rep(0, dim_init_tab[2])
            init_table_new <- rbind(rws_before, new_rw, rws_after)
            rownames(init_table_new)[0:(dd-1)] <- rownames(init_table)[0:(dd-1)]
            rownames(init_table_new)[dd] <- colz[dd]
            rownames(init_table_new)[(dd+1):nrow(init_table_new)] <- rownames(init_table)[dd:nrow(init_table)]
            init_table <- init_table_new
          }
        }
      }
    } else {
      for (dd in 1:dim_init_tab[1]) {
        if (!(rowz[dd] %in% colz)) {
          cls_before = init_table[, 0:(dd-1)]
          if (dd > ncol(init_table)) {
            new_cl = rep(0, dim_init_tab[1])
            init_table <- cbind(cls_before, new_cl)
            colnames(init_table)[dd] <- rowz[dd]
          } else {
            cls_after = init_table[,dd:ncol(init_table)]
            new_cl = rep(0, dim_init_tab[1])
            init_table_new <- cbind(cls_before, new_cl, cls_after)
            colnames(init_table_new)[0:(dd-1)] <- colnames(init_table)[0:(dd-1)]
            colnames(init_table_new)[dd] <- rowz[dd]
            colnames(init_table_new)[(dd+1):ncol(init_table_new)] <- colnames(init_table)[dd:ncol(init_table)]
            init_table <- init_table_new
          }
        }
      }
    }
    acc_out <- sum(diag(init_table)) / length(truth)
  }
  
  rez <- list("accuracy"=acc_out, "conf_mat"=init_table)
  
  return(rez)
}

svm_func <- function(all_data, samp_data, name_label) {
  
  resp_colno <- which(colnames(all_data)==paste0("datacl_", name_label))
  resp_colno_samp <- which(colnames(samp_data)=="response")
  true_colno <- which(colnames(all_data)==paste0("truthcl_", name_label))
  
  datmat <- samp_data[, c(1,2,resp_colno_samp)]
  predmat <- all_data[, c(1,2,true_colno)]
  colnames(predmat)[3] <- "response"
  datmat$response <- as.factor(datmat$response)
  predmat$response <- as.factor(predmat$response)
  
  svm_start <- Sys.time()
  #create svm
  svmfit <- svm(response~., data=datmat, kernel ="radial", gamma=1, cost=1, probability=TRUE)
  svm_time <- Sys.time() - svm_start
  # Get predicted probabilities, classes and store
  fit_svm <- predict(svmfit, datmat, probability=TRUE)
  pred_svm <- predict(svmfit, predmat, probability=TRUE)
  pred_p_svm <- as.data.frame(attr(pred_svm, 'probabilities'), row.names=NULL)
  
  classez <- c("A", "B", "C")
  fit_cl_svm <- as.vector(fit_svm)
  pred_cl_svm <- as.vector(pred_svm)
  resid_true_svm <- pred_cl_svm == all_data[,true_colno]
  
  svm_out <- cbind(pred_p_svm$A, pred_p_svm$B, pred_p_svm$C)
  svm_out <- as.data.frame(svm_out)
  svm_out <- cbind(svm_out, pred_cl_svm, resid_true_svm)
  colnames(svm_out)[1:3] <- str.colnames(paste0("predp_svm_",name_label), 3)
  colnames(svm_out)[4] <- paste0("predcl_svm_", name_label)
  colnames(svm_out)[5] <- paste0("resid_true_svm_", name_label)
  
  acc_truth_svm_all <- acc_func(pred_cl_svm, all_data[,true_colno])
  acc_data_svm_all <- acc_func(pred_cl_svm, all_data[,resp_colno])
  acc_truth_svm_fit <- acc_func(fit_cl_svm, samp_data[,true_colno])
  acc_data_svm_fit <- acc_func(fit_cl_svm, samp_data[,resp_colno])
  
  conf_mat <- cbind(acc_truth_svm_all$conf_mat, acc_data_svm_all$conf_mat, 
                    acc_truth_svm_fit$conf_mat, acc_data_svm_fit$conf_mat)
  
  dat_svm <- as.data.frame(matrix(c(svm_time, acc_truth_svm_all$accuracy, acc_data_svm_all$accuracy, acc_truth_svm_fit$accuracy, acc_data_svm_fit$accuracy), ncol=5))
  colnames(dat_svm) <- c(
    paste0("time_svm_", name_label), 
    paste0("acc_truth_svm_all_", name_label),
    paste0("acc_data_svm_all_", name_label),
    paste0("acc_truth_svm_fit_", name_label),
    paste0("acc_data_svm_fit_", name_label)
  )
  
  list_svm <- list("output_mat"=svm_out, "output_dat"=dat_svm, "model"=svm_out, "conf_mat"=conf_mat)
  
  return(list_svm)
}


file_in_dir <- file_save_dir


if (file.exists(paste0(file_in_dir, "dataset_rows1.csv"))) {
  dataset_rows1 <- read.csv(paste0(file_in_dir, "dataset_rows1.csv"))
  dataset_rows2 <- read.csv(paste0(file_in_dir, "dataset_rows2.csv"))
  dataset_rows3 <- read.csv(paste0(file_in_dir, "dataset_rows3.csv"))
  dataset_rows1i <- read.csv(paste0(file_in_dir, "dataset_rows1_island.csv"))
  dataset_rows2i <- read.csv(paste0(file_in_dir, "dataset_rows2_island.csv"))
  dataset_rows3i <- read.csv(paste0(file_in_dir, "dataset_rows3_island.csv"))
}

if (file.exists(paste0(file_in_dir, "location_grid.csv"))) { 
  grid <- read.csv(paste0(file_in_dir, "location_grid.csv"))
} 
if (file.exists(paste0(file_in_dir, "location_grid_island.csv"))) { 
  grid_i <- read.csv(paste0(file_in_dir, "location_grid_island.csv"))
}



data_in = "simp"
input_name <- paste0(file_in_dir, data_in, "_input_data.csv")
inputdata_class_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_data.csv")
inputdata_numeric_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_datanum.csv")
inputdata_class_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_data.csv")
inputdata_numeric_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_datanum.csv")
outmat_name <- paste0(file_save_dir, "svm_", data_in, "_mat2.csv")
outstr_name <- paste0(file_save_dir, "svm_", data_in, "_dat2.csv")
confmat_name <- paste0(file_save_dir, "svm_", data_in, "_conf2.csv")
valz_name <- paste0(file_save_dir, "svm_", data_in, "_valz2.csv")
listout_name <- paste0(file_save_dir, "svm_", data_in, "_deets_list2.csv")
matout_name <- paste0(file_save_dir, "svm_", data_in, "_deets_mat2.csv")
coefout_name <- paste0(file_save_dir, "svm_", data_in, "_coefs.csv")

store_state <- read.csv(paste0(file_in_dir, "store_state_mgis.csv"))

data_out <- read.csv(input_name)
data_lowF <- read.csv(inputdata_class_lw)
data_highF <- read.csv(inputdata_class_hi)
data_lowN <- read.csv(inputdata_numeric_lw)
data_highN <- read.csv(inputdata_numeric_hi)

repz = 100

# Set model values
fitmeas <- "BIC"
sksk <- 6
mnk <- 2
mxk <- 20
gpgp <- 0
knkn <- 300

data_out$x.pos <- data_out$xx
data_out$y.pos <- data_out$yy

for (ii in 1:repz) {
  print(ii)
  sim_data <- cbind(data_out[,c("xx", "yy", "x.pos", "y.pos", "truthcl_simp_lw", "truthcl_simp_hi")], 
                    data_lowF[,ii], data_highF[,ii], data_lowN[,ii], data_highN[,ii])
  colnames(sim_data) <- c("xx", "yy", "x.pos", "y.pos", "truthcl_simp_lw", "truthcl_simp_hi", 
                          "datacl_simp_lw", "datacl_simp_hi", "dataclnum_simp_lw", "dataclnum_simp_hi")
  
  samp_rows1 <- unlist(dataset_rows1[ii,])
  samp_rows2 <- unlist(dataset_rows2[ii,])
  samp_rows3 <- unlist(dataset_rows3[ii,])
  samp_data1 <- sim_data[samp_rows1,]
  samp_data2 <- sim_data[samp_rows2,]
  samp_data3 <- sim_data[samp_rows3,]
  nsamp1 <- length(samp_rows1)
  nsamp2 <- length(samp_rows2)
  nsamp3 <- length(samp_rows3)
  
  # data
  dat_cl_naml <- paste0("datacl_", data_in, "_lw")
  dat_cl_namh <- paste0("datacl_", data_in, "_hi")
  samp_data1$response <- samp_data1[,dat_cl_naml]
  samp_data <- samp_data1
  print("low1")
  svm_low1 <- svm_func(sim_data, samp_data, paste0(data_in, "_lw"))
  samp_data2$response <- samp_data2[,dat_cl_naml]
  samp_data <- samp_data2
  print("low2")
  svm_low2 <- svm_func(sim_data, samp_data, paste0(data_in, "_lw"))
  samp_data3$response <- samp_data3[,dat_cl_naml]
  samp_data <- samp_data3
  print("low3")
  svm_low3 <- svm_func(sim_data, samp_data, paste0(data_in, "_lw"))
  samp_data1$response <- samp_data1[,dat_cl_namh]
  samp_data <- samp_data1
  print("high1")
  svm_hi1 <- svm_func(sim_data, samp_data, paste0(data_in, "_hi"))
  samp_data2$response <- samp_data2[,dat_cl_namh]
  samp_data <- samp_data2
  print("high2")
  svm_hi2 <- svm_func(sim_data, samp_data, paste0(data_in, "_hi"))
  samp_data3$response <- samp_data3[,dat_cl_namh]
  samp_data <- samp_data3
  print("high3")
  svm_hi3 <- svm_func(sim_data, samp_data, paste0(data_in, "_hi"))
  
  # write out results 
  samp_out1 <- cbind(rep(paste0(ii, "_1"), 10000), svm_low1$output_mat, svm_hi1$output_mat)
  colnames(samp_out1)[1] <- c("sim_no")
  samp_out2 <- cbind(rep(paste0(ii, "_2"), 10000), svm_low2$output_mat, svm_hi2$output_mat)
  colnames(samp_out2)[1] <- c("sim_no")
  samp_out3 <- cbind(rep(paste0(ii, "_3"), 10000), svm_low3$output_mat, svm_hi3$output_mat)
  colnames(samp_out3)[1] <- c("sim_no")
  samp_out <- rbind(samp_out1, samp_out2, samp_out3)
  sim_res1 <- cbind(paste0(ii, "_1"), svm_low1$output_dat, svm_hi1$output_dat)
  colnames(sim_res1)[1] <- "sim_no"
  sim_res2 <- cbind(paste0(ii, "_2"), svm_low2$output_dat, svm_hi2$output_dat)
  colnames(sim_res2)[1] <- "sim_no"
  sim_res3 <- cbind(paste0(ii, "_3"), svm_low3$output_dat, svm_hi3$output_dat)
  colnames(sim_res3)[1] <- "sim_no"
  sim_res <- rbind(sim_res1, sim_res2, sim_res3)
  
  models_list <- list(svm_low1$model, svm_hi1$model,
                      svm_low2$model, svm_hi2$model,
                      svm_low3$model, svm_hi3$model)
  names(models_list) <- c(paste0("svmlow1_", ii), paste0("svmhigh1_", ii),
                          paste0("svmlow2_", ii), paste0("svmhigh2_", ii),
                          paste0("svmlow3_", ii), paste0("svmhigh3_", ii))
  
  conf_mat <- rbind(svm_low1$conf_mat, svm_low2$conf_mat, svm_low3$conf_mat,
                    svm_hi1$conf_mat, svm_hi2$conf_mat, svm_hi3$conf_mat)
  
  if (ii == 1) {
    output_matrix <- samp_out
    output_store <- sim_res
    models_out <- models_list
    conf_mat_svm <- conf_mat
  } else {
    output_matrix <- rbind(output_matrix, samp_out)
    output_store <- rbind(output_store, sim_res)
    models_out <- append(models_out, models_list)
    conf_mat_svm <- conf_mat_svm + conf_mat
  }
}

write.csv(output_matrix, outmat_name, row.names=FALSE)
write.csv(output_store, outstr_name, row.names=FALSE)
write.csv(conf_mat_svm, confmat_name, row.names=FALSE)
valz <- c(repz, fitmeas, sksk, mnk, mxk, gpgp, knkn)
write.csv(valz, valz_name, row.names=FALSE)
write.csv(store_state, paste0(file_in_dir, "store_state_svsi.csv"), row.names=FALSE)

out_deets <- list()

for (md in 1:length(models_out)){
  nam <- names(models_out)[md]
  mod <- models_out[[nam]]
  coefz <- mod$coefficients
  if (md==1){
    namez <- nam
    coefz_out <- coefz
  } else {
    namez <- c(namez, nam)
    coefz_out <- c(coefz_out, coefz)
  }
}

write.csv(coefz_out, coefout_name, row.names=FALSE)


#######################################################################


data_in <- "comp"
input_name <- paste0(file_save_dir, data_in, "_input_data.csv")
inputdata_class_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_data.csv")
inputdata_numeric_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_datanum.csv")
inputdata_class_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_data.csv")
inputdata_numeric_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_datanum.csv")
outmat_name <- paste0(file_save_dir, "svm_", data_in, "_mat2.csv")
outstr_name <- paste0(file_save_dir, "svm_", data_in, "_dat2.csv")
confmat_name <- paste0(file_save_dir, "svm_", data_in, "_conf2.csv")
valz_name <- paste0(file_save_dir, "svm_", data_in, "_valz2.csv")
listout_name <- paste0(file_save_dir, "svm_", data_in, "_deets_list2.csv")
matout_name <- paste0(file_save_dir, "svm_", data_in, "_deets_mat2.csv")
coefout_name <- paste0(file_save_dir, "svm_", data_in, "_coefs.csv")

store_state <- read.csv(paste0(file_in_dir, "store_state_svsi.csv"))

data_out <- read.csv(input_name)
data_lowF <- read.csv(inputdata_class_lw)
data_highF <- read.csv(inputdata_class_hi)
data_lowN <- read.csv(inputdata_numeric_lw)
data_highN <- read.csv(inputdata_numeric_hi)

repz = 100

# Set model values
fitmeas <- "BIC"
sksk <- 6
mnk <- 2
mxk <- 20
gpgp <- 0
knkn <- 300

data_out$x.pos <- data_out$xx
data_out$y.pos <- data_out$yy

for (ii in 1:repz) {
  print(ii)
  sim_data <- cbind(data_out[,c("xx", "yy", "x.pos", "y.pos", "truthcl_comp_lw", "truthcl_comp_hi")], 
                    data_lowF[,ii], data_highF[,ii], data_lowN[,ii], data_highN[,ii])
  colnames(sim_data) <- c("xx", "yy", "x.pos", "y.pos", "truthcl_comp_lw", "truthcl_comp_hi", 
                          "datacl_comp_lw", "datacl_comp_hi", "dataclnum_comp_lw", "dataclnum_comp_hi")
  
  samp_rows1 <- unlist(dataset_rows1[ii,])
  samp_rows2 <- unlist(dataset_rows2[ii,])
  samp_rows3 <- unlist(dataset_rows3[ii,])
  samp_data1 <- sim_data[samp_rows1,]
  samp_data2 <- sim_data[samp_rows2,]
  samp_data3 <- sim_data[samp_rows3,]
  nsamp1 <- length(samp_rows1)
  nsamp2 <- length(samp_rows2)
  nsamp3 <- length(samp_rows3)
  
  # data
  dat_cl_naml <- paste0("datacl_", data_in, "_lw")
  dat_cl_naml <- "datacl_comp_lw"
  dat_cl_namh <- paste0("datacl_", data_in, "_hi")
  samp_data1$response <- samp_data1[,dat_cl_naml]
  samp_data <- samp_data1
  svm_low1 <- svm_func(sim_data, samp_data, paste0(data_in, "_lw"))
  samp_data2$response <- samp_data2[,dat_cl_naml]
  samp_data <- samp_data2
  svm_low2 <- svm_func(sim_data, samp_data, paste0(data_in, "_lw"))
  samp_data3$response <- samp_data3[,dat_cl_naml]
  samp_data <- samp_data3
  svm_low3 <- svm_func(sim_data, samp_data, paste0(data_in, "_lw"))
  samp_data1$response <- samp_data1[,dat_cl_namh]
  samp_data <- samp_data1
  svm_hi1 <- svm_func(sim_data, samp_data, paste0(data_in, "_hi"))
  samp_data2$response <- samp_data2[,dat_cl_namh]
  samp_data <- samp_data2
  svm_hi2 <- svm_func(sim_data, samp_data, paste0(data_in, "_hi"))
  samp_data3$response <- samp_data3[,dat_cl_namh]
  samp_data <- samp_data3
  svm_hi3 <- svm_func(sim_data, samp_data, paste0(data_in, "_hi"))
  
  # write out results 
  samp_out1 <- cbind(rep(paste0(ii, "_1"), 10000), svm_low1$output_mat, svm_hi1$output_mat)
  colnames(samp_out1)[1] <- c("sim_no")
  samp_out2 <- cbind(rep(paste0(ii, "_2"), 10000), svm_low2$output_mat, svm_hi2$output_mat)
  colnames(samp_out2)[1] <- c("sim_no")
  samp_out3 <- cbind(rep(paste0(ii, "_3"), 10000), svm_low3$output_mat, svm_hi3$output_mat)
  colnames(samp_out3)[1] <- c("sim_no")
  samp_out <- rbind(samp_out1, samp_out2, samp_out3)
  sim_res1 <- cbind(paste0(ii, "_1"), svm_low1$output_dat, svm_hi1$output_dat)
  colnames(sim_res1)[1] <- "sim_no"
  sim_res2 <- cbind(paste0(ii, "_2"), svm_low2$output_dat, svm_hi2$output_dat)
  colnames(sim_res2)[1] <- "sim_no"
  sim_res3 <- cbind(paste0(ii, "_3"), svm_low3$output_dat, svm_hi3$output_dat)
  colnames(sim_res3)[1] <- "sim_no"
  sim_res <- rbind(sim_res1, sim_res2, sim_res3)
  
  models_list <- list(svm_low1$model, svm_hi1$model,
                      svm_low2$model, svm_hi2$model,
                      svm_low3$model, svm_hi3$model)
  names(models_list) <- c(paste0("svmlow1_", ii), paste0("svmhigh1_", ii),
                          paste0("svmlow2_", ii), paste0("svmhigh2_", ii),
                          paste0("svmlow3_", ii), paste0("svmhigh3_", ii))
  
  conf_mat <- rbind(svm_low1$conf_mat, svm_low2$conf_mat, svm_low3$conf_mat,
                    svm_hi1$conf_mat, svm_hi2$conf_mat, svm_hi3$conf_mat)
  
  if (ii == 1) {
    output_matrix <- samp_out
    output_store <- sim_res
    models_out <- models_list
    conf_mat_svm <- conf_mat
  } else {
    output_matrix <- rbind(output_matrix, samp_out)
    output_store <- rbind(output_store, sim_res)
    models_out <- append(models_out, models_list)
    conf_mat_svm <- conf_mat_svm + conf_mat
  }
}

write.csv(output_matrix, outmat_name, row.names=FALSE)
write.csv(output_store, outstr_name, row.names=FALSE)
write.csv(conf_mat_svm, confmat_name, row.names=FALSE)
valz <- c(repz, fitmeas, sksk, mnk, mxk, gpgp, knkn)
write.csv(valz, valz_name, row.names=FALSE)
write.csv(store_state, paste0(file_in_dir, "store_state_svco.csv"), row.names=FALSE)

out_deets <- list()

for (md in 1:length(models_out)){
  nam <- names(models_out)[md]
  mod <- models_out[[nam]]
  coefz <- mod$coefficients
  if (md==1){
    namez <- nam
    coefz_out <- coefz
  } else {
    namez <- c(namez, nam)
    coefz_out <- c(coefz_out, coefz)
  }
}

write.csv(coefz_out, coefout_name, row.names=FALSE)


#######################################################################


data_in = "island"
input_name <- paste0(file_in_dir, data_in, "_input_data.csv")
inputdata_class_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_data.csv")
inputdata_numeric_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_datanum.csv")
inputdata_class_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_data.csv")
inputdata_numeric_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_datanum.csv")
outmat_name <- paste0(file_save_dir, "svm_", data_in, "_mat2_new.csv")
outstr_name <- paste0(file_save_dir, "svm_", data_in, "_dat2_new.csv")
confmat_name <- paste0(file_save_dir, "svm_", data_in, "_conf2_new.csv")
valz_name <- paste0(file_save_dir, "svm_", data_in, "_valz2_new.csv")
listout_name <- paste0(file_save_dir, "svm_", data_in, "_deets_list2_new.csv")
matout_name <- paste0(file_save_dir, "svm_", data_in, "_deets_mat2_new.csv")
dist_name <- paste0(file_in_dir, data_in, "_geodistmat.csv")
coefout_name <- paste0(file_save_dir, "svm_", data_in, "_coefs_new.csv")

store_state <- read.csv(paste0(file_in_dir, "store_state_svco.csv"))

data_out <- read.csv(input_name)
data_lowF <- read.csv(inputdata_class_lw)
data_highF <- read.csv(inputdata_class_hi)
data_lowN <- read.csv(inputdata_numeric_lw)
data_highN <- read.csv(inputdata_numeric_hi)

if (!exists(dist_name)) {
  GeoDistMat <- read.csv(dist_name)
}

GeoDistMat <- as.matrix(GeoDistMat)

repz = 100

# Set model values
fitmeas <- "BIC"
sksk <- 6
mnk <- 2
mxk <- 20
gpgp <- 0
knkn <- 300

data_out$x.pos <- data_out$xx
data_out$y.pos <- data_out$yy

for (ii in 1:repz) {
  print(paste0("rep", ii))
  sim_data <- cbind(data_out[,c("xx", "yy", "x.pos", "y.pos", "truthcl_island_lw", "truthcl_island_hi")], 
                    data_lowF[,ii], data_highF[,ii], data_lowN[,ii], data_highN[,ii])
  colnames(sim_data) <- c("xx", "yy", "x.pos", "y.pos", "truthcl_island_lw", "truthcl_island_hi", 
                          "datacl_island_lw", "datacl_island_hi", "dataclnum_island_lw", "dataclnum_island_hi")
  
  samp_rows1 <- unlist(dataset_rows1i[ii,])
  samp_rows2 <- unlist(dataset_rows2i[ii,])
  samp_rows3 <- unlist(dataset_rows3i[ii,])
  samp_data1 <- sim_data[samp_rows1,]
  samp_data2 <- sim_data[samp_rows2,]
  samp_data3 <- sim_data[samp_rows3,]
  nsamp1 <- length(samp_rows1)
  nsamp2 <- length(samp_rows2)
  nsamp3 <- length(samp_rows3)
  GeoDistMat_samp1 <- GeoDistMat[samp_rows1,]
  GeoDistMat_samp2 <- GeoDistMat[samp_rows2,]
  GeoDistMat_samp3 <- GeoDistMat[samp_rows3,]
  
  # data
  dat_cl_naml <- paste0("datacl_", data_in, "_lw")
  dat_cl_namh <- paste0("datacl_", data_in, "_hi")
  print("lw1")
  samp_data1$response <- samp_data1[,dat_cl_naml]
  samp_data <- samp_data1
  svm_low1 <- suppressWarnings(svm_func(sim_data, samp_data, paste0(data_in, "_lw")))
  print("lw2")
  samp_data2$response <- samp_data2[,dat_cl_naml]
  samp_data <- samp_data2
  svm_low2 <- suppressWarnings(svm_func(sim_data, samp_data, paste0(data_in, "_lw")))
  print("lw3")
  samp_data3$response <- samp_data3[,dat_cl_naml]
  samp_data <- samp_data3
  svm_low3 <- suppressWarnings(svm_func(sim_data, samp_data3, paste0(data_in, "_lw")))
  print("hi1")
  samp_data1$response <- samp_data1[,dat_cl_namh]
  samp_data <- samp_data1
  svm_hi1 <- suppressWarnings(svm_func(sim_data, samp_data1, paste0(data_in, "_hi")))
  print("hi2")
  samp_data2$response <- samp_data2[,dat_cl_namh]
  samp_data <- samp_data2
  svm_hi2 <- suppressWarnings(svm_func(sim_data, samp_data2, paste0(data_in, "_hi")))
  print("hi3")
  samp_data3$response <- samp_data3[,dat_cl_namh]
  samp_data <- samp_data3
  svm_hi3 <- suppressWarnings(svm_func(sim_data, samp_data3, paste0(data_in, "_hi")))
  
  # write out results 
  samp_out1 <- cbind(rep(paste0(ii, "_1"), nrow(data_out)), svm_low1$output_mat, svm_hi1$output_mat)
  colnames(samp_out1)[1] <- c("sim_no")
  samp_out2 <- cbind(rep(paste0(ii, "_2"), nrow(data_out)), svm_low2$output_mat, svm_hi2$output_mat)
  colnames(samp_out2)[1] <- c("sim_no")
  samp_out3 <- cbind(rep(paste0(ii, "_3"), nrow(data_out)), svm_low3$output_mat, svm_hi3$output_mat)
  colnames(samp_out3)[1] <- c("sim_no")
  samp_out <- rbind(samp_out1, samp_out2, samp_out3)
  sim_res1 <- cbind(paste0(ii, "_1"), svm_low1$output_dat, svm_hi1$output_dat)
  colnames(sim_res1)[1] <- "sim_no"
  sim_res2 <- cbind(paste0(ii, "_2"), svm_low2$output_dat, svm_hi2$output_dat)
  colnames(sim_res2)[1] <- "sim_no"
  sim_res3 <- cbind(paste0(ii, "_3"), svm_low3$output_dat, svm_hi3$output_dat)
  colnames(sim_res3)[1] <- "sim_no"
  sim_res <- rbind(sim_res1, sim_res2, sim_res3)
  
  models_list <- list(svm_low1$model, svm_hi1$model,
                      svm_low2$model, svm_hi2$model,
                      svm_low3$model, svm_hi3$model)
  names(models_list) <- c(paste0("svmlow1_", ii), paste0("svmhigh1_", ii),
                          paste0("svmlow2_", ii), paste0("svmhigh2_", ii),
                          paste0("svmlow3_", ii), paste0("svmhigh3_", ii))
  
  conf_mat <- rbind(svm_low1$conf_mat, svm_low2$conf_mat, svm_low3$conf_mat,
                    svm_hi1$conf_mat, svm_hi2$conf_mat, svm_hi3$conf_mat)
  
  if (ii == 1) {
    output_matrix <- samp_out
    output_store <- sim_res
    models_out <- models_list
    conf_mat_svm <- conf_mat
  } else {
    output_matrix <- rbind(output_matrix, samp_out)
    output_store <- rbind(output_store, sim_res)
    models_out <- append(models_out, models_list)
    conf_mat_svm <- conf_mat_svm + conf_mat
  }
}

write.csv(output_matrix, outmat_name, row.names=FALSE)
write.csv(output_store, outstr_name, row.names=FALSE)
write.csv(conf_mat_svm, confmat_name, row.names=FALSE)
valz <- c(repz, fitmeas, sksk, mnk, mxk, gpgp, knkn)
write.csv(valz, valz_name, row.names=FALSE)
write.csv(store_state, paste0(file_in_dir, "store_state_svis.csv"), row.names=FALSE)

out_deets <- list()

for (md in 1:length(models_out)){
  nam <- names(models_out)[md]
  mod <- models_out[[nam]]
  coefz <- mod$coefficients
  if (md==1){
    namez <- nam
    coefz_out <- coefz
  } else {
    namez <- c(namez, nam)
    coefz_out <- c(coefz_out, coefz)
  }
}

write.csv(coefz_out, coefout_name, row.names=FALSE)

