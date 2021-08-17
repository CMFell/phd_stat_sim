# This file fits CReSS based models to the simulated data
# This is in 3 sections for simple, complex and island separated by lines of ########
# WARNING - try running each section separately each may take more than 24 hours to run
# set the directory locations in lines 6-13 then run whole file.

# base directory to save all output
base_dir_on_comp <- "C:/Users/christina/"
base_dir_sim <- base_dir_on_comp + "multinomial_simulations/"
# subdirectories to save csv output and generated images
file_save_dir <- paste0(base_dir_sim, "csv_files/")
img_save_dir <- paste0(base_dir_sim, "image_files/")
# code currently in multinomial branch of MRSea so need to download package from github to base dir on comp 
# then switch to multinomial branch
devtools::load_all(paste0(base_dir_on_comp, "MRsea"))

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

#sim_data, samp_data, paste0(data_in, "_lw"), d2k_samp1, salsa2dlist1, k2k_samp1, d2k_mat1

cress_func <- function(all_data, samp_data, name_label, sampd2k, salsalist, k2k, all_dist, chserad=FALSE) {
  
  resp_colno <- which(colnames(all_data)==paste0("datacl_", name_label))
  true_colno <- which(colnames(all_data)==paste0("truthcl_", name_label))
  
  cress_start <- Sys.time()
  # Create intercept only model
  mod_int_only <- vglm(
    response~1, 
    family=multinomial, 
    data=samp_data
  )
  # Fit multinomial model
  mod_mn_sim <- runSALSA2Dmn(
    mod_int_only, 
    salsalist,
    vglmdatain=samp_data,
    d2k=sampd2k, 
    k2k=k2k,
    suppress.printout=TRUE,
    chooserad=chserad,
    hdetest=TRUE
  )
  cress_time <- Sys.time() - cress_start
  # Get predicted probabilities, classes and store
  fit_p_cress <- fitted(mod_mn_sim$bestModel)
  if (class(all_dist) == "matrix") {
    all_dist_mat = all_dist
  } else {
    all_dist_mat = all_dist$dataDist
  }
  pred_p_cress_all <- predict.vglmMRSea(newdata=all_data[,c("x.pos","y.pos")], object=mod_mn_sim$bestModel, newdists=all_dist_mat, type="response", conf_int=T)
  pred_p_cress <- as.data.frame(pred_p_cress_all$predictions, row.names=NULL)
  pred_p_cress_lw <- as.data.frame(pred_p_cress_all$lower_limit, row.names=NULL)
  pred_p_cress_hi <- as.data.frame(pred_p_cress_all$higher_limit, row.names=NULL)
  ci_comb <- cbind.data.frame(pred_p_cress_lw, pred_p_cress_hi, pred_p_cress)
  colnames(ci_comb) <- c("lowA", "lowB", "lowC", "highA", "highB", "highC", "meanA", "meanB", "meanC")
  
  classez <- c("A", "B", "C")
  fit_cl_cress <- classez[apply(fit_p_cress,1,which.max)]
  pred_cl_cress <- classez[apply(pred_p_cress,1,which.max)]
  resid_true_cress <- pred_cl_cress == all_data[,true_colno]
  
  cress_out <- cbind(pred_p_cress, pred_cl_cress, resid_true_cress)
  colnames(cress_out)[1:3] <- str.colnames(paste0("predp_cress_",name_label), 3)
  colnames(cress_out)[4] <- paste0("predcl_cress_", name_label)
  colnames(cress_out)[5] <- paste0("resid_true_cress_", name_label)
  
  acc_truth_cress_all <- acc_func(pred_cl_cress, all_data[,true_colno])
  acc_data_cress_all <- acc_func(pred_cl_cress, all_data[,resp_colno])
  acc_truth_cress_fit <- acc_func(fit_cl_cress, samp_data[,true_colno])
  acc_data_cress_fit <- acc_func(fit_cl_cress, samp_data[,resp_colno])
  
  conf_mat <- cbind(acc_truth_cress_all$conf_mat, acc_data_cress_all$conf_mat, 
                    acc_truth_cress_fit$conf_mat, acc_data_cress_fit$conf_mat)
  
  dat_cress <- as.data.frame(matrix(c(cress_time, acc_truth_cress_all$accuracy, acc_data_cress_all$accuracy, acc_truth_cress_fit$accuracy, acc_data_cress_fit$accuracy), ncol=5))
  colnames(dat_cress) <- c(
    paste0("time_cress_", name_label), 
    paste0("acc_truth_cress_all_", name_label),
    paste0("acc_data_cress_all_", name_label),
    paste0("acc_truth_cress_fit_", name_label),
    paste0("acc_data_cress_fit_", name_label)
  )
  
  list_cress <- list("output_mat"=cress_out, "output_dat"=dat_cress, "model"=mod_mn_sim, "conf_mat"=conf_mat, "ci_out"=ci_comb)
  
  return(list_cress)
}


data_in = "simp"
input_name <- paste0(file_in_dir, data_in, "_input_data.csv")
inputdata_class_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_data.csv")
inputdata_numeric_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_datanum.csv")
inputdata_class_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_data.csv")
inputdata_numeric_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_datanum.csv")
outmat_name <- paste0(file_save_dir, "cress_", data_in, "_mat2_new_cv.csv")
outstr_name <- paste0(file_save_dir, "cress_", data_in, "_dat2_new_cv.csv")
confmat_name <- paste0(file_save_dir, "cress_", data_in, "_conf2_new_cv.csv")
valz_name <- paste0(file_save_dir, "cress_", data_in, "_valz2_new_cv.csv")
listout_name <- paste0(file_save_dir, "cress_", data_in, "_deets_list2_new_cv.csv")
matout_name <- paste0(file_save_dir, "cress_", data_in, "_deets_mat2_new_cv.csv")
coefout_name <- paste0(file_save_dir, "cress_", data_in, "_coefs_new_cv.csv")
coefstd_name <- paste0(file_save_dir, "cress_", data_in, "_coefstd_new_cv.csv")
cov_name <- paste0(file_save_dir, "cress_", data_in, "_cov_new_cv.csv")
ci_name <- paste0(file_save_dir, "cress_", data_in, "_cipred_new_cv.csv")

if (file.exists(paste0(file_in_dir, "store_state_island.csv"))) { 
  store_state <- read.csv(paste0(file_in_dir, "store_state_island.csv"))
}

data_out <- read.csv(input_name)
data_lowF <- read.csv(inputdata_class_lw)
data_highF <- read.csv(inputdata_class_hi)
data_lowN <- read.csv(inputdata_numeric_lw)
data_highN <- read.csv(inputdata_numeric_hi)

repz = 100

# Set model values
fitmeas <- "cv.gamMRSea"
sksk <- 6
mnk <- 2
mxk <- 20
gpgp <- 0
knkn <- 300

data_out$x.pos <- data_out$xx
data_out$y.pos <- data_out$yy

store_state <- cbind(store_state, .Random.seed)
colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_1_1")

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
  
  knotgrid1 <- cbind(samp_data1$xx, samp_data1$yy)
  samp_dst1 <- makeDists(
    cbind(samp_data1$xx, samp_data1$yy),
    na.omit(knotgrid1)
  )
  all_dst1 <- makeDists(
    cbind(data_out$xx, data_out$yy),
    na.omit(knotgrid1)
  )
  #store_state <- cbind(store_state, .Random.seed)
  #colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_2_", ii)
  knotgrid2 <- getKnotgrid(coordData = cbind(samp_data2$xx, samp_data2$yy), numKnots = knkn)
  samp_dst2 <- makeDists(
    cbind(samp_data2$xx, samp_data2$yy),
    na.omit(knotgrid2)
  )
  all_dst2 <- makeDists(
    cbind(data_out$xx, data_out$yy),
    na.omit(knotgrid2)
  )
  #store_state <- cbind(store_state, .Random.seed)
  #colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_3_", ii)
  knotgrid3 <- getKnotgrid(coordData = cbind(samp_data3$xx, samp_data3$yy), numKnots = knkn)
  samp_dst3 <- makeDists(
    cbind(samp_data3$xx, samp_data3$yy),
    na.omit(knotgrid3)
  )
  all_dst3 <- makeDists(
    cbind(data_out$xx, data_out$yy),
    na.omit(knotgrid3)
  )
  # create salsa input for all runs
  knotgrid=knotgrid1
  salsa2dlist1 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid1, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10
  )
  knotgrid=knotgrid2
  salsa2dlist2 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid2, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10
  )
  knotgrid=knotgrid3
  salsa2dlist3 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid3, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10
  )
  
  # data
  dat_cl_naml <- paste0("datacl_", data_in, "_lw")
  dat_cl_namh <- paste0("datacl_", data_in, "_hi")
  samp_data1$response <- samp_data1[,dat_cl_naml]
  samp_data <- samp_data1
  print("simp low 1")
  cress_low1 <- cress_func(sim_data, samp_data, paste0(data_in, "_lw"), samp_dst1$dataDist, salsa2dlist1, samp_dst1$knotDist, all_dst1)
  samp_data2$response <- samp_data2[,dat_cl_naml]
  samp_data <- samp_data2
  print("simp low 2")
  cress_low2 <- cress_func(sim_data, samp_data, paste0(data_in, "_lw"), samp_dst2$dataDist, salsa2dlist2, samp_dst2$knotDist, all_dst2)
  samp_data3$response <- samp_data3[,dat_cl_naml]
  samp_data <- samp_data3
  print("simp low 3")
  cress_low3 <- cress_func(sim_data, samp_data, paste0(data_in, "_lw"), samp_dst3$dataDist, salsa2dlist3, samp_dst3$knotDist, all_dst3)
  samp_data1$response <- samp_data1[,dat_cl_namh]
  samp_data <- samp_data1
  print("simp high 1")
  cress_hi1 <- cress_func(sim_data, samp_data, paste0(data_in, "_hi"), samp_dst1$dataDist, salsa2dlist1, samp_dst1$knotDist, all_dst1)
  samp_data2$response <- samp_data2[,dat_cl_namh]
  samp_data <- samp_data2
  print("simp high 2")
  cress_hi2 <- cress_func(sim_data, samp_data, paste0(data_in, "_hi"), samp_dst2$dataDist, salsa2dlist2, samp_dst2$knotDist, all_dst2)
  samp_data3$response <- samp_data3[,dat_cl_namh]
  samp_data <- samp_data3
  print("simp high 3")
  cress_hi3 <- cress_func(sim_data, samp_data, paste0(data_in, "_hi"), samp_dst3$dataDist, salsa2dlist3, samp_dst3$knotDist, all_dst3)
  
  # write out results 
  samp_out1 <- cbind(rep(paste0(ii, "_1"), 10000), cress_low1$output_mat, cress_hi1$output_mat)
  colnames(samp_out1)[1] <- c("sim_no")
  samp_out2 <- cbind(rep(paste0(ii, "_2"), 10000),cress_low2$output_mat, cress_hi2$output_mat)
  colnames(samp_out2)[1] <- c("sim_no")
  samp_out3 <- cbind(rep(paste0(ii, "_3"), 10000),cress_low3$output_mat, cress_hi3$output_mat)
  colnames(samp_out3)[1] <- c("sim_no")
  samp_out <- rbind(samp_out1, samp_out2, samp_out3)
  sim_res1 <- cbind(paste0(ii, "_1"), cress_low1$output_dat, cress_hi1$output_dat)
  colnames(sim_res1)[1] <- "sim_no"
  sim_res2 <- cbind(paste0(ii, "_2"), cress_low2$output_dat, cress_hi2$output_dat)
  colnames(sim_res2)[1] <- "sim_no"
  sim_res3 <- cbind(paste0(ii, "_3"), cress_low3$output_dat, cress_hi3$output_dat)
  colnames(sim_res3)[1] <- "sim_no"
  sim_res <- rbind(sim_res1, sim_res2, sim_res3)
  ci_low1 <- cbind.data.frame(rep(paste0("cresslow1_", ii), nrow(cress_low1$ci_out)), cress_low1$ci_out)
  colnames(ci_low1)[1] <- "sim_no"
  ci_low2 <- cbind.data.frame(rep(paste0("cresslow2_", ii), nrow(cress_low2$ci_out)), cress_low2$ci_out)
  colnames(ci_low2)[1] <- "sim_no"
  ci_low3 <- cbind.data.frame(rep(paste0("cresslow3_", ii), nrow(cress_low3$ci_out)), cress_low3$ci_out)
  colnames(ci_low3)[1] <- "sim_no"
  ci_high1 <- cbind.data.frame(rep(paste0("cresshigh1_", ii), nrow(cress_hi1$ci_out)), cress_hi1$ci_out)
  colnames(ci_high1)[1] <- "sim_no"
  ci_high2 <- cbind.data.frame(rep(paste0("cresshigh2_", ii), nrow(cress_hi2$ci_out)), cress_hi2$ci_out)
  colnames(ci_high2)[1] <- "sim_no"
  ci_high3 <- cbind.data.frame(rep(paste0("cresshigh3_", ii), nrow(cress_hi3$ci_out)), cress_hi3$ci_out)
  colnames(ci_high3)[1] <- "sim_no"
  ci_comb <- rbind.data.frame(ci_low1, ci_low2, ci_low3, ci_high1, ci_high2, ci_high3)
  
  models_list <- list(cress_low1$model, cress_hi1$model,
                      cress_low2$model, cress_hi2$model,
                      cress_low3$model, cress_hi3$model)
  names(models_list) <- c(paste0("cresslow1_", ii), paste0("cresshigh1_", ii),
                          paste0("cresslow2_", ii), paste0("cresshigh2_", ii),
                          paste0("cresslow3_", ii), paste0("cresshigh3_", ii))
  
  conf_mat <- rbind(cress_low1$conf_mat, cress_low2$conf_mat, cress_low3$conf_mat,
                    cress_hi1$conf_mat, cress_hi2$conf_mat, cress_hi3$conf_mat)
  
  if (ii == 1) {
    output_matrix <- samp_out
    output_store <- sim_res
    models_out <- models_list
    knotgrid1list <- list(knotgrid1)
    names(knotgrid1list) <- "kgrid1_1"
    knotgrid2list <- list(knotgrid2)
    names(knotgrid2list) <- "kgrid2_1"
    knotgrid3list <- list(knotgrid3)
    names(knotgrid3list) <- "kgrid3_1"
    conf_mat_cress <- conf_mat
    ci_out <- ci_comb
  } else {
    output_matrix <- rbind(output_matrix, samp_out)
    output_store <- rbind(output_store, sim_res)
    models_out <- append(models_out, models_list)
    knotgrid1list <- append(knotgrid1list, knotgrid1)
    names(knotgrid1list)[ii] <- paste0("kgrid1_", ii)
    knotgrid2list <- append(knotgrid2list, knotgrid2)
    names(knotgrid2list)[ii] <- paste0("kgrid2_", ii)
    knotgrid3list <- append(knotgrid3list, knotgrid3)
    names(knotgrid3list)[ii] <- paste0("kgrid3_", ii)
    conf_mat_cress <- conf_mat_cress + conf_mat
    ci_out <- rbind.data.frame(ci_out, ci_comb)
  }
}

write.csv(output_matrix, outmat_name, row.names=FALSE)
write.csv(output_store, outstr_name, row.names=FALSE)
write.csv(conf_mat_cress, confmat_name, row.names=FALSE)
valz <- c(repz, fitmeas, sksk, mnk, mxk, gpgp, knkn)
write.csv(valz, valz_name, row.names=FALSE)
write.csv(store_state, paste0(file_in_dir, "store_state_crsi.csv"), row.names=FALSE)
write.csv(ci_out, ci_name, row.names=FALSE)

kgrid1_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid1.csv")
kgrid2_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid2.csv")
kgrid3_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid3.csv")

write.csv(knotgrid1list, kgrid1_name, row.names=FALSE)
write.csv(knotgrid2list, kgrid2_name, row.names=FALSE)
write.csv(knotgrid3list, kgrid3_name, row.names=FALSE)


out_deets <- list()

for (md in 1:length(models_out)){
  nam <- names(models_out)[md]
  mod <- models_out[[nam]]
  knotpoz <- mod$splineParams[[1]]$knotPos
  knotz <- mod$splineParams[[1]]$knotgrid[knotpoz,]
  radz <- mod$splineParams[[1]]$radii[mod$splineParams[[1]]$radiusIndices]
  nknotz <- length(knotpoz)
  coefz <- mod$bestModel@coefficients
  list_out <- list("knots"=knotz, "no_knots"=nknotz, "radii"=radz)
  out_deets[[nam]] <- list_out
  summary_out <- summary(mod$bestModel)
  coefstd <- summary_out@coef3[,2]
  covz <- c(summary_out@cov.unscaled)
  if (md==1){
    namez <- nam
    nknt_out <- nknotz
    knotz_out <- knotz
    radz_out <- radz
    coefz_out <- coefz
    coefstd_out <- coefstd
    covz_out <- covz
  } else {
    namez <- c(namez, nam)
    nknt_out <- c(nknt_out, nknotz)
    knotz_out <- rbind(knotz_out, knotz)
    radz_out <- c(radz_out, radz)
    coefz_out <- c(coefz_out, coefz)
    coefstd_out <- c(coefstd_out, coefstd)
    covz_out <- c(covz_out, covz)
  }
}

listout <- cbind(namez, nknt_out)
matout <- cbind(knotz_out, radz_out)
colnames(matout) <- c("Xknot", "Yknot", "radknot")
write.csv(listout, listout_name, row.names=FALSE)
write.csv(matout, matout_name, row.names=FALSE)
write.csv(coefz_out, coefout_name, row.names=FALSE)
write.csv(coefstd_out, coefstd_name, row.names=FALSE)
write.csv(covz_out, cov_name, row.names=FALSE)


#######################################################################


data_in <- "comp"
input_name <- paste0(file_save_dir, data_in, "_input_data.csv")
inputdata_class_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_data.csv")
inputdata_numeric_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_datanum.csv")
inputdata_class_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_data.csv")
inputdata_numeric_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_datanum.csv")
outmat_name <- paste0(file_save_dir, "cress_", data_in, "_mat2_new_cv.csv")
outstr_name <- paste0(file_save_dir, "cress_", data_in, "_dat2_new_cv.csv")
confmat_name <- paste0(file_save_dir, "cress_", data_in, "_conf2_new_cv.csv")
valz_name <- paste0(file_save_dir, "cress_", data_in, "_valz2_new_cv.csv")
listout_name <- paste0(file_save_dir, "cress_", data_in, "_deets_list2_new_cv.csv")
matout_name <- paste0(file_save_dir, "cress_", data_in, "_deets_mat2_new_cv.csv")
coefout_name <- paste0(file_save_dir, "cress_", data_in, "_coefs_new_cv.csv")
coefstd_name <- paste0(file_save_dir, "cress_", data_in, "_coefstd_new_cv.csv")
cov_name <- paste0(file_save_dir, "cress_", data_in, "_cov_new_cv.csv")
ci_name <- paste0(file_save_dir, "cress_", data_in, "_cipred_new_cv.csv")

if (file.exists(paste0(file_in_dir, "store_state_crsi.csv"))) { 
  store_state <- read.csv(paste0(file_in_dir, "store_state_crsi.csv"))
}

data_out <- read.csv(input_name)
data_lowF <- read.csv(inputdata_class_lw)
data_highF <- read.csv(inputdata_class_hi)
data_lowN <- read.csv(inputdata_numeric_lw)
data_highN <- read.csv(inputdata_numeric_hi)

repz = 100

# Set model values
fitmeas <- "cv.gamMRSea"
sksk <- 6
mnk <- 2
mxk <- 20
gpgp <- 0
knkn <- 300

data_out$x.pos <- data_out$xx
data_out$y.pos <- data_out$yy

store_state <- cbind(store_state, .Random.seed)
colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_1_1")

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
  
  knotgrid1 <- cbind(samp_data1$xx, samp_data1$yy)
  samp_dst1 <- makeDists(
    cbind(samp_data1$xx, samp_data1$yy),
    na.omit(knotgrid1)
  )
  all_dst1 <- makeDists(
    cbind(data_out$xx, data_out$yy),
    na.omit(knotgrid1)
  )
  #store_state <- cbind(store_state, .Random.seed)
  #colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_2_", ii)
  knotgrid2 <- getKnotgrid(coordData = cbind(samp_data2$xx, samp_data2$yy), numKnots = knkn)
  samp_dst2 <- makeDists(
    cbind(samp_data2$xx, samp_data2$yy),
    na.omit(knotgrid2)
  )
  all_dst2 <- makeDists(
    cbind(data_out$xx, data_out$yy),
    na.omit(knotgrid2)
  )
  #store_state <- cbind(store_state, .Random.seed)
  #colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_3_", ii)
  knotgrid3 <- getKnotgrid(coordData = cbind(samp_data3$xx, samp_data3$yy), numKnots = knkn)
  samp_dst3 <- makeDists(
    cbind(samp_data3$xx, samp_data3$yy),
    na.omit(knotgrid3)
  )
  all_dst3 <- makeDists(
    cbind(data_out$xx, data_out$yy),
    na.omit(knotgrid3)
  )
  # create salsa input for all runs
  knotgrid=knotgrid1
  salsa2dlist1 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid1, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10
  )
  knotgrid=knotgrid2
  salsa2dlist2 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid2, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10
  )
  knotgrid=knotgrid3
  salsa2dlist3 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid3, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10
  )
  
  # data
  dat_cl_naml <- paste0("datacl_", data_in, "_lw")
  dat_cl_namh <- paste0("datacl_", data_in, "_hi")
  samp_data1$response <- samp_data1[,dat_cl_naml]
  samp_data <- samp_data1
  print("comp low 1")
  cress_low1 <- cress_func(sim_data, samp_data, paste0(data_in, "_lw"), samp_dst1$dataDist, salsa2dlist1, samp_dst1$knotDist, all_dst1)
  samp_data2$response <- samp_data2[,dat_cl_naml]
  samp_data <- samp_data2
  print("comp low 2")
  cress_low2 <- cress_func(sim_data, samp_data, paste0(data_in, "_lw"), samp_dst2$dataDist, salsa2dlist2, samp_dst2$knotDist, all_dst2)
  samp_data3$response <- samp_data3[,dat_cl_naml]
  samp_data <- samp_data3
  print("comp low 3")
  cress_low3 <- cress_func(sim_data, samp_data, paste0(data_in, "_lw"), samp_dst3$dataDist, salsa2dlist3, samp_dst3$knotDist, all_dst3)
  samp_data1$response <- samp_data1[,dat_cl_namh]
  samp_data <- samp_data1
  print("comp high 1")
  cress_hi1 <- cress_func(sim_data, samp_data, paste0(data_in, "_hi"), samp_dst1$dataDist, salsa2dlist1, samp_dst1$knotDist, all_dst1)
  samp_data2$response <- samp_data2[,dat_cl_namh]
  samp_data <- samp_data2
  print("comp high 2")
  cress_hi2 <- cress_func(sim_data, samp_data, paste0(data_in, "_hi"), samp_dst2$dataDist, salsa2dlist2, samp_dst2$knotDist, all_dst2)
  samp_data3$response <- samp_data3[,dat_cl_namh]
  samp_data <- samp_data3
  print("comp high 3")
  cress_hi3 <- cress_func(sim_data, samp_data, paste0(data_in, "_hi"), samp_dst3$dataDist, salsa2dlist3, samp_dst3$knotDist, all_dst3)
  
  # write out results 
  samp_out1 <- cbind(rep(paste0(ii, "_1"), 10000), cress_low1$output_mat, cress_hi1$output_mat)
  colnames(samp_out1)[1] <- c("sim_no")
  samp_out2 <- cbind(rep(paste0(ii, "_2"), 10000), cress_low2$output_mat, cress_hi2$output_mat)
  colnames(samp_out2)[1] <- c("sim_no")
  samp_out3 <- cbind(rep(paste0(ii, "_3"), 10000), cress_low3$output_mat, cress_hi3$output_mat)
  colnames(samp_out3)[1] <- c("sim_no")
  samp_out <- rbind(samp_out1, samp_out2, samp_out3)
  sim_res1 <- cbind(paste0(ii, "_1"), cress_low1$output_dat, cress_hi1$output_dat)
  colnames(sim_res1)[1] <- "sim_no"
  sim_res2 <- cbind(paste0(ii, "_2"), cress_low2$output_dat, cress_hi2$output_dat)
  colnames(sim_res2)[1] <- "sim_no"
  sim_res3 <- cbind(paste0(ii, "_3"), cress_low3$output_dat, cress_hi3$output_dat)
  colnames(sim_res3)[1] <- "sim_no"
  sim_res <- rbind(sim_res1, sim_res2, sim_res3)
  ci_low1 <- cbind.data.frame(rep(paste0("cresslow1_", ii), nrow(cress_low1$ci_out)), cress_low1$ci_out)
  colnames(ci_low1)[1] <- "sim_no"
  ci_low2 <- cbind.data.frame(rep(paste0("cresslow2_", ii), nrow(cress_low2$ci_out)), cress_low2$ci_out)
  colnames(ci_low2)[1] <- "sim_no"
  ci_low3 <- cbind.data.frame(rep(paste0("cresslow3_", ii), nrow(cress_low3$ci_out)), cress_low3$ci_out)
  colnames(ci_low3)[1] <- "sim_no"
  ci_high1 <- cbind.data.frame(rep(paste0("cresshigh1_", ii), nrow(cress_hi1$ci_out)), cress_hi1$ci_out)
  colnames(ci_high1)[1] <- "sim_no"
  ci_high2 <- cbind.data.frame(rep(paste0("cresshigh2_", ii), nrow(cress_hi2$ci_out)), cress_hi2$ci_out)
  colnames(ci_high2)[1] <- "sim_no"
  ci_high3 <- cbind.data.frame(rep(paste0("cresshigh3_", ii), nrow(cress_hi3$ci_out)), cress_hi3$ci_out)
  colnames(ci_high3)[1] <- "sim_no"
  ci_comb <- rbind.data.frame(ci_low1, ci_low2, ci_low3, ci_high1, ci_high2, ci_high3)
  
  models_list <- list(cress_low1$model, cress_hi1$model,
                      cress_low2$model, cress_hi2$model,
                      cress_low3$model, cress_hi3$model)
  names(models_list) <- c(paste0("cresslow1_", ii), paste0("cresshigh1_", ii),
                          paste0("cresslow2_", ii), paste0("cresshigh2_", ii),
                          paste0("cresslow3_", ii), paste0("cresshigh3_", ii))
  
  conf_mat <- rbind(cress_low1$conf_mat, cress_low2$conf_mat, cress_low3$conf_mat,
                    cress_hi1$conf_mat, cress_hi2$conf_mat, cress_hi3$conf_mat)
  
  if (ii == 1) {
    output_matrix <- samp_out
    output_store <- sim_res
    models_out <- models_list
    knotgrid1list <- list(knotgrid1)
    names(knotgrid1list) <- "kgrid1_1"
    knotgrid2list <- list(knotgrid2)
    names(knotgrid2list) <- "kgrid2_1"
    knotgrid3list <- list(knotgrid3)
    names(knotgrid3list) <- "kgrid3_1"
    conf_mat_cress <- conf_mat
    ci_out <- ci_comb
  } else {
    output_matrix <- rbind(output_matrix, samp_out)
    output_store <- rbind(output_store, sim_res)
    models_out <- append(models_out, models_list)
    knotgrid1list <- append(knotgrid1list, knotgrid1)
    names(knotgrid1list)[ii] <- paste0("kgrid1_", ii)
    knotgrid2list <- append(knotgrid2list, knotgrid2)
    names(knotgrid2list)[ii] <- paste0("kgrid2_", ii)
    knotgrid3list <- append(knotgrid3list, knotgrid3)
    names(knotgrid3list)[ii] <- paste0("kgrid3_", ii)
    conf_mat_cress <- conf_mat_cress + conf_mat
    ci_out <- rbind.data.frame(ci_out, ci_comb)
  }

}

write.csv(output_matrix, outmat_name, row.names=FALSE)
write.csv(output_store, outstr_name, row.names=FALSE)
write.csv(conf_mat_cress, confmat_name, row.names=FALSE)
valz <- c(repz, fitmeas, sksk, mnk, mxk, gpgp, knkn)
write.csv(valz, valz_name, row.names=FALSE)
write.csv(store_state, paste0(file_in_dir, "store_state_crco.csv"), row.names=FALSE)
write.csv(ci_out, ci_name, row.names=FALSE)

kgrid1_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid1.csv")
kgrid2_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid2.csv")
kgrid3_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid3.csv")

write.csv(knotgrid1list, kgrid1_name, row.names=FALSE)
write.csv(knotgrid2list, kgrid2_name, row.names=FALSE)
write.csv(knotgrid3list, kgrid3_name, row.names=FALSE)

out_deets <- list()

for (md in 1:length(models_out)){
  nam <- names(models_out)[md]
  mod <- models_out[[nam]]
  knotpoz <- mod$splineParams[[1]]$knotPos
  knotz <- mod$splineParams[[1]]$knotgrid[knotpoz,]
  radz <- mod$splineParams[[1]]$radii[mod$splineParams[[1]]$radiusIndices]
  nknotz <- length(knotpoz)
  coefz <- mod$bestModel@coefficients
  list_out <- list("knots"=knotz, "no_knots"=nknotz, "radii"=radz)
  out_deets[[nam]] <- list_out
  summary_out <- summary(mod$bestModel)
  coefstd <- summary_out@coef3[,2]
  covz <- c(summary_out@cov.unscaled)
  if (md==1){
    namez <- nam
    nknt_out <- nknotz
    knotz_out <- knotz
    radz_out <- radz
    coefz_out <- coefz
    coefstd_out <- coefstd
    covz_out <- covz
  } else {
    namez <- c(namez, nam)
    nknt_out <- c(nknt_out, nknotz)
    knotz_out <- rbind(knotz_out, knotz)
    radz_out <- c(radz_out, radz)
    coefz_out <- c(coefz_out, coefz)
    coefstd_out <- c(coefstd_out, coefstd)
    covz_out <- c(covz_out, covz)
  }
}

listout <- cbind(namez, nknt_out)
matout <- cbind(knotz_out, radz_out)
colnames(matout) <- c("Xknot", "Yknot", "radknot")
write.csv(listout, listout_name, row.names=FALSE)
write.csv(matout, matout_name, row.names=FALSE)
write.csv(coefz_out, coefout_name, row.names=FALSE)
write.csv(coefstd_out, coefstd_name, row.names=FALSE)
write.csv(covz_out, cov_name, row.names=FALSE)


#######################################################################


data_in = "island"
input_name <- paste0(file_in_dir, data_in, "_input_data.csv")
inputdata_class_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_data.csv")
inputdata_numeric_lw <- paste0(file_in_dir, "low_s2n_", data_in, "_datanum.csv")
inputdata_class_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_data.csv")
inputdata_numeric_hi <- paste0(file_in_dir, "high_s2n_", data_in, "_datanum.csv")
outmat_name <- paste0(file_save_dir, "cress_", data_in, "_mat2_new_cv.csv")
outstr_name <- paste0(file_save_dir, "cress_", data_in, "_dat2_new_cv.csv")
confmat_name <- paste0(file_save_dir, "cress_", data_in, "_conf2_new_cv.csv")
valz_name <- paste0(file_save_dir, "cress_", data_in, "_valz2_new_cv.csv")
listout_name <- paste0(file_save_dir, "cress_", data_in, "_deets_list2_new_cv.csv")
matout_name <- paste0(file_save_dir, "cress_", data_in, "_deets_mat2_new_cv.csv")
dist_name <- paste0(file_in_dir, data_in, "_geodistmat.csv")
coefout_name <- paste0(file_save_dir, "cress_", data_in, "_coefs_new_cv.csv")
coefstd_name <- paste0(file_save_dir, "cress_", data_in, "_coefstd_new_cv.csv")
cov_name <- paste0(file_save_dir, "cress_", data_in, "_cov_new_cv.csv")
ci_name <- paste0(file_save_dir, "cress_", data_in, "_cipred_new_cv.csv")

if (file.exists(paste0(file_in_dir, "store_state_crco.csv"))) { 
  store_state <- read.csv(paste0(file_in_dir, "store_state_crco.csv"))
}

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
fitmeas <- "cv.gamMRSea"
sksk <- 6
mnk <- 2
mxk <- 20
gpgp <- 0
knkn <- 300

data_out$x.pos <- data_out$xx
data_out$y.pos <- data_out$yy

store_state <- cbind(store_state, .Random.seed)
colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_1_1")

for (ii in 1:repz) {
#for (ii in 1:3) {
  print(paste0("rep", ii))
  sim_data <- cbind(data_out[,c("xx", "yy", "x.pos", "y.pos", "truthcl_island_lw", "truthcl_island_hi")], 
                    data_lowF[,ii], data_highF[,ii], data_lowN[,ii], data_highN[,ii])
  colnames(sim_data) <- c("xx", "yy", "x.pos", "y.pos", "truthcl_island_lw", "truthcl_island_hi", 
                          "datacl_island_lw", "datacl_island_hi", "dataclnum_island_lw", 
                          "dataclnum_island_hi")
  
  samp_rows1 <- unlist(dataset_rows1i[ii,])
  samp_rows2 <- unlist(dataset_rows2i[ii,])
  samp_rows3 <- unlist(dataset_rows3i[ii,])
  samp_data1 <- sim_data[samp_rows1,]
  samp_data2 <- sim_data[samp_rows2,]
  samp_data3 <- sim_data[samp_rows3,]
  nsamp1 <- length(samp_rows1)
  nsamp2 <- length(samp_rows2)
  nsamp3 <- length(samp_rows3)
  #GeoDistMat_samp1 <- GeoDistMat[samp_rows1,samp_rows1]
  #GeoDistMat_samp2 <- GeoDistMat[samp_rows2,samp_rows2]
  #GeoDistMat_samp3 <- GeoDistMat[samp_rows3,samp_rows3]
  
  knotgrid1 <- cbind(samp_data1$xx, samp_data1$yy)
  # first need to sample rows / cols for this dataset - then knots are all those
  #d2k_samp1 <- as.matrix(GeoDistMat_samp1[,samp_rows1])
  d2k_samp1 <- as.matrix(GeoDistMat[samp_rows1,samp_rows1])
  #k2k_samp1 <- as.matrix(GeoDistMat_samp1[,samp_rows1])
  k2k_samp1 <- GeoDistMat[samp_rows1,samp_rows1]
  d2k_mat1 <- as.matrix(GeoDistMat[, samp_rows1])
  k2k_mat1 <- as.matrix(GeoDistMat[samp_rows1, samp_rows1])
  #store_state <- cbind(store_state, .Random.seed)
  #colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_2_", ii)
  knotgrid2 <- getKnotgrid(coordData = cbind(samp_data2$xx, samp_data2$yy), numKnots = knkn)
  knot_ind2 <- attr(knotgrid2, "points.selected")
  # first need to sample rows / cols for this dataset - then need to sample knots from those
  d2k_samp2 <- as.matrix(GeoDistMat[samp_rows2, samp_rows2])
  d2k_samp2 <- as.matrix(d2k_samp2[, knot_ind2])
  k2k_samp2 <- as.matrix(GeoDistMat[samp_rows2, samp_rows2])
  k2k_samp2 <- as.matrix(k2k_samp2[knot_ind2, knot_ind2])
  d2k_mat2 <- as.matrix(GeoDistMat[, samp_rows2])
  d2k_mat2 <- as.matrix(d2k_mat2[, knot_ind2])
  k2k_mat2 <- as.matrix(GeoDistMat[samp_rows2, samp_rows2])
  k2k_mat2 <- as.matrix(GeoDistMat[knot_ind2, knot_ind2])
  #store_state <- cbind(store_state, .Random.seed)
  #colnames(store_state)[ncol(store_state)] <- paste0("kgrid_", data_in, "_3_", ii)
  knotgrid3 <- getKnotgrid(coordData = cbind(samp_data3$xx, samp_data3$yy), numKnots = knkn)
  knot_ind3 <- attr(knotgrid3, "points.selected")
  # first need to sample rows / cols for this dataset - then need to sample knots from those
  d2k_samp3 <- as.matrix(GeoDistMat[samp_rows3, samp_rows3])
  d2k_samp3 <- as.matrix(d2k_samp3[, knot_ind3])
  k2k_samp3 <- as.matrix(GeoDistMat[samp_rows3, samp_rows3])
  k2k_samp3 <- as.matrix(k2k_samp3[knot_ind3, knot_ind3])
  d2k_mat3 <- as.matrix(GeoDistMat[, samp_rows3])
  d2k_mat3 <- as.matrix(d2k_mat3[, knot_ind3])
  k2k_mat3 <- as.matrix(GeoDistMat[samp_rows3, samp_rows3])
  k2k_mat3 <- as.matrix(k2k_mat3[knot_ind3, knot_ind3])
  # create salsa input for all runs
  salsa2dlist1 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid1, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10,
    max.iter = 100,
    cv.opts=list(K=3)
  )
  salsa2dlist2 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid2, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10,
    max.iter = 100,
    cv.opts=list(K=5)
  )
  salsa2dlist3 <- list(
    fitnessMeasure = fitmeas, 
    knotgrid = knotgrid3, 
    startKnots=sksk, 
    minKnots=mnk,
    maxKnots=mxk, 
    gap=gpgp,
    noradii=10,
    max.iter = 100,
    cv.opts=list(K=5)
  )
  
  # data
  dat_cl_naml <- paste0("datacl_", data_in, "_lw")
  dat_cl_namh <- paste0("datacl_", data_in, "_hi")
  samp_data1$response <- samp_data1[,dat_cl_naml]
  samp_data <- samp_data1
  print("island low 1")
  set.seed(3004)
  cress_low1 <- suppressWarnings(cress_func(sim_data, samp_data, paste0(data_in, "_lw"), d2k_samp1, salsa2dlist1, k2k_samp1, d2k_mat1))
  samp_data2$response <- samp_data2[,dat_cl_naml]
  samp_data <- samp_data2
  print("island low 2")
  set.seed(3004)
  cress_low2 <- suppressWarnings(cress_func(sim_data, samp_data, paste0(data_in, "_lw"), d2k_samp2, salsa2dlist2, k2k_samp2, d2k_mat2))
  samp_data3$response <- samp_data3[,dat_cl_naml]
  samp_data <- samp_data3
  print("island low 3")
  set.seed(3004)
  cress_low3 <- suppressWarnings(cress_func(sim_data, samp_data3, paste0(data_in, "_lw"), d2k_samp3, salsa2dlist3, k2k_samp3, d2k_mat3))
  samp_data1$response <- samp_data1[,dat_cl_namh]
  samp_data <- samp_data1
  print("island high 1")
  set.seed(3004)
  cress_hi1 <- suppressWarnings(cress_func(sim_data, samp_data1, paste0(data_in, "_hi"), d2k_samp1, salsa2dlist1, k2k_samp1, d2k_mat1))
  samp_data2$response <- samp_data2[,dat_cl_namh]
  samp_data <- samp_data2
  print("island high 2")
  set.seed(3004)
  cress_hi2 <- suppressWarnings(cress_func(sim_data, samp_data2, paste0(data_in, "_hi"), d2k_samp2, salsa2dlist2, k2k_samp2, d2k_mat2))
  samp_data3$response <- samp_data3[,dat_cl_namh]
  samp_data <- samp_data3
  print("island high 3")
  set.seed(3004)
  cress_hi3 <- suppressWarnings(cress_func(sim_data, samp_data3, paste0(data_in, "_hi"), d2k_samp3, salsa2dlist3, k2k_samp3, d2k_mat3))
  
  # write out results 
  samp_out1 <- cbind(rep(paste0(ii, "_1"), nrow(data_out)), cress_low1$output_mat, cress_hi1$output_mat)
  colnames(samp_out1)[1] <- c("sim_no")
  samp_out2 <- cbind(rep(paste0(ii, "_2"), nrow(data_out)), cress_low2$output_mat, cress_hi2$output_mat)
  colnames(samp_out2)[1] <- c("sim_no")
  samp_out3 <- cbind(rep(paste0(ii, "_3"), nrow(data_out)), cress_low3$output_mat, cress_hi3$output_mat)
  colnames(samp_out3)[1] <- c("sim_no")
  samp_out <- rbind(samp_out1, samp_out2, samp_out3)
  sim_res1 <- cbind(paste0(ii, "_1"), cress_low1$output_dat, cress_hi1$output_dat)
  colnames(sim_res1)[1] <- "sim_no"
  sim_res2 <- cbind(paste0(ii, "_2"), cress_low2$output_dat, cress_hi2$output_dat)
  colnames(sim_res2)[1] <- "sim_no"
  sim_res3 <- cbind(paste0(ii, "_3"), cress_low3$output_dat, cress_hi3$output_dat)
  colnames(sim_res3)[1] <- "sim_no"
  sim_res <- rbind(sim_res1, sim_res2, sim_res3)
  ci_low1 <- cbind.data.frame(rep(paste0("cresslow1_", ii), nrow(cress_low1$ci_out)), cress_low1$ci_out)
  colnames(ci_low1)[1] <- "sim_no"
  ci_low2 <- cbind.data.frame(rep(paste0("cresslow2_", ii), nrow(cress_low2$ci_out)), cress_low2$ci_out)
  colnames(ci_low2)[1] <- "sim_no"
  ci_low3 <- cbind.data.frame(rep(paste0("cresslow3_", ii), nrow(cress_low3$ci_out)), cress_low3$ci_out)
  colnames(ci_low3)[1] <- "sim_no"
  ci_high1 <- cbind.data.frame(rep(paste0("cresshigh1_", ii), nrow(cress_hi1$ci_out)), cress_hi1$ci_out)
  colnames(ci_high1)[1] <- "sim_no"
  ci_high2 <- cbind.data.frame(rep(paste0("cresshigh2_", ii), nrow(cress_hi2$ci_out)), cress_hi2$ci_out)
  colnames(ci_high2)[1] <- "sim_no"
  ci_high3 <- cbind.data.frame(rep(paste0("cresshigh3_", ii), nrow(cress_hi3$ci_out)), cress_hi3$ci_out)
  colnames(ci_high3)[1] <- "sim_no"
  ci_comb <- rbind.data.frame(ci_low1, ci_low2, ci_low3, ci_high1, ci_high2, ci_high3)
  
  models_list <- list(cress_low1$model, cress_hi1$model,
                      cress_low2$model, cress_hi2$model,
                      cress_low3$model, cress_hi3$model)
  names(models_list) <- c(paste0("cresslow1_", ii), paste0("cresshigh1_", ii),
                          paste0("cresslow2_", ii), paste0("cresshigh2_", ii),
                          paste0("cresslow3_", ii), paste0("cresshigh3_", ii))
  
  conf_mat <- rbind(cress_low1$conf_mat, cress_low2$conf_mat, cress_low3$conf_mat,
                    cress_hi1$conf_mat, cress_hi2$conf_mat, cress_hi3$conf_mat)
  
  if (ii == 1) {
    output_matrix <- samp_out
    output_store <- sim_res
    models_out <- models_list
    knotgrid1list <- as.data.frame(knotgrid1)
    colnames(knotgrid1list) <- c("x.pos", "y.pos")
    knotgrid1list$sim <- rep("kgrid1_1", nrow(knotgrid1))
    knotgrid2list <- as.data.frame(knotgrid2)
    colnames(knotgrid2list) <- c("x.pos", "y.pos")
    knotgrid2list$sim <- rep("kgrid2_1", nrow(knotgrid2))
    knotgrid3list <- as.data.frame(knotgrid3)
    colnames(knotgrid3list) <- c("x.pos", "y.pos")
    knotgrid3list$sim <- rep("kgrid3_1", nrow(knotgrid3))
    conf_mat_cress <- conf_mat
    ci_out <- ci_comb
  } else {
    output_matrix <- rbind(output_matrix, samp_out)
    output_store <- rbind(output_store, sim_res)
    models_out <- append(models_out, models_list)
    knotgrid1sim <- as.data.frame(knotgrid1)
    colnames(knotgrid1sim) <- c("x.pos", "y.pos")
    knotgrid1sim$sim <- rep(paste0("kgrid1_", ii), nrow(knotgrid1))
    knotgrid1list <- rbind.data.frame(knotgrid1list, knotgrid1sim)
    knotgrid2sim <- as.data.frame(knotgrid2)
    colnames(knotgrid2sim) <- c("x.pos", "y.pos")
    knotgrid2sim$sim <- rep(paste0("kgrid2_", ii), nrow(knotgrid2))
    knotgrid2list <- rbind.data.frame(knotgrid2list, knotgrid2sim)
    knotgrid3sim <- as.data.frame(knotgrid3)
    colnames(knotgrid3sim) <- c("x.pos", "y.pos")
    knotgrid3sim$sim <- rep(paste0("kgrid3_", ii), nrow(knotgrid3))
    knotgrid3list <- rbind.data.frame(knotgrid3list, knotgrid3sim)
    conf_mat_cress <- conf_mat_cress + conf_mat
    ci_out <- rbind.data.frame(ci_out, ci_comb)
  }
}

write.csv(output_matrix, outmat_name, row.names=FALSE)
write.csv(output_store, outstr_name, row.names=FALSE)
write.csv(conf_mat_cress, confmat_name, row.names=FALSE)
valz <- c(repz, fitmeas, sksk, mnk, mxk, gpgp, knkn)
write.csv(valz, valz_name, row.names=FALSE)
write.csv(store_state, paste0(file_in_dir, "store_state_cris.csv"), row.names=FALSE)
write.csv(ci_out, ci_name, row.names=FALSE)

kgrid1_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid1.csv")
kgrid2_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid2.csv")
kgrid3_name <- paste0(file_save_dir, "cress_", data_in, "_kgrid3.csv")

write.csv(knotgrid1list, kgrid1_name, row.names=FALSE)
write.csv(knotgrid2list, kgrid2_name, row.names=FALSE)
write.csv(knotgrid3list, kgrid3_name, row.names=FALSE)

out_deets <- list()

grid_island <- data_out[,c("xx", "yy")]

for (md in 1:length(models_out)){
  nam <- names(models_out)[md]
  mod <- models_out[[nam]]
  knotpoz <- mod$splineParams[[1]]$knotPos
  knotz <- mod$splineParams[[1]]$knotgrid[knotpoz,]
  #knot_pt <- attr(mod$splineParams[[1]]$knotgrid, "points.selected")[knotpoz]
  #knotz <- grid_i[knot_pt,]
  radz <- mod$splineParams[[1]]$radii[mod$splineParams[[1]]$radiusIndices]
  nknotz <- length(knotpoz)
  coefz <- mod$bestModel@coefficients
  summary_out <- summary(mod$bestModel)
  coefstd <- summary_out@coef3[,2]
  covz <- c(summary_out@cov.unscaled)
  #knotz_radz <- cbind(knotz, radz)
  #for (kn in 1:length(knotpoz)){
  #  basis_col <- create_basis(grid_i, GeoDistMat, knotz_radz[kn,])
  #  if (kn == 1){
  #    basis_cols <- basis_col
  #  } else {
  #    basis_cols <- cbind(basis_cols, basis_col)
  #  }
  #}
  #basis_cols <- cbind(rep(1, nrow(basis_cols)), basis_cols)
  #rcoefs <- rmvnorm(1000, coefz, summary_out@cov.unscaled)
  #quant.func<- function(x){quantile(x, probs=c(0.025, 0.975))}
  #cis <- apply(rcoefs, 2, quant.func)
  #lw_ci <- matrix(cis[1,], ncol=2, byrow=T)
  #hi_ci <- matrix(cis[2,], ncol=2, byrow=T)
  #low_lim <- basis_cols %*% lw_ci
  #high_lim <- basis_cols %*% hi_ci
  #low_resp_cr <- inverse_link_func(low_lim)
  #high_resp_cr <- inverse_link_func(high_lim)
  #ci_comb <- cbind.data.frame(rep(nam, nrow(low_resp_cr)), low_resp_cr, high_resp_cr)
  #colnames(ci_comb) <- c("model", "lowA", "lowB", "lowC", "highA", "highB", "highC")
  list_out <- list("knots"=knotz, "no_knots"=nknotz, "radii"=radz)
  out_deets[[nam]] <- list_out
  if (md==1){
    namez <- nam
    nknt_out <- nknotz
    knotz_out <- knotz
    radz_out <- radz
    coefz_out <- coefz
    coefstd_out <- coefstd
    covz_out <- covz
  } else {
    namez <- c(namez, nam)
    nknt_out <- c(nknt_out, nknotz)
    knotz_out <- rbind(knotz_out, knotz)
    radz_out <- c(radz_out, radz)
    coefz_out <- c(coefz_out, coefz)
    coefstd_out <- c(coefstd_out, coefstd)
    covz_out <- c(covz_out, covz)
  }
}

listout <- cbind(namez, nknt_out)
matout <- cbind(knotz_out, radz_out)
colnames(matout) <- c("Xknot", "Yknot", "radknot")
write.csv(listout, listout_name, row.names=FALSE)
write.csv(matout, matout_name, row.names=FALSE)
write.csv(coefz_out, coefout_name, row.names=FALSE)
write.csv(coefstd_out, coefstd_name, row.names=FALSE)
write.csv(covz_out, cov_name, row.names=FALSE)


