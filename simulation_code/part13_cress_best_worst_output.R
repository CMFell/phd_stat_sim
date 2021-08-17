# This file generates images of results of a single simulation
# Used to show results on simulations that do well and badly

# base directory to save all output
base_dir_on_comp <- "C:/Users/christina/"
base_dir_sim <- base_dir_on_comp + "multinomial_simulations/"
# subdirectories to save csv output and generated images
file_save_dir <- paste0(base_dir_sim, "csv_files/")
img_save_dir <- paste0(base_dir_sim, "image_files/")

library(tidyverse)
library(viridis)
library(hexbin)
library(car)
library(mvtnorm)
library(patchwork)
library(scales)

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output <- sub("1\\*", "", output) # Removes 1*10 but leaves anything else raised to 10 
  parse(text=output)
}


#cbbPalette <- c("#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")
cbbPalette_colours <- c(dark_blue="#004949",mid_green="#65aa70",mid_pink="#ff6db6",light_pink="#ffb6db",
                        dark_purple="#490092","#006ddb",mid_blue="#b66dff",light_blue="#afeeee",
                        dark_red="#920000",brown="#924900",orange="#db6d00",green="#24ff24",yellow="#ffff6d", 
                        black="#000000")
cbb_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (cbbPalette_colours)
  
  cbbPalette_colours[cols]
}
cbb_palettes <- list(three_class = cbb_cols("dark_red", "mid_green", "light_blue"))
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


file_in_dir <- file_save_dir
img_out_dir <- paste0(img_save_dir,"best_worst_simp/")

# Worst and Best results for CRESS

tab_out_disagg <- read.csv(paste0(img_save_dir, "disagg_comparison_table_new.csv"))
tab_out_disagg_rss <- read.csv(paste0(img_save_dir, "rss_disagg_comparison_table_new.csv"))

tab_out_disagg_simp_lw_small <- tab_out_disagg %>% filter(dataset == "simp", noise == "lw", noPoints == 100)
tab_sls_mgcv <- tab_out_disagg_simp_lw_small %>% filter(model == "mgcv")
tab_sls_cress <- tab_out_disagg_simp_lw_small %>% filter(model == "cress")
tab_sls_svm <- tab_out_disagg_simp_lw_small %>% filter(model == "svm")
worst_cr_sls <- tab_sls_mgcv$accuracy_to_truth_all - tab_sls_cress$accuracy_to_truth_all
sim_sls <- which.max(worst_cr_sls) # 17 example of cress doing badly
sim_sls_mn <- which.min(worst_cr_sls) # 24 example of cress doing well
print(paste(tab_sls_mgcv$accuracy_to_truth_all[sim_sls], tab_sls_cress$accuracy_to_truth_all[sim_sls], tab_sls_svm$accuracy_to_truth_all[sim_sls]))
print(paste(tab_sls_mgcv$accuracy_to_truth_all[sim_sls_mn], tab_sls_cress$accuracy_to_truth_all[sim_sls_mn], tab_sls_svm$accuracy_to_truth_all[sim_sls_mn]))


tab_out_disagg_simp_hi_small <- tab_out_disagg %>% filter(dataset == "simp", noise == "hi", noPoints == 100)
tab_shs_mgcv <- tab_out_disagg_simp_hi_small %>% filter(model == "mgcv")
tab_shs_cress <- tab_out_disagg_simp_hi_small %>% filter(model == "cress")
tab_shs_svm <- tab_out_disagg_simp_hi_small %>% filter(model == "svm")
best_cr_shs <- tab_shs_cress$accuracy_to_truth_all - tab_shs_mgcv$accuracy_to_truth_all
sim_shs <- which.max(best_cr_shs) # 93 example of cress doing well
sim_shs_mn <- which.min(best_cr_shs) # 63 example of cress doing badly
print(paste(tab_shs_mgcv$accuracy_to_truth_all[sim_shs], tab_shs_cress$accuracy_to_truth_all[sim_shs], tab_shs_svm$accuracy_to_truth_all[sim_shs]))
print(paste(tab_shs_mgcv$accuracy_to_truth_all[sim_shs_mn], tab_shs_cress$accuracy_to_truth_all[sim_shs_mn], tab_shs_svm$accuracy_to_truth_all[sim_shs_mn]))


tab_out_disagg_simp_hi_large <- tab_out_disagg %>% filter(dataset == "simp", noise == "hi", noPoints == 5000)
tab_shl_mgcv <- tab_out_disagg_simp_hi_small %>% filter(model == "mgcv")
tab_shl_cress <- tab_out_disagg_simp_hi_small %>% filter(model == "cress")
tab_shl_svm <- tab_out_disagg_simp_hi_small %>% filter(model == "svm")
worst_cr_shl <- tab_shl_mgcv$accuracy_to_truth_all - tab_shl_cress$accuracy_to_truth_all
sim_shl <- which.max(worst_cr_shl) # 63 example of cress doing badly
sim_shl_mn <- which.min(worst_cr_shl) # 93 example of cress doing badly
print(paste(tab_shl_mgcv$accuracy_to_truth_all[sim_shl_mn], tab_shl_cress$accuracy_to_truth_all[sim_shl_mn], tab_shl_svm$accuracy_to_truth_all[sim_shl_mn]))


dataset_rows1 <- read.csv(paste0(file_save_dir,"dataset_rows1.csv"))
dataset_rows3 <- read.csv(paste0(file_save_dir,"dataset_rows3.csv"))

deets_simp <- read.csv(paste0(file_in_dir,"cress_simp_deets_list2_new_cv.csv"))

outmat_nam_save <- ""

# cipred_matrix <- read.csv(paste0(file_save_dir, "cress_simp_cipred_new_cv.csv"), sep=",", header=TRUE)

draw_bases <- FALSE

sim_nos <- c(sim_sls, sim_shs, sim_shl, sim_sls_mn, sim_shs_mn, sim_shl_mn)
samp_sz <- c(1, 1, 3, 1, 1, 3)
noize <- c("lw", "hi", "hi", "lw", "hi", "hi")
no_pts <- c(100, 100, 5000, 100, 100, 5000)

ci_dset_save <- "simp"

for (sm in 1:length(sim_nos)){
  sim <- sim_nos[sm]
  
  grid <- read.csv(paste0(file_save_dir,"location_grid.csv"))
  
  file_dir <- file_in_dir
  data_out <- read.csv(paste0(file_save_dir, "simp_input_data.csv"), sep=",", header=TRUE)

  models <- c("cress", "mgcv", "svm")
  for (mod in models){
    if (mod == 'cress') {
      outmat_name <- paste0(file_dir, mod, "_simp_mat2_new_cv.csv")
    } else {
      outmat_name <- paste0(file_dir, mod, "_simp_mat2_new.csv")
    }
    output_matrix <- read.csv(outmat_name, sep=",", header=TRUE)
    
    sim_ex <- output_matrix[output_matrix$sim_no == paste0(sim, "_", samp_sz[sm]),]
    
    predp_col_namA <- paste0("predp_", mod, "_simp_", noize[sm], "_1")
    predp_col_namB <- paste0("predp_", mod, "_simp_", noize[sm], "_2")
    predp_col_namC <- paste0("predp_", mod, "_simp_", noize[sm], "_3")
    trup_col_namA <- paste0("truep_simp_", noize[sm], "_1")
    trup_col_namB <- paste0("truep_simp_", noize[sm], "_2")
    trup_col_namC <- paste0("truep_simp_", noize[sm], "_3")
    tru_col_nam_cl <- paste0("truthcl_simp_", noize[sm])
    biasA <- sim_ex[,predp_col_namA] - data_out[,trup_col_namA]
    biasB <- sim_ex[,predp_col_namB] - data_out[,trup_col_namB]
    biasC <- sim_ex[,predp_col_namC] - data_out[,trup_col_namC]
    pred_cl <- c("A", "B", "C")[max.col(sim_ex[,c(predp_col_namA, predp_col_namB, predp_col_namC)])]
    misclass <- pred_cl == data_out[,tru_col_nam_cl]
    meanerror <- (abs(biasA) + abs(biasB) + abs(biasB)) / 3
    
    df_to_plot <- cbind.data.frame(grid$xx, grid$yy, biasA, biasB, biasC, meanerror, 
                                   sim_ex[,predp_col_namA], sim_ex[,predp_col_namB], sim_ex[,predp_col_namC], 
                                  pred_cl, misclass)
    colnames(df_to_plot) <- c("xx", "yy", "bias_A", "bias_B", "bias_C", "meanerror",
                             "pred_A", "pred_B", "pred_C", "pred_cl", "misclass")

    if (samp_sz[sm] == 1) {
      samplerows <- unlist(dataset_rows1[sim,])
    } else {
      samplerows <- unlist(dataset_rows3[sim,])
    }
    if (noize[sm] == "lw") {
      simcolname <- paste0("X",sim)
      simdata <- data_out[,simcolname]
    } else{
      simcolname <- paste0("X",sim, ".1")
      simdata <- data_out[,simcolname]
    }
    simdf <- cbind.data.frame(grid$xx, grid$yy, simdata)
    simdf <- simdf[samplerows,]
    colnames(simdf) <- c("xx", "yy", "simdata")
    
    # Data points to fit
    sim_pts <- ggplot(simdf, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      geom_point(aes(colour=simdata), shape=15, size=2) +
      scale_color_cbb(guide = guide_legend(title=" data for \n fitting"), palette="three_class") +
      theme(legend.position="right", text = element_text(size=20))
    sim_pts_name <- paste0("FittedPtsSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, sim_pts_name))
    print(sim_pts)
    dev.off()
    
    # Predicted probability class A
    graph_A <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_summary_2d(aes(z=pred_A), fun=max, binwidth = c(0.02, 0.02)) +
      scale_fill_distiller(palette = "YlOrBr", direction=1, name = "probability", limits=c(0,1), oob=squish) +
      theme(legend.position="right", text = element_text(size=20))
    graph_A_name <- paste0("PredProbASimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, graph_A_name))
    print(graph_A)
    dev.off()
    
    # Predicted probability class B
    graph_B <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_summary_2d(aes(z=pred_B), fun=max, binwidth = c(0.02, 0.02)) +
      scale_fill_distiller(palette = "YlOrBr", direction=1, name = "probability", limits=c(0,1), oob=squish) +
      theme(legend.position="right", text = element_text(size=20))
    graph_B_name <- paste0("PredProbBSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, graph_B_name))
    print(graph_B)
    dev.off()
    
    # Predicted probability class C
    graph_C <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_summary_2d(aes(z=pred_C), fun=max, binwidth = c(0.02, 0.02)) +
      scale_fill_distiller(palette = "YlOrBr", direction=1, name = "probability", limits=c(0,1), oob=squish) +
      theme(legend.position="right", text = element_text(size=20))
    graph_C_name <- paste0("PredProbCSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, graph_C_name))
    print(graph_C)
    dev.off()
    
    # Predicted data class
    pred_pts <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      geom_point(aes(colour=pred_cl), shape=15, size=2) +
      scale_color_cbb(guide = guide_legend(title=" predicted \n data class"), palette="three_class") +
      theme(legend.position="right", text = element_text(size=20))
    pred_pts_name <- paste0("PredPtsSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, pred_pts_name))
    print(pred_pts)
    dev.off()
    
    # Error class A
    error_A <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_summary_2d(aes(z=bias_A), fun=max, binwidth = c(0.02, 0.02)) +
      scale_fill_distiller(palette = "YlOrBr", direction=1, name = "probability", limits=c(0,0.25), oob=squish) +
      theme(legend.position="right", text = element_text(size=20))
    error_A_name <- paste0("ErrorASimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, error_A_name))
    print(error_A)
    dev.off()    
    
    # Error class B
    error_B <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_summary_2d(aes(z=bias_B), fun=max, binwidth = c(0.02, 0.02)) +
      scale_fill_distiller(palette = "YlOrBr", direction=1, name = "probability", limits=c(0,0.25), oob=squish) +
      theme(legend.position="right", text = element_text(size=20))
    error_B_name <- paste0("ErrorBSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, error_B_name))
    print(error_B)
    dev.off() 
    
    # Error class C
    error_C <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_summary_2d(aes(z=bias_C), fun=max, binwidth = c(0.02, 0.02)) +
      scale_fill_distiller(palette = "YlOrBr", direction=1, name = "probability", limits=c(0,0.25), oob=squish) +
      theme(legend.position="right", text = element_text(size=20))
    error_C_name <- paste0("ErrorCSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, error_C_name))
    print(error_C)
    dev.off() 
    
    # Misclassifications
    mc_graph <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      geom_point(aes(colour=misclass), shape=15, size=2) +
      scale_color_cbb(guide = guide_legend(title=""), palette="three_class", labels=c("incorrect", "correct")) +
      theme(legend.position="right", text = element_text(size=20))
    mc_graph_name <- paste0("MisClSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, mc_graph_name))
    print(mc_graph)
    dev.off() 
    
    # Mean error
    mean_error <- ggplot(df_to_plot, aes(xx, yy)) +
      theme_bw() + theme(legend.key=element_blank()) +
      coord_fixed(ratio = 1) +
      xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
      theme(plot.title = element_text(hjust = 0.5)) + 
      stat_summary_2d(aes(z=meanerror), fun=max, binwidth = c(0.02, 0.02)) +
      scale_fill_distiller(palette = "YlOrBr", direction=1, name = "probability", limits=c(0,0.5), oob=squish) +
      theme(legend.position="right", text = element_text(size=20))
    mean_error_name <- paste0("MeanErrSimNo", sim, "NoiseLev", noize[sm], "NPts", no_pts[sm], "Model", mod, ".png")
    png(filename=paste0(img_out_dir, mean_error_name))
    print(mean_error)
    dev.off()
  
  }
  
}


