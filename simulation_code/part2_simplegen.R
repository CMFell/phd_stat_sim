# This file creates the data for the simple dataset
# In particular it generates 100 simulations of simple data
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


grid <- read.csv(paste0(file_save_dir, "location_grid.csv"))
dataset_rows1 <- read.csv(paste0(file_save_dir, "dataset_rows1.csv"))
dataset_rows2 <- read.csv(paste0(file_save_dir, "dataset_rows2.csv"))
dataset_rows3 <- read.csv(paste0(file_save_dir, "dataset_rows3.csv"))
store_state <- read.csv(paste0(file_save_dir, "store_state.csv"))

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


dist2point <- function(xx, yy, xknot, yknot) {
  distout <- sqrt((xx - xknot)^2 + (yy - yknot)^2)
  return(distout)
}

radgaussbasis <- function(distz, rad) {
  zzz <- exp(-(distz*rad)**2)
  return(zzz)
}

simp2knots <- function(xx, yy, rd, nz, label_str, nsimz, store_state=NULL) {
  
  kt1x <- 0.2
  kt1y <- 0.2
  kt2x <- 0.8
  kt2y <- 0.8
  
  ## Create design matrix
  distK1 <- dist2point(xx, yy, kt1x, kt1y)
  zz_1 <- radgaussbasis(distK1, rd)
  distK2 <- dist2point(xx, yy, kt2x, kt2y)
  zz_2 <- radgaussbasis(distK2, rd)
  
  ## calculate probabilities
  mxpp <- round(max(rowSums(cbind(zz_1+zz_2))), 4) + 0.0001
  if (mxpp > 1) {
    pp1 <- zz_1 / mxpp 
    pp2 <- zz_2 / mxpp
  } else {
    pp1 <- zz_1
    pp2 <- zz_2
  }
  # bodge so can do high or low noise
  if (nz > 1.0) {
    pp1 <- pp1 / nz + 1 / nz
    pp2 <- pp2 / nz + 1 / nz
  }
  pp3 = 1 - pp1 - pp2
  peas <- as.data.frame(cbind(pp1, pp2, pp3))
  colnames(peas) <- str.colnames(paste0("truep_",label_str), 3)
  
  ## calculate true classes
  classez <- c("A", "B", "C")
  peas <- cbind(peas, classez[apply(peas,1,which.max)])
  colnames(peas)[ncol(peas)] <- paste0("truthcl_",label_str)
  
  # Create data classes from truth
  datacl <- matrix(rep(NA, nrow(peas)*nsimz), ncol=nsimz)
  if (!is.null(store_state)) {
    new_store <- .Random.seed
  }
  
  # create dependent vector
  for (sm in 1:nsimz){
    for(j in 1:nrow(peas)){
      data <- rmultinom(1,1,c(peas[j,1:3]))
      datacl[j, sm] <- classez[which.max(data)]
    }
  }
  
  if (!is.null(store_state)) {
    new_store <- cbind(new_store, .Random.seed)
    new_store <- as.data.frame(new_store)
    colnames(new_store) <- c(paste0("dataclgens_start_",label_str), paste0("dataclgens_end_",label_str))
    store_state <<- cbind(store_state, new_store)
  }
  # add first data gen to peas to keep backward compatibility
  peas <- cbind(peas, datacl[,1])
  colnames(peas)[ncol(peas)] <- paste0("datacl_",label_str)
  
  # Create data s numeric for mgcv
  dataclnum <- as.numeric(as.factor(datacl[,1])) - 1
  for (cl in 2:nsimz){
    dataclnum <- cbind(dataclnum, as.numeric(as.factor(datacl[,cl])) - 1)
  }

  peas <- cbind(peas, dataclnum[,1])
  colnames(peas)[ncol(peas)] <- paste0("dataclnum_",label_str)

  return(list(s2n=peas, data_fac_mat=datacl, data_num_mat=dataclnum))  
}

low_s2n <- simp2knots(grid$xx, grid$yy, 3.2, 1.0, "simp_lw", 100, store_state)
high_s2n <- simp2knots(grid$xx, grid$yy, 3.2, 4.0, "simp_hi", 100, store_state)

write.csv(low_s2n$s2n, paste0(file_save_dir, "low_s2n_simp.csv"), row.names=F)
write.csv(high_s2n$s2n, paste0(file_save_dir, "high_s2n_simp.csv"), row.names=F)
write.csv(low_s2n$data_fac_mat, paste0(file_save_dir, "low_s2n_simp_data.csv"), row.names=F)
write.csv(high_s2n$data_fac_mat, paste0(file_save_dir, "high_s2n_simp_data.csv"), row.names=F)
write.csv(low_s2n$data_num_mat, paste0(file_save_dir, "low_s2n_simp_datanum.csv"), row.names=F)
write.csv(high_s2n$data_num_mat, paste0(file_save_dir, "high_s2n_simp_datanum.csv"), row.names=F)
input_data <- cbind(grid, low_s2n$s2n, high_s2n$s2n, low_s2n$data_fac_mat, high_s2n$data_fac_mat, low_s2n$data_num_mat, high_s2n$data_num_mat)
write.csv(input_data, paste0(file_save_dir, "simp_input_data.csv"), row.names=FALSE)
write.csv(store_state, paste0(file_save_dir, "store_state_simp.csv"), row.names=F)

low_s2n <- list()
low_s2n$s2n <- read.csv(paste0(file_save_dir, "low_s2n_simp.csv"))

low_s2n_grid <- cbind.data.frame(low_s2n$s2n, grid)

graph_lo_A_simp <- ggplot(low_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=truep_simp_lw_1), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, name = "", limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_lo_B_simp <- ggplot(low_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=truep_simp_lw_2), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, name = "", limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_lo_C_simp <- ggplot(low_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=truep_simp_lw_3), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, name = "", limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "GraphLoASimp.png"))
graph_lo_A_simp
dev.off()
png(filename=paste0(img_save_dir, "GraphLoBSimp.png"))
graph_lo_B_simp
dev.off()
png(filename=paste0(img_save_dir, "GraphLoCSimp.png"))
graph_lo_C_simp
dev.off()

high_s2n <- list()
high_s2n$s2n <- read.csv(paste0(file_save_dir, "high_s2n_simp.csv"))

high_s2n_grid <- cbind.data.frame(high_s2n$s2n, grid)

graph_hi_A_simp <- ggplot(high_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=truep_simp_hi_1), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, name = "", limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_hi_B_simp <- ggplot(high_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=truep_simp_hi_2), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, name = "", limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

graph_hi_C_simp <- ggplot(high_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  # ggtitle(paste("GFRC Distance to any water")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_summary_2d(aes(z=truep_simp_hi_3), fun=max, binwidth = c(0.02, 0.02)) +
  scale_fill_distiller(palette = "YlOrBr", direction=1, name = "", limits=c(0,1)) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "GraphHiASimp.png"))
graph_hi_A_simp
dev.off()
png(filename=paste0(img_save_dir, "GraphHiBSimp.png"))
graph_hi_B_simp
dev.off()
png(filename=paste0(img_save_dir, "GraphHiCSimp.png"))
graph_hi_C_simp
dev.off()

graph_lo_trucl_simp <- ggplot(low_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  ggtitle(paste("True class for simple case with low noise")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=truthcl_simp_lw), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title="true class")) +
  theme(legend.position="right", text = element_text(size=20)) 

graph_lo_datcl_simp <- ggplot(low_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  #ggtitle(paste("Example of simulated data \n simple case with low noise")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_lw), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20)) 

png(filename=paste0(img_save_dir, "GraphLoTruClSimp.png"))
graph_lo_trucl_simp
dev.off()
png(filename=paste0(img_save_dir, "GraphLoDatClSimp.png"))
graph_lo_datcl_simp
dev.off()

graph_hi_trucl_simp <- ggplot(high_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") +  
  ggtitle(paste("True class for simple case with high noise")) +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=truthcl_simp_hi), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title="true class")) +
  theme(legend.position="right", text = element_text(size=20))

graph_hi_datcl_simp <- ggplot(high_s2n_grid, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  # ggtitle(paste("Example of simulated data \n simple case with high noise")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_hi), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "GraphHiTruClSimp.png"))
graph_hi_trucl_simp
dev.off()
png(filename=paste0(img_save_dir, "GraphHiDatClSimp.png"))
graph_hi_datcl_simp
dev.off()


file_in_dir <- paste0(base_dir_on_comp, "OneDrive - University of St Andrews/PhD/Notes/stats_files/")

if (file.exists(paste0(file_in_dir, "dataset_rows1.csv"))) {
  dataset_rows1 <- read.csv(paste0(file_in_dir, "dataset_rows1.csv"))
  dataset_rows2 <- read.csv(paste0(file_in_dir, "dataset_rows2.csv"))
  dataset_rows3 <- read.csv(paste0(file_in_dir, "dataset_rows3.csv"))
}

# example sample sizes
low_s2n_grid_samp1 = low_s2n_grid[unlist(dataset_rows1[1,]),]
low_s2n_grid_samp2 = low_s2n_grid[unlist(dataset_rows2[1,]),]
low_s2n_grid_samp3 = low_s2n_grid[unlist(dataset_rows3[1,]),]
high_s2n_grid_samp1 = high_s2n_grid[unlist(dataset_rows1[1,]),]
high_s2n_grid_samp2 = high_s2n_grid[unlist(dataset_rows2[1,]),]
high_s2n_grid_samp3 = high_s2n_grid[unlist(dataset_rows3[1,]),]

graph_lo_datcl_simp1 <- ggplot(low_s2n_grid_samp1, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  #ggtitle(paste("Example of sample of simulated data \n simple case with low noise 100 datapoints")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_lw), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20)) 
graph_lo_datcl_simp2 <- ggplot(low_s2n_grid_samp2, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  #ggtitle(paste("Example of sample of simulated data \n simple case with low noise 500 datapoints")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_lw), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20))
graph_lo_datcl_simp3 <- ggplot(low_s2n_grid_samp3, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  #ggtitle(paste("Example of sample of simulated data \n simple case with low noise 5000 datapoints")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_lw), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20))
graph_hi_datcl_simp1 <- ggplot(high_s2n_grid_samp1, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  ggtitle(paste("Example of sample of simulated data \n simple case with high noise 100 datapoints")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_hi), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20))
graph_hi_datcl_simp2 <- ggplot(high_s2n_grid_samp2, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  ggtitle(paste("Example of sample of simulated data \n simple case with high noise 500 datapoints")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_hi), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20))
graph_hi_datcl_simp3 <- ggplot(high_s2n_grid_samp3, aes(xx, yy)) +
  theme_bw() + theme(legend.key=element_blank()) +
  coord_fixed(ratio = 1) +
  xlim(0, 1) + ylim(0, 1) + xlab("X coordinate") + ylab("Y coordinate") + 
  ggtitle(paste("Example of sample of simulated data \n simple case with high noise 5000 datapoints")) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(colour=datacl_simp_hi), shape=15, size=2) +
  scale_color_cbb(guide = guide_legend(title=" simulated \n data class")) +
  theme(legend.position="right", text = element_text(size=20))

png(filename=paste0(img_save_dir, "GraphLoDatClSimp1.png"))
graph_lo_datcl_simp1
dev.off()
png(filename=paste0(img_save_dir, "GraphLoDatClSimp2.png"))
graph_lo_datcl_simp2
dev.off()
png(filename=paste0(img_save_dir, "GraphLoDatClSimp3.png"))
graph_lo_datcl_simp3
dev.off()
png(filename=paste0(img_save_dir, "GraphHiDatClSimp1.png"))
graph_hi_datcl_simp1
dev.off()
png(filename=paste0(img_save_dir, "GraphHiDatClSimp2.png"))
graph_hi_datcl_simp2
dev.off()
png(filename=paste0(img_save_dir, "GraphHiDatClSimp3.png"))
graph_hi_datcl_simp3
dev.off()

