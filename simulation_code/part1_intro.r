# This file creates introductory data used by all following code
# In particular it creates grids of points at which to generate simulations
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
library(sf)
library(Rfast)

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "*10^", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output <- sub("1\\*", "", output) # Removes 1*10 but leaves anything else raised to 10 
  parse(text=output)
}


cbbPalette <- c("#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")



if (file.exists(paste0(file_save_dir, "location_grid.csv"))) { 
  grid <- read.csv(paste0(file_save_dir, "location_grid.csv"))
} else { 
  n <- seq(0.005, 0.995, length.out=100)
  xx <- rep(n, times=100)
  yy <- rep(n, each=100)
  grid <- cbind(xx,yy)
  grid <- as.data.frame(grid)
  write.csv(grid, paste0(file_save_dir, "location_grid.csv"), row.names=F)
}

png(filename=paste0(img_save_dir, "Grid.png"))
par(pty="s")
plot(grid$xx,grid$yy, type='p', pch=4, cex=.2, main='Base grid', asp=1, xlab="X coordinate", ylab="Y coordinate")
dev.off()

linestr_create <- function(pt1, pt2){
  linestr <- st_linestring(rbind(st_coordinates(pt1), st_coordinates(pt2)))
  return(linestr)
}

dist_exclude_vec <- function(vec, pt, ez) {
  distpt <- st_distance(vec, pt)
  linechk <- st_sfc(lapply(vec, linestr_create, pt2=pt))
  intchk <- as.vector(st_intersects(ez, linechk, sparse=FALSE))
  distpt <- ifelse(intchk, Inf, distpt)
  return(distpt)
}

if (file.exists(paste0(file_save_dir,"location_grid_island.csv"))) { 
  grid_i <- read.csv(paste0(file_save_dir, "location_grid_island.csv"))
} else { 
  nnn <- nrow(grid)
  pointz <- cbind(grid$xx, grid$yy)
  pt1 <- st_point(pointz[1,])
  pt2 <- st_point(pointz[2,])
  ptlist <- list(pt1, pt2)
  for (pt in 3:nnn){
    pt1 <- st_point(pointz[pt,])
    ptlist[[pt]] <- pt1
  }
  sfpointz <- st_sfc(ptlist)
  #create exclusion zone and remove points in exclusion zone
  #island_in <- st_polygon(list(rbind(c(0.1, 0.3), c(0.15, 0.3), c(0.8, 0.55), 
  #                                   c(0.75, 0.75), c(0.1, 0.4), c(0.1, 0.3))))
  island_in <- st_polygon(list(rbind(c(0.0, 0.3), c(0.15, 0.3), c(0.8, 0.55), 
                                     c(0.75, 0.75), c(0.0, 0.4), c(0.0, 0.3))))
  pointz_keep <- st_disjoint(island_in, sfpointz)
  sfpointz <- sfpointz[unlist(pointz_keep)]
  nnn <- length(sfpointz)
  # create distmat
  distmat1 <- dist_exclude_vec(sfpointz[2:nnn], pt=sfpointz[1], ez=island_in)
  distmat1 <- c(0, distmat1)
  # need to do the first two so that there is a row for the distmat[ii,] line to work 
  distcol <- dist_exclude_vec(sfpointz[3:nnn], sfpointz[2], island_in)
  distcol <- c(distmat1[2], 0, distcol)
  distmat1 <- cbind(distmat1, distcol)
  for (ii in 3:(nnn-1)) {
    print(ii)
    distcol <- dist_exclude_vec(sfpointz[(ii + 1):nnn], sfpointz[ii], island_in)
    distcol <- c(distmat1[ii,], 0, distcol)
    distmat1 <- cbind(distmat1, distcol)
  }
  distmat1 <- cbind(distmat1, c(distmat1[nnn,],0))
  write.csv(distmat1, paste0(file_save_dir, "island_distmat1.csv"), row.names=F)
  positionMat <- st_coordinates(sfpointz)
  xxxx <- positionMat[, 1]
  yyyy <- positionMat[, 2]
  GeoDistMat <- floyd(as.matrix(distmat1))
  write.csv(GeoDistMat, paste0(file_save_dir, "island_geodistmat.csv"), row.names=FALSE)
  grid_i <- as.data.frame(cbind(xxxx, yyyy))
  colnames(grid_i) <- c("xx", "yy")
  write.csv(grid_i, paste0(file_save_dir, "location_grid_island.csv"), row.names=FALSE)
}

png(filename=paste0(img_save_dir, "GridIsland.png"))
par(pty="s")
plot(grid_i$xx,grid_i$yy, type='p', pch=4, cex=.2, main='Points with exclusion zone', asp=1, xlab="X coordinate", ylab="Y coordinate")
dev.off()

nsim <- 100
nnn <- dim(grid)[1]

set.seed(3004)
store_state <- .Random.seed
store_state <- as.data.frame(store_state)
colnames(store_state) <- "Start"

if (file.exists(paste0(file_save_dir, "dataset_rows1.csv"))) {
  dataset_rows1 <- read.csv(paste0(file_save_dir, "dataset_rows1.csv"))
} else {
  nsamp1 <- 100
  dataset_rows1 <- matrix(rep(NA, nsim*nsamp1), nrow=nsim)
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- paste0("predatasetgen")
  # data generation
  for (ii in 1:nsim) {
    dataset_rows1[ii,] <- sample(nnn, nsamp1)
  }
  write.csv(dataset_rows1, paste0(file_save_dir,"dataset_rows1.csv"), row.names=F)
}

if (file.exists(paste0(file_save_dir, "dataset_rows2.csv"))) {
  dataset_rows2 <- read.csv(paste0(file_save_dir, "dataset_rows2.csv"))
} else {
  nsamp2 <- 500
  dataset_rows2 <- matrix(rep(NA, nsim*nsamp2), nrow=nsim)
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- paste0("predatasetgen2")
  # data generation
  for (ii in 1:nsim) {
    dataset_rows2[ii,] <- sample(nnn, nsamp2)
  }
  write.csv(dataset_rows2, paste0(file_save_dir,"dataset_rows2.csv"), row.names=F)
}

if (file.exists(paste0(file_save_dir, "dataset_rows3.csv"))) {
  dataset_rows3 <- read.csv(paste0(file_save_dir, "dataset_rows3.csv"))
} else {
  nsamp3 <- 5000
  dataset_rows3 <- matrix(rep(NA, nsim*nsamp3), nrow=nsim)
  # data generation
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- paste0("predatasetgen3")
  for (ii in 1:nsim) {
    dataset_rows3[ii,] <- sample(nnn, nsamp3)
  }
  write.csv(dataset_rows3, paste0(file_save_dir,"dataset_rows3.csv"), row.names=F)
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- "postdatasetgen"
  write.csv(store_state, paste0(file_save_dir,"store_state.csv"), row.names=F)
}
  
nnn_i <- dim(grid_i)[1]

if (file.exists(paste0(file_save_dir, "dataset_rows1_island.csv"))) {
  dataset_rows1 <- read.csv(paste0(file_save_dir, "dataset_rows1_island.csv"))
} else {
  nsamp1 <- 100
  dataset_rows1 <- matrix(rep(NA, nsim*nsamp1), nrow=nsim)
  # data generation
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- paste0("predatasetgen1island")
  for (ii in 1:nsim) {
    dataset_rows1[ii,] <- sample(nnn_i, nsamp1)
  }
  write.csv(dataset_rows1, paste0(file_save_dir,"dataset_rows1_island.csv"), row.names=F)
}

if (file.exists(paste0(file_save_dir, "dataset_rows2_island.csv"))) {
  dataset_rows2 <- read.csv(paste0(file_save_dir, "dataset_rows2_island.csv"))
} else {
  nsamp2 <- 500
  dataset_rows2 <- matrix(rep(NA, nsim*nsamp2), nrow=nsim)
  # data generation
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- paste0("predatasetgen2island")
  for (ii in 1:nsim) {
    dataset_rows2[ii,] <- sample(nnn_i, nsamp2)
  }
  write.csv(dataset_rows2, paste0(file_save_dir,"dataset_rows2_island.csv"), row.names=F)
}

if (file.exists(paste0(file_save_dir, "dataset_rows3_island.csv"))) {
  dataset_rows3 <- read.csv(paste0(file_save_dir, "dataset_rows3_island.csv"))
} else {
  nsamp3 <- 5000
  dataset_rows3 <- matrix(rep(NA, nsim*nsamp3), nrow=nsim)
  # data generation
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- paste0("predatasetgen3island")
  for (ii in 1:nsim) {
    dataset_rows3[ii,] <- sample(nnn_i, nsamp3)
  }
  write.csv(dataset_rows3, paste0(file_save_dir,"dataset_rows3_island.csv"), row.names=F)
  store_state <- cbind(store_state, .Random.seed)
  colnames(store_state)[ncol(store_state)] <- "postdatasetgenisland"
  write.csv(store_state, paste0(file_save_dir,"store_state.csv"), row.names=F)
}


