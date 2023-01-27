library(vegan)
library(tidyverse)
library(ggplot2)
library(mobr)
library(spatialbiodiv)
library(scam)
library(matrixStats)

#' GET_sSBR_DIST()
#' the parameters of the function:
#' @param data the data frame we want to look at, e.g. stav's output
#' @param frag_level we currently look at frag = 0.1 and 0.9, has to be in the data frame
#' @param timestep 41 is currently directly after frag and 100 is at the end of sim

  get_sSBR_dist = function(dataset, frag_level_low, frag_level_high, sample_size) {
 
    # taking only n samples from the dataset, filtering for frag level
    
    simdat_test_low <- dplyr::slice_sample(
      filter(dataset, fragmentation == frag_level_low), n = sample_size)
  
    simdat_test_high <- dplyr::slice_sample(
      filter(dataset, fragmentation == frag_level_high), n = sample_size)
  
    # saving coordinates in new variable for sSBR function
    
    simdat_xy_low <- simdat_test_low %>% dplyr::select(loc_x, loc_y)
    simdat_xy_high <- simdat_test_high %>% dplyr::select(loc_x, loc_y)
  
    sSBR_test_low <- sSBR(comm = simdat_test_low, xy_coords = simdat_xy_low)
    sSBR_test_low$sSBR_data$fragmentation <- frag_level_low
 
    sSBR_test_high <- sSBR(comm = simdat_test_high, xy_coords = simdat_xy_high)
    sSBR_test_high$sSBR_data$fragmentation <- frag_level_high
    
    # creating global ! new_dist data for the curve
  
    new_dist <<- data.frame(distance = seq(0, 99, by = 1))
  
    scam_low_frag <- scam::scam(S ~ s(distance, bs = "mpi"),
                                data = sSBR_test_low$sSBR_data, family = "poisson")
    scam_high_frag <- scam::scam(S ~ s(distance, bs = "mpi"),
                                 data = sSBR_test_high$sSBR_data, family = "poisson")
  
    pred_low_frag <- predict(scam_low_frag, new_dist, type = "response")
    pred_high_frag <- predict(scam_high_frag, new_dist, type = "response")
  
    new_dist$diff_sSBR <<- pred_high_frag - pred_low_frag
    
    # returning the global variable (optional)
  
    return(new_dist) 
    
    }
  
  #' GET_NULL_MODEL()
  #' the parameters of the function (same as above):
  #' @param data the data frame we want to look at, e.g. stav's output
  #' @param frag_level we currently look at frag = 0.1 and 0.9, has to be in the data frame
  #' @param timestep 41 is currently directly after frag and 100 is at the end of sim
  
  get_null_model = function(dataset, frag_level_low, frag_level_high, 
                        sample_size, permutations) {
    
    nd <- matrix(0,
                            nrow = 100,
                            ncol = permutations)
    
    #' this runs the sSBR for n permutations
    
    for(i in 1:permutations) {
      
      data_null <- dataset
      
      frag <- data_null$fragmentation
      
      data_null$fragmentation <- sample(frag, size = nrow(data_null), replace = FALSE)
      
      simdat_test_low <- dplyr::slice_sample(
        filter(data_null, fragmentation == frag_level_low), n = sample_size)
      
      simdat_test_high <- dplyr::slice_sample(
        filter(data_null, fragmentation == frag_level_high), n = sample_size)
      
      simdat_xy_low <- simdat_test_low %>% dplyr::select(loc_x, loc_y)
      simdat_xy_high <- simdat_test_high %>% dplyr::select(loc_x, loc_y)
      
      sSBR_test_low <- sSBR(comm = simdat_test_low, xy_coords = simdat_xy_low)
      sSBR_test_low$sSBR_data$fragmentation <- frag_level_low
      
      sSBR_test_high <- sSBR(comm = simdat_test_high, xy_coords = simdat_xy_high)
      sSBR_test_high$sSBR_data$fragmentation <- frag_level_high
      
      # new_dist doesn't have to be global because ribbonbase will be the output
      
      new_dist <- data.frame(distance = seq(0, 99, by = 1))
      
      scam_low_frag <- scam::scam(S ~ s(distance, bs = "mpi"),
                                  data = sSBR_test_low$sSBR_data, family = "poisson")
      scam_high_frag <- scam::scam(S ~ s(distance, bs = "mpi"),
                                   data = sSBR_test_high$sSBR_data, family = "poisson")
      
      pred_low_frag <- predict(scam_low_frag, new_dist, type = "response")
      pred_high_frag <- predict(scam_high_frag, new_dist, type = "response")
      
      new_dist$diff_sSBR <- pred_high_frag - pred_low_frag
      
      #' i'm creating a table that grows with each permutation so i can later use quantile()
      #' on all permutations for distance x
      #' TIP: in a table [x,y] stands for [row,column]      

      nd[,i] <- new_dist[,2]
      
    }
    
    ribbonbase <<- data.frame(matrix(0,
                                    nrow = nrow(nd),
                                    ncol = 5))
    
    colnames(ribbonbase) <<- c("distance","ymin", "ymax", "minline", "maxline")
    
    #' this saves the distances as one column because geom_ribbon needs an x value
    
    ribbonbase[,1] <<- new_dist[,1]
  
    #' former quantiles() loop now replaced by matrixStats::rowQuantiles()
  
      ribbonbase[,2] <<- rowQuantiles(nd, probs = 0.025)
      ribbonbase[,3] <<- rowQuantiles(nd, probs = 0.975)
      ribbonbase[,4] <<- rowQuantiles(nd, probs = 0) # should do the same as min
      ribbonbase[,5] <<- rowQuantiles(nd, probs = 1) # should do the same as max
    
    # returning ribbonbase (optional)
    
    return(ribbonbase)
  }
  
# run the function, template :   
# get_null_model(dataset = simdat0, frag_level_low = 0.1, frag_level_high = 0.9, timestep = 41, permutations = 50, sample_size = 100)
