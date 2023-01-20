#' library(vegan)
#' library(tidyverse)
#' library(ggplot2)
#' library(mobr)
#' library(spatialbiodiv)
#' library(scam)

#' setwd("C:/Users/aylin/Desktop/frag_mod_scale-main/simulation_1")
#' simdat0 <- read.csv("sSBR7_rep_1_output_sample.csv", header=TRUE, sep = ",")

#' next step: isolate the sSBR function into an additional one so it's cleaner!
#'

#' sSBR_curve = function(dataset, frag_level_low, frag_level_high, timestep,
                      #' permutations, sample_size)

  #' the parameters of the function:
  #' @param data the data frame we want to look at, e.g. stav's output
  #' @param frag_level we currently look at frag = 0.1 and 0.9, has to be in the data frame
  #' @param timestep 41 is currently directly after frag and 100 is at the end of sim

  null_model = function(dataset, frag_level_low, frag_level_high, timestep,
                        permutations, sample_size) {

    p <- ggplot()

    #' this is what will become the ribbon
    #' when testing parts of the code, use these presets:
    #'    dataset <- simdat0
    #'    frag_level_low <- 0.1
    #'    frag_level_high <- 0.9
    #'    timestep <- 41
    #'    sample_size <- 50

    nd <- data.frame(matrix(0,
                            nrow = 100,
                            ncol = 101)

    colnames(nd) <- c("distance",1:permutations)

    #' this runs the sSBR for n amount of permutations

    for(i in 1:100) {

      data_null <- dataset

      frag <- data_null$fragmentation

      data_null$fragmentation <- sample(frag, size = nrow(data_null), replace = FALSE)

      simdat_test_low <- dplyr::slice_sample(
        filter(data_null, fragmentation == frag_level_low & step == timestep), n = sample_size)

      simdat_test_high <- dplyr::slice_sample(
        filter(data_null, fragmentation == frag_level_high & step == timestep), n = sample_size)

      # simdat_test_abund_low <- colSums(dplyr::select(
      #   simdat_test_low, sp_1:sp_1000))
      #
      # simdat_test_abund_high <- colSums(dplyr::select(
      #   simdat_test_high, sp_1:sp_1000))

      simdat_xy_low <- simdat_test_low %>% dplyr::select(loc_x, loc_y)
      simdat_xy_high <- simdat_test_high %>% dplyr::select(loc_x, loc_y)

      sSBR_test_low <- sSBR(comm = simdat_test_low, xy_coords = simdat_xy_low)
      sSBR_test_low$sSBR_data$fragmentation <- frag_level_low

      sSBR_test_high <- sSBR(comm = simdat_test_high, xy_coords = simdat_xy_high)
      sSBR_test_high$sSBR_data$fragmentation <- frag_level_high

      #' the min_dist code from Felix's Github
      #' we don't actually need this at the moment so i commented it out
      #'
      #'      min_dist <- round(min(c(max(sSBR_test_high$sSBR_data$distance),
      #'                              max(sSBR_test_high$sSBR_data$distance))))
      #'
      #' instead we're using new_dist to have the same distances for all permutations

      new_dist <- data.frame(distance = seq(0, 99, by = 1)) # by = 1  or 0.5 instead of length



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

      nd[,1] <- new_dist[,1]
      nd[,i+1] <- new_dist[,2]

    }

    #' this is the actual curve that we compare to the nullmodel

    simdat_test_low <- dplyr::slice_sample(
      filter(dataset, fragmentation == frag_level_low & step == timestep), n = sample_size)

    simdat_test_high <- dplyr::slice_sample(
      filter(dataset, fragmentation == frag_level_high & step == timestep), n = sample_size)

    simdat_test_abund_low <- colSums(dplyr::select(
      simdat_test_low, sp_1:sp_1000))

    simdat_test_abund_high <- colSums(dplyr::select(
      simdat_test_high, sp_1:sp_1000))

    simdat_xy_low <- simdat_test_low %>% dplyr::select(loc_x, loc_y)
    simdat_xy_high <- simdat_test_high %>% dplyr::select(loc_x, loc_y)

    sSBR_test_low <- sSBR(comm = simdat_test_low, xy_coords = simdat_xy_low)
    sSBR_test_low$sSBR_data$fragmentation <- frag_level_low

    sSBR_test_high <- sSBR(comm = simdat_test_high, xy_coords = simdat_xy_high)
    sSBR_test_high$sSBR_data$fragmentation <- frag_level_high

    min_dist <- round(min(c(max(sSBR_test_high$sSBR_data$distance),
                            max(sSBR_test_high$sSBR_data$distance))))

    new_dist <- data.frame(distance = seq(0, 99, by = 1))

    scam_low_frag <- scam::scam(S ~ s(distance, bs = "mpi"),
                                data = sSBR_test_low$sSBR_data, family = "poisson")
    scam_high_frag <- scam::scam(S ~ s(distance, bs = "mpi"),
                                 data = sSBR_test_high$sSBR_data, family = "poisson")

    pred_low_frag <- predict(scam_low_frag, new_dist, type = "response")
    pred_high_frag <- predict(scam_high_frag, new_dist, type = "response")

    new_dist$diff_sSBR <- pred_high_frag - pred_low_frag

    #' initially, i created the ribbon with the min and max value of the diff_sSBR throughout
    #' all permutations, but in order to solve Felix's request for showing the 95% and 5%
    #' percentile, i now use the 95% of all permutations as the ymax and the 5% as the ymin
    #' and r wants me to transform nd to a matrix instead of data.table
    #' TIP: in a table [x,y] stands for [row,column]

    nd <- as.matrix(nd)

    ribbonbase <- data.frame(matrix(0,
                                    nrow = nrow(nd),
                                    ncol = 5))

    colnames(ribbonbase) <- c("distance","ymin", "ymax", "minline", "maxline")

    #' this saves the distances as one column because geom_ribbon needs an x value

    ribbonbase[,1] <- nd[,1]

    #' looping through all permutations to create the data.set used for my ribbon
    #'
    #' I think the for loop is not necessary, because quantile can deal with vectores
    for (i in 1:nrow(ribbonbase)) {
      ribbonbase[i,2] <- quantile(nd[i,2:101], probs = c(0.05)) # 0.025
      ribbonbase[i,3] <- quantile(nd[i,2:101], probs = c(0.95)) # 0.975
      ribbonbase[i,4] <- min(nd[i,])
      ribbonbase[i,5] <- max(nd[i,])
    }

    #' this is now the plot: we have the blue rarefaction curve of the sample we want to test
    #' against the nullmodel and the transparent green null ribbon, with additional lines for the
    #' min and max value of the nullmodel (i don't know if these are necessary?)

    p <- p + geom_point(data = new_dist, mapping = aes(distance, diff_sSBR), color = "darkturquoise") +
      geom_line(data = new_dist, mapping = aes(distance, diff_sSBR), color = "darkturquoise") +
      geom_ribbon(data = ribbonbase, mapping = aes(ymin = ymin, ymax = ymax, x=distance),
                  fill = "darkgrey", color= "darkgray", alpha = 0.5) +
      geom_line(data = ribbonbase, mapping = aes(distance, minline),
                color = "darkgray", linetype = 2) +
      geom_line(data = ribbonbase, mapping = aes(distance, maxline),
                color = "darkgray", linetype = 2)

    #' printing the result!

    print(p)
    ggsave("plot.jpg")
    print(nd)
    print(ribbonbase)

  }

# run the function, template :
# null_model(dataset = simdat0, frag_level_low = 0.1, frag_level_high = 0.9, timestep = 41, permutations = 50, sample_size = 100)
