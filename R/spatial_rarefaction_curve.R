#' Spatially-explicit sample-based species rarefaction curve
#'
#' This function calculates the spatial and sample-based species rarefaction curve
#' (sSBR after McGlinn et al. 2019), also called spatially-constrained
#' rarefaction curve (SCR) (Chiarucci et al. 2009).
#'
#' @param comm Community matrix with samples in rows and species in columns.
#'             The matrix can contain presence/absence data (coded with 0/1) or
#'             abundance data with counts of individuals.
#'
#' @param xy_coords Two-column table with x-y coordinates of the samples.
#'
#' @details In contrast to McGlinn et al. 2019, this function uses accumulated distances
#' among samples as x-axis and not the accumulated number of samples.
#'
#' @references McGlinn, D. J. et al. 2019. Measurement of Biodiversity (MoB):
#'  A method to separate the scale-dependent effects of species abundance distribution,
#'   density, and aggregation on diversity change. - Methods in Ecology and Evolution 10: 258–269.
#'
#'  Chiarucci, A. et al. 2009. Spatially constrained rarefaction:
#'  incorporating the autocorrelated structure of biological communities into
#'   sample-based rarefaction. - Community Ecology 10: 209–214.
#'
#' @return The function returns a list of two dataframes. The first, \code{sSBR_data} includes
#'         all accumulated distances and their respective accumulated species richness for
#'         the same number of rarefaction curves as there are samples in the data set.
#'         Each samples is used as a starting point and then samples are accumulated
#'         by distances to the k nearest neighbour.
#'         The second dataframe, \code{ssBR_smooth} includes a non-linear smoother
#'        \code{\link[stats]{loess}} for the relationship between cimulative distance and
#'        cumulative species richness.
#' @export
#'

sSBR <- function(comm,
                xy_coords) {
  comm <- (comm > 0) * 1 # change to presence-absence matrix
  n <- nrow(comm)
  # drop species with no observations
  comm <- comm[ , colSums(comm) > 0]
  scr_mat <- matrix(0, n, n)
  dist_mat <- matrix(0, n, n)
  pair_dist <- as.matrix(stats::dist(xy_coords))

  for (i in 1:n) {
    dist_to_site <- pair_dist[i, ]
    # Shuffle plots, so that tied grouping is not biased by original order.
    new_order <- sample(1:n)
    dist_new <- dist_to_site[new_order]
    new_order <- new_order[order(dist_new)]
    # Move focal site to the front
    new_order <- c(i, new_order[new_order != i])
    comm_ordered <- comm[new_order, ]
    # 1 for absence, 0 for presence
    comm_bool <- as.data.frame((comm_ordered == 0) * 1)
    rich <- cumprod(comm_bool)

    scr_mat[i, ] <- as.numeric(ncol(comm) - rowSums(rich))
    dist_mat[i, ] <- dist_to_site[order(dist_to_site)]
  }

  # This aggregation is likely inappropriate, because it averages distances
  # of pairs of locations that can have very different distances

  # out_mean <- data.frame(n = 1:nrow(scr_mat),
  #                        mean_dist = colMeans(dist_mat),
  #                        mean_S    = colMeans(scr_mat),
  #                        low_S     = apply(scr_mat, 2, quantile, prob = 0.025),
  #                        up_S      = apply(scr_mat, 2, quantile, prob = 0.975))

  out_dat <- data.frame(id       = rep(1:nrow(comm), times = nrow(comm)),
                        distance = as.vector(dist_mat),
                        S        = as.vector(scr_mat))

  out_dat <- out_dat[order(out_dat$id, out_dat$distance),]

  loess1 <-  stats::loess(S ~ distance, data = out_dat)

  out_pred <- data.frame(distance = seq(min(out_dat$distance),
                                        max(out_dat$distance),
                                        length = 200),
                         S        = NA)

  out_pred$S <- stats::predict(loess1, newdata = out_pred)

  return(list(sSBR_data = out_dat, sSBR_smooth = out_pred))

  #return(list(sSBR_mean = out_mean, dist_mat = dist_mat, S_mat = scr_mat))
}
