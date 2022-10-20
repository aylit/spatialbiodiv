#' Distance decay of similarity
#'
#' Estimate pairwise similarities of communities in subplots as
#' function of distance
#'
#' @param comm \code{\link{community}} object
#' @param prop_area Subplot size as proportion of the total area
#' @param n_samples Number of randomly located subplots
#' @param method Choice of (dis)similarity index. See \code{\link[vegan]{vegdist}}
#' @param binary Perform presence/absence standardization before analysis?
#' See \code{\link[vegan]{vegdist}}
#'
#' @return Dataframe with distances between subplot pairs and the respective
#' similarity indices
#'
#' @examples
#' sim_com1 <- sim_thomas_community(100, 10000, sigma = 0.1, mother_points = 2)
#' dd1 <- dist_decay(sim_com1, prop_area = 0.005, n_samples = 20)
#' plot(dd1)
#'
#'@export
#'
dist_decay <- function(comm, xy_coords, method = "bray", binary = F)
{
  require(vegan)

  dist_mat <- stats::dist(xy_coords) # spatial Euclidean distance

  similarity <- 1 - vegan::vegdist(comm, method = method,
                                   binary = binary) # similarity with index defined in vegan
  similarity[!is.finite(similarity)] <- NA

  dat_out <- data.frame(distance = as.numeric(dist_mat),
                        similarity = as.numeric(similarity))

  # order by increasing distance
  out_dat <- dat_out[order(dat_out$distance), ]

  out_pred <- data.frame(distance   = seq(min(out_dat$distance),
                                        max(out_dat$distance),
                                        length = 200),
                         similarity = NA)

  loess1 <-  stats::loess(similarity ~ distance, data = out_dat)

  out_pred$similarity <- stats::predict(loess1, newdata = out_pred)

  return(list(dd_data = out_dat, dd_smooth = out_pred))
}


