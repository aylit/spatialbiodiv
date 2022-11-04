#' Distance-decay of similarity
#'
#' Estimate pairwise similarities of communities in samples as
#' function of distance
#'
#' @param comm Community matrix with samples in rows and species in columns.
#'             The matrix can contain presence/absence data (coded with 0/1) or
#'             abundance data with counts of individuals.
#'
#' @param xy_coords Two-column table with x-y coordinates of the samples.

#' @param method Choice of (dis)similarity index. See \code{\link[vegan]{vegdist}}
#'
#' @param binary Perform presence/absence standardization before analysis?
#' See \code{\link[vegan]{vegdist}}
#'
#' @return The function returns a list of two dataframes. The first, \code{dd_data} includes
#'         all pairwise distances and their respective similarities.
#'         The second, \code{dd_smooth} includes a fitted GAM
#'        \code{\link[mgcv]{gam}} for the relationship between distance and
#'        similarity.
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

  # Fit model - GAM with monotonously increasing constraint
  gam1 <- mgcv::gam(similarity ~ s(distance), data = out_dat)

  # Predictions - SCAM
  pred <- stats::predict(gam1, out_pred, se = T)

  out_pred$similarity <- pred$fit
  out_pred$simi_low <- pred$fit - 2*pred$se.fit
  out_pred$simi_high <- pred$fit + 2*pred$se.fit

  return(list(dd_data = out_dat, dd_smooth = out_pred))
}


