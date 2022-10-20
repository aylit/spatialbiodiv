# Spatially-constrained rarefaction curve --------------

# Diese Funktion errechnet die räumliche und Stichproben-basierte Rarefaction curve
# (auch spatial, sample-based rarefaction curve, sSBR after McGlinn et al. 2019)
# oder auch spatially-constrained rarefaction curve (SCR), after Chiarucci et al. 2009
# Im Unterschied zu McGlinn und Chairucci wird dabei aber nicht die Anzahl an Stichproben
# auf der x-Achse darstellt, sondern die mittlere Distanz zwischen einer
# gegebenen Anzahl von Stichproben
#
# Argumente:
# comm ... community matrix mit den Stichproben als Zeilen und den Arten als Spalten
#          Die Matrix kann Abundanzen enthalten oder auch nur 0 für Art kommt nicht vor
#          und 1 fuer Art kommt vor
# xy_coords ... die x und y Koordinaten der Stichproben
#
# Ausgabe der Funktion:
# data.frame mit den folgenden Spalten:
# n ... Anzahl der Stichproben
# mean_dist ... mittlere Distanz zwischen den n-Stichproben
# mean_S    ... mittlere Artenzahl in den n-Stichproben
# low_S, up_S ... obere und untere Grenze des 95%-Konfidenzintervalls er Artenzahl

SCR <- function(comm,
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
