#
# # fit and predict with SCAM
# sSBR_frag_low <- filter(sSBR_data, fragmentation == "Low")
# scam_log_frag <- scam(S ~ s(distance, bs = "mpi"),
#                       data = sSBR_frag_low, family = poisson)
# newdat1 <- data.frame(distance = seq(min(sSBR_frag_low$distance),
#                                      max(sSBR_frag_low$distance),
#                                      length = 100))
# newdat1$predS <- predict(scam_log_frag, newdata = newdat1, type = "response")
#
# ggplot(sSBR_frag_low, aes(distance, S)) +
#   geom_point() +
#   geom_smooth(method = "gam") +
#   geom_line(aes(x = distance, y = predS), data = newdat1,
#             color = "red")
#
# # Constrained GAM from
# # https://stackoverflow.com/questions/66648829/force-gam-model-fit-to-be-monotonic-and-go-through-a-fixed-point-x0-y0-with-r
#
# # ----
# library(mgcv)
# set.seed(1)
# x <- sort(runif(100) * 4 - 1)
# f <- exp(4*x)/(1+exp(4*x))
# y <- f + rnorm(100) * 0.1
# dat <- data.frame(x=x, y=y)
#
# k <- 13
#
#
# #fit0 <- gam(y ~ s(x, k = k, bs = "cr"), data = dat)
# fit0 <- gam(y ~ s(x, bs = "cr"), data = dat)
#
# # predict from unconstrained GAM fit
# newdata <- data.frame(x = seq(-1, 3, length.out = 1000))
# newdata$y_pred_fit0 <- predict(fit0, newdata = newdata)
#
# # Show regular spline fit (and save fitted object)
# #f.ug <- gam(y~s(x,k=k,bs="cr"))
# f.ug <- gam(y~s(x,bs="cr"))
#
#
# # explicitly construct smooth term's design matrix
# #sm <- smoothCon(s(x,k=k,bs="cr"),dat,knots=NULL)[[1]]
# sm <- smoothCon(s(x,k=k,bs="cr"),dat,knots=NULL)[[1]]
# # find linear constraints sufficient for monotonicity of a cubic regression spline
# # it assumes "cr" is the basis and its knots are provided as input
# F <- mono.con(sm$xp)
#
# G <- list(
#   X=sm$X,
#   C=matrix(0,0,0),  # [0 x 0] matrix (no equality constraints)
#   sp=f.ug$sp,       # smoothing parameter estimates (taken from unconstrained model)
#   p=sm$xp,          # array of feasible initial parameter estimates
#   y=y,
#   w= dat$y * 0 + 1  # weights for data
# )
# G$Ain <- F$A        # matrix for the inequality constraints
# G$bin <- F$b        # vector for the inequality constraints
# G$S <- sm$S         # list of penalty matrices; The first parameter it penalizes is given by off[i]+1
# G$off <- 0          # Offset values locating the elements of M$S in the correct location within each penalty coefficient matrix.  (Zero offset implies starting in first location)
#
# p <- pcls(G);       # fit spline (using smoothing parameter estimates from unconstrained fit)
#
# # predict
# newdata$y_pred_fit2 <- Predict.matrix(sm, data.frame(x = newdata$x)) %*% p
# # plot
# plot(y ~ x, data = dat)
# lines(y_pred_fit0 ~ x, data = newdata, col = 2, lwd = 2)
# lines(y_pred_fit2 ~ x, data = newdata, col = 4, lwd = 2)
