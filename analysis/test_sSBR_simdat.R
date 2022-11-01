# Read data
library(tidyverse)
#library(broom)
library(modelr)
library(scam)

# Read and prepare the data ----------------------------------------------------

simdat1 <- read.table("data_raw/GeDo_test5_rep_10_output_sample.txt", header = TRUE, sep = ",")
simdat1 <- simdat1 %>% select(-X)

frag_high <- simdat1 %>%
  filter(step == 31 & fragmentation == 0.5)

frag05_spec <- frag_high %>%
  select(sp_1:sp_200)

xy_frag05 <-  frag_high %>%
  select(loc_x, loc_y)

ggplot(xy_frag05, aes(loc_x, loc_y)) +
  geom_point()

# Spatial sample-based rarefaction curves (sSBR) -------------------------------

sSBR_frag05 <- sSBR(comm = frag05_spec, xy_coords = xy_frag05)

ggplot(sSBR_frag05$sSBR_data, aes(distance, S, color = as.factor(id))) +
  geom_point() +
  geom_line()

# Single SCAM model for all curves with CI band --------------------------------

# Fit model
scam1 <- scam(S ~ s(distance, bs = "mpi"),
                    data = sSBR_frag05$sSBR_data, family = "poisson")

# New data
dat_new <- data_grid(sSBR_frag05$sSBR_data,
                     distance = seq_range(distance, 100))

# Predictions - SCAM
pred <- predict(scam1, dat_new, se = T, type = "response")

dat_new$pred <- pred$fit
dat_new$pred_se <- pred$se.fit

dat_new <- dat_new %>%
  mutate(low1 = pred - 2*pred_se,
         up1 = pred + 2*pred_se)

# # Predictions LOESS
# loess1 <-  stats::loess(S ~ distance, data = sSBR_frag05$sSBR_data)
#
# pred_loess <- predict(loess1, dat_new, se = T)
#
# dat_new$pred_loess <- pred_loess$fit
# dat_new$pred_se_loess <- pred_loess$se.fit
#
# dat_new <- dat_new %>%
#   mutate(low2 = pred_loess - 2*pred_se_loess,
#          up2 = pred_loess + 2*pred_se_loess)


geom_line(aes(x = distance, y = pred_resp),
          color = "black", size = 1, data = dat_new) +


# single sSBR
#sSBR1 <- filter(sSBR_frag05$sSBR_data, id == 1)


#newdat1 <- data.frame(distance = seq(min(sSBR1$distance),
#                                     max(sSBR1$distance),
#                                     length = 100))
#newdat1$predS <- predict(scam1, newdata = newdat1, type = "response")
#plot(predS ~ distance, data = newdat1)

# One scam for each single curve
# After
# https://r4ds.had.co.nz/many-models.html#many-models

by_id <- sSBR_frag05$sSBR_data %>%
  group_by(id) %>%
  nest()

scam_model <- function(df){
  scam(S ~ s(distance, bs = "mpi"), data = df, family = poisson)
}


# pred_scam <- function(df){
#   mod1 <- scam(S ~ s(distance, bs = "mpi"), data = df, family = poisson)
#   newdat1 <- data_grid(df, distance = seq_range(distance, 100))
#   newdat1 <- newdat1 %>% add_predictions(mod1, type = "response")
#   return(newdat1)
# }

# models <- map(by_id$data, scam_model)

by_id <- by_id %>%
  mutate(model = map(data, scam_model))

# map(by_id$model, function(model){predict(model,
#                                          newdata = data.frame(distance = seq_range(sSBR_frag05$sSBR_data, 100))
#                                          type = response)}
# )


#pred_id <- unnest(by_id_pred, pred)

# ggplot(sSBR_frag05$sSBR_data, aes(distance, S, color = as.factor(id))) +
#   #geom_point() +
#   #geom_line(linetype = 2) +
#   geom_line(data = pred_id, aes(y = pred))


#by_id_pred <- map(by_id$data, pred_scam)

# Create table for predictions
dat_new <- data_grid(sSBR_frag05$sSBR_data, id, distance = seq_range(distance, 100))

#by_id_2 <- by_id %>% left_join(dat_new)

dat_new_by_id <- dat_new %>%
  group_by(id) %>%
  nest()

dat_new_by_id <- dat_new_by_id %>% left_join(select(by_id, id, model))

dat_new_by_id <- dat_new_by_id %>%
  mutate(predS = map2(data, model, function(x,y){predict(y, newdata = x, type = "response")}))

dat_new_pred <- dat_new_by_id %>%
  select(id, data, predS) %>%
  unnest(c(data, predS))

ggplot(sSBR_frag05$sSBR_data, aes(distance, S, color = as.factor(id))) +
  geom_point() +
  geom_line(linetype = 2) +
  geom_line(data = dat_new_pred, aes(y = predS))

# also try:
# models <- mtcars %>%
#   split(.$cyl) %>%
#   map(function(df) lm(mpg ~ wt, data = df))

# models <- mtcars %>%
#   split(.$cyl) %>%
#   map(~lm(mpg ~ wt, data = .))

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
