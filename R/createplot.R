#' this is now the plot: we have the blue rarefaction curve of the sample we want to test
#' against the nullmodel and the transparent green null ribbon, with additional lines for the
#' min and max value of the nullmodel (i don't know if these are necessary?)
#' 
  ggplot() + 
  geom_point(data = new_dist, mapping = aes(distance, diff_sSBR), color = "darkturquoise") +
  geom_line(data = new_dist, mapping = aes(distance, diff_sSBR), color = "darkturquoise") +
  geom_ribbon(data = ribbonbase, mapping = aes(ymin = ymin, ymax = ymax, x=distance), 
              fill = "darkgrey", color= "darkgray", alpha = 0.5) +
  geom_line(data = ribbonbase, mapping = aes(distance, minline), 
            color = "darkgray", linetype = 2) +
  geom_line(data = ribbonbase, mapping = aes(distance, maxline), 
            color = "darkgray", linetype = 2)