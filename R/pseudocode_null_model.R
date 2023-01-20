# ------------------------------------------------------------------------------
# sSBR <- function()
# input: samples with coordinates, maybe distances for sSBR evaluation
# output: single sSBR curve

# ------------------------------------------------------------------------------
# get_sSBR_dist <- function
# input: sample_with_coordinates_1, sample_with_coordinates_2,
# output: differenz between two sSBR curves (dependent on distance of course)

# ToDo in function body:
# - calculate distances for sSBR evaluation
# - calculate sSBR1
# - calculate sSBR2
# - calculate difference

# ------------------------------------------------------------------------------
# get_null_interval <- function()
# input: sample_with_coordinates_1, sample_with_coordinates_2
# output: upper and lower interval for null expectation that the two sSBR curves and
# observed difference
# do not differ

# Body
# choose number of permutation N
# Create big table OUT that will include all permutation: N rows, as many columns as distances in sSBR
# combine both data sets (with rbind or similar)
# for i in 1:N {
# permute labels between the two studies
# call get_sSBR_dist
# save distance in OUT
# }
# get 2.5% and 97.5% quantiles from every column of OUT


