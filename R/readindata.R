data <- read.csv("sSBR9_rep_1_output_sample.csv", header=TRUE, sep = ",")
data <- filter(data, step == 41) # filtering for timestep 41

source("nullmodel.R")

get_sSBR_dist(dataset = data, frag_level_low = 0.1, frag_level_high = 0.9, 
              sample_size = 100)
get_null_model(dataset = data, frag_level_low = 0.1, frag_level_high = 0.9, 
               sample_size = 100, permutations = 100)
