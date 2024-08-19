library(tidyverse)
library(here)

df <- read_csv(here('day2_walkthrough_script/sfo_trees.csv'))

set.seed(42)
df_sample <- slice_sample(df, prop = .1)

write_csv(df_sample, here('day2_walkthrough_script/sfo_trees_sample.csv'))