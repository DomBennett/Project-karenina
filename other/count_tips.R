# Count number of added fossil tips to trees
# (In response to reviewer 3 comment)

# Lib ----
library(treeman)

# Vars ----
input_dir <- '1_pin'
parent <- 'Mammalia'
pinfiles <- list.files(file.path(input_dir, paste0(parent, '_rndm')))
pinfiles <- file.path(input_dir, paste0(parent, '_rndm'), pinfiles)
mammal_tree <- readTree(file.path('0_data', 'mammalia.tre'))
ntips_before <- mammal_tree@ntips

# Count
fossil_tips_count <- function(pinfile) {
  load(pinfile)
  res <- tree@ntips - ntips_before
  res
}
fossil_counts <- vapply(X = pinfiles, FUN = fossil_tips_count,
                        FUN.VALUE = integer(1))
mean(fossil_counts)
sd(fossil_counts)