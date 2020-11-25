rm(list = ls())

# load test database
files = list.files(recursive = T)
split_files = lapply(files, function(f){strsplit(f, split = '/')} )
mask = unlist(lapply(split_files, function(file)file[[1]][[length(file[[1]])]] == 'agdb_h3_small.json'))
db_path = files[which.max(mask)]

library(jsonlite)
agdb = read.agdb(db_path)


# gene_inheritance

giHA = gene_inheritance(agdb[[1]], 'HA')
giNA = gene_inheritance(agdb[[1]], 'NA')
giBB = gene_inheritance(agdb[[1]], 'backbone')

test_that("gene_inheritance", {
  expect_identical(giHA, giNA)
  expect_identical(giBB[[1]], agdb[[4]])
  expect_identical(giBB[[2]], giHA[[2]])
  expect_equal(giBB[[1]]$long, "A/Puerto Rico/8/1934")
})
