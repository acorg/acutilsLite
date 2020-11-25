rm(list = ls())

# load test database
files = list.files(recursive = T)
split_files = lapply(files, function(f){strsplit(f, split = '/')} )
mask = unlist(lapply(split_files, function(file)file[[1]][[length(file[[1]])]] == 'agdb_h3_small.json'))
db_path = files[which.max(mask)]

library(jsonlite)
agdb = read.agdb(db_path)



test_that("ag.attribute", {
  expect_identical(ag.attribute(agdb[[1]], attribute = 'isolation'), agdb[[3]]$isolation)
  expect_equal(is.null(ag.attribute(agdb[[1]], attribute = 'isolation', inherit = F)), TRUE)
  expect_identical(ag.attribute(agdb[[1]], attribute = 'isolation', inherit.gene = 'BB'), agdb[[4]]$isolation)
})
