rm(list = ls())

# load test database
files = list.files(recursive = T)
split_files = lapply(files, function(f){strsplit(f, split = '/')} )
mask = unlist(lapply(split_files, function(file)file[[1]][[length(file[[1]])]] == 'agdb_h3_small.json'))
db_path = files[which.max(mask)]

library(jsonlite)
db_list = read_json(db_path)



# I assume that all the checker functions are correct, and enforce structure as intended


# agdb.new
db <- agdb.new(db_list)
test_that('agdb.new', {
  expect_equal(agdb.checkagdb(db), T)
  expect_equal(db[[3]]$long,  "A/TOKYO/UT-IMS2-1/2014")
  })

# agdb.id
id <- agdb.id(db)
test_that('agdb.id', {expect_identical(agdb.checkid(id), T)})

# agdb.isolation
isolation <- agdb.isolation('42', '2014', 'egg', 'Cambridge', 'Europe')
test_that('agdb.isolation', {expect_identical(agdb.checkisolation(isolation), T)})

# agdb.genes
HAgene <- agdb.genes(gene = 'HA', sequence = sample(aavalues(), 484))
NAgene <- agdb.genes(gene = 'NA')
genes <- list(HAgene, NAgene)
test_that('agdb.gene',{expect_identical(agdb.checkgenes(genes), T)})

# agdb.alterations
alteration <- list(agdb.alterations('HA', parent_id = 'BNDWF2', substitutions = c('F123S')))
test_that('agdb.alteration', {expect_identical(agdb.checkalterations(alteration), F)})

# agdb.passage
passage = agdb.passage(history = c('MDCK', 'SIAT'), cell = 'egg')
test_that('agdb.passage', agdb.checkpassage(passage))

# agdb.ag
ag = agdb.ag(id = 'ABCDEF',
             strain = 'strain_name',
             long = 'long strain_name',
             aliases = 'alias_name',
             wildtype = T,
             type = 'A',
             subtype = 'H3N2',
             lineage = 'lineage_name',
             isolation = isolation,
             genes = genes,
             parent_id = 'BNDWF2',
             alterations = alterations,
             passage = passage,
             groups = c('mutant', 'gen 1 root'),
             agdb = db,
             .parent = NULL)

test_that('agdb.ag',expect_identical(agdb.checkAG(ag), T))



# agdb.find
test_that('agdb.find', expext_equal(agdb.find(db, wildtype = T), c(3,4)) )


# acdb.search


test_that('acdb.search', {
  search_res = acdb.search(db, wildtype = T)
  expect_equal(search_res, c(3,4))
  expect_equal(search_res[[2]]$long , "A/Puerto Rico/8/1934")
  })


# agdb.extract
test_that('acdb.extract', {
  err = acdb.extract(db, wildtype = T)
  expect_equal(class(err), "try-error")
  expect_identical(acdb.extract(db, long = "A/TOKYO/UT-IMS2-1/2014") , db[[1]])
})
