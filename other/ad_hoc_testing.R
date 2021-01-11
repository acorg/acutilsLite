setwd('~/Desktop/phd/code/acutilsLite/')
devtools::unload('acutilsLite')
devtools::document()
devtools::install()

rm(list = ls())

# load test database
setwd('~/Desktop/phd/code/databases/h3n2/')

library(jsonlite)
library(acutilsLite)
library(tidyverse)

agdb  <- read.agdb("antigens.json")
srdb  <- read.srdb("sera.json")
expdb <- read.expdb("experiments/h3_mutants/results.json")


# Extract the data relating to the single nucs exp
exp <- acdb.search(expdb, name = "Alaska single nuc mutants")[[1]]
exp$results <- exp$results[1:3] # remove two incorrectly included results


# get square titer table
exp %>% exper.merge(expect_repeats = T, threshold = 40) -> mergedTables
mergedTables %>%squaretiters.addNames() -> mergedTables.named

# get long titer table

longTiters.unmerged = exper.toLongTibble(exp) %>% longtiters.addRecords()



# add plotdata
longTiters.unmerged %>%
  longtiters.plotdata() %>%
  longtiters.order() %>%
  mutate(ag_substitutions = agdb.substitutions(.$ag_records)) %>%
  longtiters.splitSubstitutions() ->  longTiters.unmerged

Rprof()
longTiters.merged = longtiters.merge(longTiters.unmerged, columns = c('ag', 'sr'))
Rprof(NULL)
summaryRprof()

titers.getThresholds(longTiters.unmerged$titer)$less

longtiters.merge(filter(longTiters.unmerged, ag %in% c('LNH52L','8AO2PX')))
filter(longTiters.unmerged, ag %in% c('LNH52L','8AO2PX'))
filter(longTiters.merged, ag %in% c('LNH52L','8AO2PX'))

filter(longTiters.merged, ag %in% c('8AO2PX'))

