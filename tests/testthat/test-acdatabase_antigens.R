rm(list = ls())


# load database

files = list.files(recursive = T)
split_files = lapply(files, function(f){strsplit(f, split = '/')} )
mask = unlist(lapply(split_files, function(file)file[[1]][[length(file[[1]])]] == 'agdb_h3_small.json'))
db_path = files[which.max(mask)]

library(jsonlite)
db_list = read_json(db_path)
db <- agdb.new(db_list)



# agdb.ag

ag = agdb.ag(id = 'ABCDEF',
             strain = 'strain_name',
             long = 'long strain_name',
             aliases = 'alias_name',
             wildtype = T,
             type = 'A',
             subtype = 'H3N2',
             lineage = 'lineage_name',
             isolation = list(location = 'Cambridge', date = '2020', cell = 'egg'),
             genes = list(gene = 'HA', sequence = 'ATCGATCG'),
             parent_id = 'BNDWF2',
             alterations = list(gene = 'HA', substitutions = 'X123Y'),
             passage = c('SIAT', 'MDCK'),
             groups = c('mutant', 'gen 1 root'),
             agdb = db,
             .parent = NULL)

