##########################
#
#    HELPER
#
##########################

#
# If databases are not supplied, search for agdb / srdb / expdb in global frame
#


get_agdb <- function(){
  get("agdb", envir = parent.frame(n = 2))
}

get_srdb <- function(){
  get("srdb", envir = parent.frame(n = 2))
}

get_expdb <- function(){
  get("expdb", envir = parent.frame(n = 2))
}



#' List to vector
#'
#' Converts list of length 1 vectors to vector. Throws error if an element has length > 1
#'
#' @param x list
#'
#' @return vector
#' @export
#'
#' @examples
collate <- function(x, null = NA){
  sapply(x, function(n){
    if(is.null(n))    return(null)
    if(length(n) > 1) stop("Length > 1")
    n
  })
}



##########################
#
#    EXTRACT INFO
#
##########################

#
# From single record
#


#' Get gene attribute
#'
#' Find an attribute value for an antigen. By default this is set to traverse the inheritance list for 'HA' to find the attribute value but the gene can be set to 'NA' or 'backbone', or this behaviour turned off.
#'
#' @param ag environment: an antigen
#' @param attribute char: the name of the desired attiribute
#' @param inherit: bool: whether to inherit attributes from parents
#' @param inherit.gene char: which gene to inherit along
#'
#' @return attribute value
#' @export
#'
#' @examples
ag.attribute <- function(
  ag,
  attribute,
  inherit = TRUE,
  inherit.gene = "HA"
){

  # Assemble antigen query
  query <- paste0(
    "[['",
    paste0(
      strsplit(attribute, ".", TRUE)[[1]],
      collapse = "']][['"
    ),
    "']]"
  )

  # Get antigen inheritance
  if(inherit){
    ag_inheritance <- ag.inheritance(
      ag   = ag,
      gene = inherit.gene
    )
  } else {
    ag_inheritance <- list(ag)
  }

  # Get attribute
  for(ag in rev(ag_inheritance)){
    value <- eval(parse(text = paste0("ag", query)))
    if(!is.null(value)) break
  }

  # Return the value
  value

}


#' Get serum attribute
#'
#' Find an attribute value for an serum By default this is set to check parent sera if attribute not found at lower level.
#'
#' @param sr environment: a serum
#' @param attribute char: the name of the desired attiribute
#' @param inherit: bool: whether to inherit attributes from parents
#'
#' @return attribute value
#' @export
#'
#' @examples
sr.attribute <- function(sr, attribute, inherit = T){
  return(ag.attribute(sr, attribute, inherit, inherit.gene = 'serum'))
}


#' Get parent
#'
#' @param ag environment: an antigen
#'
#' @return list
#' @export
#'
#' @examples
ag.parent <- function(ag){
  parent <- ag$.parent
  if(length(parent) == 0) NULL
  else                    parent
}

#' Get parent
#'
#' @param sr environment: a serum
#'
#' @return list
#' @export
#'
#' @examples
sr.parent <- ag.parent



#' Get inheritance list for a gene.
#'
#' Get a list of parent antigens for a gene in an antigen. The list will follow the parent_id's for the specified component ('HA', 'NA', or 'backbone') until reaching the root.
#'
#' @param ag environment: an antigen
#' @param gene char: which gene inheritance list is required (must be 'HA', 'NA', or 'backbone)
#'
#' @return list
#' @export
#'
#' @examples
ag.inheritance <- function(
  ag,
  gene,
  inheritance = list()){

  # Check input
  if(!inherits(ag, "acdatabase.entry")){
    stop(sprintf("Class of input must inherit from 'acdatabase.entry' but was '%s'", paste(class(ag), collapse = ", ")))
  }

  # Get the ha alterations
  ha_alteration <- match(gene, sapply(ag[["alterations"]], function(x) x$gene))

  # Add antigen to inheritance
  inheritance <- c(list(ag), inheritance)

  # Look for inheritance
  if(!is.na(ha_alteration) && !is.null(ag[["alterations"]][[ha_alteration]][["parent_id"]])){
    ag.inheritance(
      ag            = ag[["alterations"]][[ha_alteration]][[".parent"]],
      gene          = gene,
      inheritance   = inheritance
    )
  } else if(!is.null(ag$parent_id)){
    ag.inheritance(
      ag            = ag[['.parent']],
      gene          = gene,
      inheritance   = inheritance
    )
  } else{
    inheritance
  }

}



#' Get inheritance list for a serum.
#'
#' Get a list of parent sera (previous sera from same individual).
#'
#' @param sr environment: a serum
#'
#' @return list
#' @export
#'
#' @examples
sr.inheritance = function(sr) return(ag.inheritance(sr, gene = 'serum'))

#' Get children
#'
#' @param ag environment: an antigen
#'
#' @return list
#' @export
#'
#' @examples
ag.children <- function(ag, agdb){
  agdb_parents <- agdb%$%parent_id
  agdb[which(agdb_parents == ag$id)]
}

#' Get children
#'
#' @param sr environment: a serum
#'
#' @return list
#' @export
#'
#' @examples
sr.children <- function(sr, srdb) return(ag.children(sr, srdb))


#' Get long name for an antigen
#'
#' @param ag environment: an antigen
#' @param agdb list: an agdb
#'
#' @return list
#' @export
#'
#' @examples
ag.long <- function(ag, agdb = get_agdb()){
  acdb.Long(list(ag))
}

#' Get standardised long name for an antigen
#'
#' @param ag environment: an antigen
#' @param agdb list: an agdb
#'
#' @return list
#' @export
#'
#' @examples
ag.standardlong <- function(ag, agdb = get_agdb()){
  long = ag.long(ag)
  return(standardizeStrainNames(long))
}


#' Construct short name for an antigen
#'
#' @param ag environment: an antigen
#' @param agdb list: an agdb
#'
#' @return list
#' @export
#'
#' @examples
ag.short <- function(ag, agdb = get_agdb()){
  agdb.Short(list(ag))
}


#' Get all descendants
#'
#' Gets all descendants of an antigen (backbone-wise), preserving tree of relationships (compare to `?ag.descendants`)
#'
#' @param ag environment: an antigen
#'
#' @return list
#' @export
#'
#' @examples
ag.descendentsTree <- function(ag, agdb){
  ag_children <- ag.children(ag, agdb)
  list(ag = ag, children = lapply(ag_children, ag.descendentsTree, agdb = agdb))
}

#' Get all descendants
#'
#' Gets all descendants of an antigen (backbone-wise), flattening to list (compare to `?ag.descendantsTree`)
#'
#' @param ag environment: an antigen
#'
#' @return list
#' @export
#'
#' @examples
#' @export
ag.descendents <- function(ag, agdb){
  ag_tree <- ag.descendentsTree(ag, agdb)
  descendents <- list()
  tree_children <- function(ag_tree){
    descendents <<- c(descendents, ag_tree$ag)
    lapply(ag_tree$children, tree_children)
  }
  tree_children(ag_tree)
  descendents[-1]
}



#' Calculate gene sequence
#'
#'
#' @param ag environment: an antigen
#' @param gene char: which gene to get sequence for
#'
#' @return list
#' @export
#'
#' @examples
ag.sequence <- function(
  ag,
  gene = 'HA',
  agdb = acutils:::get_agdb(),
  inherit = T
){


  if (!inherit){
    ha_sequence <- match(gene, sapply(ag[["genes"]], function(x) x$gene))

    # Add any substitutions listed
    if(!is.na(ha_sequence)){
      sequence <- ag[["genes"]][[ha_sequence]][["sequence"]]
    } else {
      sequence <- NULL
    }

    # Return the substitutions
    return(sequence)
  }

  # Get the gene substitutions
  inheritance <- ag.inheritance(
    ag   = ag,
    gene = gene
  )

  # Go through the inheritance tree
  subs <- list()
  for(ag in inheritance){
    subs     <- c(subs, list(ag.substitutions(ag, gene)))
    sequence <- ag.sequence(ag, gene)
    if(!is.null(sequence)) break
  }
  if(is.null(sequence)) return(NULL)

  # Apply any modifications to the sequence
  for(sub in unlist(subs)){

    # Fetch positions and substitutions
    positions <- as.numeric(substr(sub, 2, nchar(sub)-1))
    subfrom   <- substr(sub, 1, 1)
    subto     <- substr(sub, nchar(sub), nchar(sub))

    # Do some checks
    if(sum(is.na(positions)) > 0){
      stop("Did not understand substitutions")
    }
    if(sum(duplicated(positions)) > 0){
      stop("Multiple substitutions at same site")
    }

    # Modify the sequence
    for(x in seq_along(positions)){

      if(substr(sequence, positions[x], positions[x]) != subfrom[x]){
        warning(
          sprintf(
            "Sequence of parent does not match the substitution designated, substitution was %s but parent has %s at position %s.",
            sub, substr(sequence, positions[x], positions[x]), as.character(positions[x])
          )
        )
      }

      substr(sequence, positions[x], positions[x]) <- subto[x]

    }

  }

  # Return the sequence
  sequence
}

#' Get gene substituions for antigen
#'
#' Gets all gene substitutions since the root in specified gene
#'
#' @param ag environment: an antigen
#' @param gene char
#' @param inherit bool: whether to check parents of gene for substitutions
#'
#' @return list
#' @export
#'
#' @examples
ag.substitutions <- function(ag, gene, substitutions = list(), inherit = TRUE){

  # Return substitutions if not inheriting
  if(!inherit){
    ag_inheritance <- list(ag)
  } else {
    ag_inheritance <- ag.inheritance(ag, gene = "HA") # should this be gene = gene?
  }

  # Get the gene alterations
  lapply(ag_inheritance, function(ag){
    gene_alteration <- match(gene, sapply(ag[["alterations"]], function(x) x$gene))
    if(!is.na(gene_alteration)){
      ag[["alterations"]][[gene_alteration]][["substitutions"]]
    } else {
      NULL
    }
  })

}



#' Get antigen cluster
#'
#' Get an antigen's cluster from id / long
#'
#' @param identifier char: an id, long, etc
#' @param agdb
#' @param from: char: which attribute to search with (long, id, etc)
#'
#' @return char
#' @export
#'
#' @examples
agid.cluster <- function(identifier, agdb, from = 'id'){

  if (length(identifier) == 6 & from != 'id') print('ARE YOU SURE THIS IS NOT AN id?')
  ag = agdb[[which(acdb.slice(agdb, from) == identifier)]]
  if (is.null(ag)) stop('ag not found')
  ag$meta$cluster$name

}

#' Get antigen cluster
#'
#' Get an antigen's cluster from ag
#'
#' @param ag
#'
#' @return char
#' @export
#'
#' @examples
ag.cluster <- function(ag){ag.attribute(ag, 'meta.cluster.name')}

#' Get serum cluster
#'
#' Gets the cluster of a serum's homologous antigen
#'
#' @param identifier char: an id, long, etc
#' @param agdb
#' @param srdb
#' @param from: char: which attribute to search with (long, id, etc)
#'
#' @return char
#' @export
#'
#' @examples
srid.cluster <- function(identifier, agdb, srdb, from = 'id'){

  srmatch <- which(acdb.slice(srdb, from) == identifier)
  if(length(srmatch) == 0) srmatch <- which(sapply(srdb%$%aliases, function(x) srname %in% x ))
  if(length(srmatch) == 0) stop(identifier, " not found in database")
  sr <- srdb[[srmatch]]
  ag <- agdb[[which(agdb%$%id == sr$strain_id)]]
  ag$meta$cluster$name

}


#' Get serum cluster
#'
#' Gets the cluster of a serum's homologous antigen
#'
#' @param identifier char: an id, long, etc
#' @param agdb
#' @param srdb
#' @param from: char: which attribute to search with (long, id, etc)
#'
#' @return char
#' @export
#'
#' @examples
sr.cluster <- function(sr, agdb){

  homologous_ag = acdb.getIDs(agdb, sr$strain_id)
  return(ag.cluster(homologous_ag))

}

#' Get antigen root
#'
#'
#' @param ag
#' @param gene
#'
#' @return char
#' @export
#'
#' @examples
ag.root <- function(ag, gene = 'HA'){
  ag.inheritance(ag, gene)[[1]]
}


#' Get antigen root clade
#'
#'
#' @param ag
#' @param gene
#'
#' @return char
#' @export
#'
#' @examples
ag.rootClade <- function(ag, gene = 'HA'){
  root = ag.root(ag, gene)
  clades = intersect(root$groups, h3_clade_grouping_heirarchy)
  if (length(clades) == 0) return(NA)
  first_matching_clade = h3_clade_grouping_heirarchy[[min(which(which(h3_clade_grouping_heirarchy %in% clades)))]]
  return(first_matching_clade)
}



#' Check if ag/sr is in groups
#'
#'
#' @param ag
#' @param ... groups to check
#'
#' @return char
#' @export
#'
#' @examples
ac.inGroups <- function(ac, ...){
  groups <- c(...)
  return(sum(!groups %in% ac[["groups"]]) == 0)
}





#
# From acdb
#


#' Get named attribute from acdb
#'
#' @param ids char
#' @param attribute char
#' @param inherit
#' @param inherit.gene
#' @param acdb
#'
#' @return char
#' @export
#'
#' @examples
acdb.idToAttr <- function(ids, attribute, inherit = T, inherit.gene = 'HA', acdb = get_agdb()){
  acs = acdb.getIDs(ids, acdb)
  attributes = lapply(acs, ag.attribute, attribute, inherit, inherit.gene)
  return(acdb.slice(acs, attribute))
}

#' Get long names from acdb
#'
#' @param ids char
#' @param acdb
#'
#' @return char
#' @export
#'
#' @examples
acdb.nameIDs <- function(ids, acdb, name = 'long'){
  if (name == 'long') return(acdb.idToAttr(ids, 'long', inherit = F, acdb = acdb))
  if (name == 'short') return(agdb.Short(acdb.getIDs(acdb, ids)))
  else stop('name must be "long" or "short"')
  }


#' Get ag clusters from acdb
#'
#' @param agdb
#'
#' @return char
#' @export
#'
#' @examples
agdb.clusters <- function(agdb){
  clusters.toFactor(
    collate(acdb.applyFunction(agdb, function(x){
      ag.attribute(x, "meta.cluster.name")
    }))
  )
}


#' Get sera clusters
#'
#' Gets the clusters of sera's homologous antigens.
#'
#' @param sera vector of identifiers
#' @param agdb
#' @param srdb
#' @param from: char: which attribute to search with (long, id, etc)
#'
#' @return char
#' @export
#'
#' @examples
srdb.clusters <- function(sera, agdb, srdb){

  sera <- as.factor(sera)
  level_clusters <- sapply(levels(sera), sr.cluster, agdb = agdb, srdb = srdb)
  level_clusters[as.numeric(sera)]

}




#' Get homologous sera
#'
#' Gets homologous sera for all antigens in an agdb by matching ids.
#'
#' @param agdb
#' @param srdb
#'
#' @return list
#' @export
#'
#' @examples
agdb.homologousSera <- function(agdb, srdb){

  # Get ids of antigens in database
  antigen_ids <- agdb%$%id

  # Get matching antigens
  srags_ids <- srdb%$%strain_id
  srdb[match(antigen_ids, srags_ids)]

}

#' Get homologous sera
#'
#' Gets homologous antigens for all sera in an srdb by matching ids.
#'
#' @param srdb
#' @param agdb
#'
#' @return list
#' @export
#'
#' @examples
srdb.homologousAntigens <- function(srdb, agdb){

  # Get ids of antigens in database
  antigen_ids <- srdb%$%strain_id

  # Get matching antigens
  agdb_ids <- agdb%$%id
  agdb[match(antigen_ids, agdb_ids)]

}


# Functions to check if antigens and sera are within certain groups
#' @export
acdb.inGroups <- function(acdb, ...){
  groups <- c(...)
  vapply(acdb, ac.inGroups, logical(1))
}




#' Generate short name from antigen entries
#'
#' @param ags
#'
#' @return Returns a vector of short names
#' @export
#'
agdb.Short <- function(agdb){
  collate(acdb.applyFunction(agdb, function(ag){

    # Get shortened place
    place <- ag.attribute(ag, "isolation.location")
    place_short <- names(place_abvs)[match(toupper(place), place_abvs)]
    if(is.na(place_short)) place_short <- place
    place_short <- toupper(place_short)

    if(place_short == "HANOI")           place_short <- "HN"
    if(place_short == "NEW_YORK")        place_short <- "NY"
    if(place_short == "SOUTH_AUSTRALIA") place_short <- "SA"
    if(place_short == "NEW_CALEDONIA")   place_short <- "NW"
    if(place_short == "FUKUOKA")         place_short <- "FK"
    if(place_short == "HONG_KONG")       place_short <- "HK"
    if(place_short == "NEW_CASTLE")      place_short <- "NC"
    if(place_short == "SOUTH_AUCKLAND")  place_short <- "SAU"

    # Get shortened year
    year_short <- ag.attribute(ag, "isolation.date")
    year_short <- substr(year_short, nchar(year_short)-1, nchar(year_short))

    # Get egg passage details
    if(passage_includes_egg(ag)) year_short <- paste0(year_short, "E")

    # Assemble the short name
    name_short <- paste(
      c(
        place_short,
        ag.attribute(ag, "isolation.id"),
        year_short
      ),
      collapse = "/"
    )

    # Add any alterations
    subs <- unlist(ag.substitutions(ag, gene = "HA"))
    name_short <- paste(c(name_short, subs), collapse = " ")

    # Compare overall inheritance with HA inheritance
    ag_base <- ag.inheritance(ag, gene = 'BB')[[1]]
    if(ag_base$id == "BNDWF2"){
      name_short <- paste("PR8", name_short)
    }

    # Return the name
    name_short

  }))
}


#' Select long name from antigen entries
#'
#' @param ags
#'
#' @return Returns a vector of long names
#' @export
#'
acdb.Long <- function(acdb){
  collate(acdb.applyFunction(acdb, function(x){ x[["long"]] }))
}




#' Get substitutions from agdb
#'
#' @param agdb
#'
#' @return
#' @export
agdb.substitutions <- function(agdb, inherit = T){
  acdb.applyFunction(agdb, function(ag){ unlist(ag.substitutions(ag, gene = "HA", inherit)) })
}

#' Get sequences from agdb
#'
#' @param agdb
#'
#' @return
#' @export
agdb.sequences <- function(agdb, gene = 'HA', inherit = T){
  acdb.applyFunction(agdb, function(ag){ ag.sequence(ag, gene, agdb, inherit) })
}

#' @export
agdb.year <- function(agdb){
  collate(acdb.applyFunction(agdb, function(x){
    as.numeric(ag.attribute(x, "isolation.date"))
  }))
}


#' @export
agdb.root <- function(agdb, gene){
  acdb.applyFunction(agdb, function(ag){
    ag.root(ag, gene)
  })
}


#
# Stuff to do with experiment - may be moved later
#



#' @export
acdb.agExperiments <- function(ags, expdb, include_passage_variants, agdb){
  acdb.applyFunction(ags, function(x){ ag.experiments(x, expdb, include_passage_variants, agdb) })
}


# Functions to check if antigens and sera are present in certain experimental results
#' @export
ag.inResults <- function(agdb, exp){
  agids <- collate(agdb%$%id)
  lapply(agids, function(id){
    which(vapply(exp$results, function(result){
      id %in% result$antigen_ids
    }, logical(1)))
  })
}

#' @export
sr.inResults <- function(srdb, exp){
  srids <- collate(srdb%$%id)
  lapply(srids, function(id){
    which(vapply(exp$results, function(result){
      id %in% result$serum_ids
    }, logical(1)))
  })
}


#' @export
acdb.nameTableIDs <- function(
  titertable,
  agdb = get_agdb(),
  srdb = get_srdb(),
  append_names = FALSE
){
  aglong <- acdb.nameIDs(rownames(titertable), agdb)
  srlong <- acdb.nameIDs(colnames(titertable), srdb)
  if(append_names){
    aglong <- paste(rownames(titertable), aglong, sep = ": ")
    srlong <- paste(colnames(titertable), srlong, sep = ": ")
  }
  rownames(titertable) <- aglong
  colnames(titertable) <- srlong
  titertable
}


#' @export
acdb.expName <- function(exps){
  collate(acdb.applyFunction(exps, function(x){ x[["name"]] }))
}

#' @export
consensus <- function(x){
  x[x == "*"] <- NA
  summary_x <- summary(as.factor(x))
  names(summary_x)[which.max(summary_x)]
}

#' @export
seq_consensus <- function(x){
  x <- x[!sapply(x, is.null)]
  if(length(x) == 1) return(x[[1]])
  x <- unlist(lapply(x, strsplit, split = ""), recursive = FALSE)
  paste(apply(do.call(rbind, x), 2, consensus), collapse = "")
}

#' @export
expdb.ags <- function(expdb){
  unlist(
    lapply(expdb, function(exp){
      lapply(exp$results, function(result){
        result$antigen_ids
      })
    })
  )
}

#' @export
expdb.sr <- function(expdb){
  unlist(
    lapply(expdb, function(exp){
      lapply(exp$results, function(result){
        result$serum_ids
      })
    })
  )
}

#' @export
ag.experiments <- function(
  ag,
  expdb,
  include_passage_variants = FALSE,
  agdb
){

  # Work out the id to search on
  if(include_passage_variants){
    ag_id <- passage_inheritance(ag)[[1]]$id
  } else {
    ag_id <- ag$id
  }

  # Keep a list of matching experiments
  ag_exps <- list()

  # Work out the experiment ids to search on
  for(exp in expdb){
    exp_ag_ids <- expdb.ags(list(exp))
    if(include_passage_variants){
      exp_ags <- acdb.getIDs(agdb, exp_ag_ids)
      exp_ag_ids <- vapply(
        exp_ags,
        function(ag){
          passage_inheritance(ag)[[1]]$id
        },
        character(1)
      )
    }

    # Add experiment to list if antigen is found
    if(ag_id %in% exp_ag_ids){
      ag_exps <- c(ag_exps, exp)
    }
  }

  # Return the list of matching experiments
  ag_exps

}

##########################
#
#    SEARCHING
#
##########################



#' Subset database based on IDs
#'
#' Throws error if some IDs not found
#'
#' @param db
#' @param ids
#'
#' @return list
#' @export
#'
#' @examples
acdb.getIDs <- function(db, ids){

  dbids <- unlist(db%$%ids)
  if(sum(!ids %in% dbids) > 0) stop("IDs not found in database")
  db[match(ids, dbids)]

}

#' Subset antigen database based on IDs
#'
#' Throws error if some IDs not found
#'
#' @param db
#' @param ids
#'
#' @return list
#' @export
#'
#' @examples
agdb.getIDs <- function(ids, agdb){
  if(missing(agdb)) agdb <- get_agdb()
  acdb.getIDs(agdb, ids)
}

#' Subset serum database based on IDs
#'
#' Throws error if some IDs not found
#'
#' @param db
#' @param ids
#'
#' @return list
#' @export
#'
#' @examples
srdb.getIDs <- function(ids, srdb){
  if(missing(srdb)) srdb <- get_srdb()
  acdb.getIDs(srdb, ids)
}



#' Match attributes
#'
#' Find positions in database which have attributes specified in '...'.  This function is same as acdb.find, and is curried to form srdb.find.
#'
#' Attribute inheritance behavior can be set with `inherit` and `inherit.gene` (see `?ag.attribute`)
#'
#' @param db list: an agdb (see `?agdb`)
#' @param ... attributes to match
#' @param inherit: bool: whether to inherit attributes from parents
#' @param inherit.gene char: which gene to inherit along
#' @param simplify bool: whether to simplify list of length 1 to `list[[1]]`
#'
#' @return list
#' @export
#'
#' @examples
acdb.find <- function(
  db,
  ...,
  inherit = TRUE,
  inherit.gene = "HA",
  simplify = TRUE
){

  # Get values and criteria
  values    <- list(...)
  criterias <- names(values)

  # Extract criteria
  db_criteria <- lapply(criterias, function(criteria){
    collate(sapply(db, function(ag){
      ag.attribute(ag, criteria, inherit = inherit, inherit.gene = inherit.gene)
    }))
  })

  # Compare values against criteria
  result <- lapply(seq_along(values[[1]]), function(i){
    matches <- do.call(cbind, lapply(seq_along(criterias), function(j){
      db_criteria[[j]] == values[[j]][i]
    }))
    which(rowSums(!matches) == 0)
  })

  # Return result, simplified if length one
  if(simplify && length(result) == 1){
    result[[1]]
  } else {
    result
  }

}






#' Match and extract attributes
#'
#' Find entries which have attributes specified in '...'. This function is same as agdb.search, and is curried to form srdb.search.
#'
#' Attribute inheritance behavior can be set with `inherit` and `inherit.gene` (see `?ag.attribute`)
#'
#' @param db list: an agdb (see `?agdb`)
#' @param ... attributes to match
#' @param inherit: bool: whether to inherit attributes from parents
#' @param inherit.gene char: which gene to inherit along
#' @param simplify bool: whether to simplify list of length 1 to `list[[1]]`
#'
#' @return list
#' @export
#'
#' @examples
acdb.search <- function(
  db,
  ...,
  inherit = TRUE,
  inherit.gene = "HA",
  simplify = TRUE
){

  if(missing(db)) db <- get_agdb()

  matches <- acdb.find(
    db, ...,
    inherit = inherit,
    inherit.gene = inherit.gene,
    simplify = FALSE
  )


  if(simplify && length(matches) == 1){
    db[matches[[1]]]
  } else {
    lapply(matches, function(x){
      db[x]
    })
  }

}


#' Match and extract a single entry from attributes
#'
#' Find an entry which has attributes specified in '...'. An error is thrown if there are multiple matches. Attribute inheritance behavior can be set with `inherit` and `inherit.gene` (see `?ag.attribute`).
#'
#' This function is the same as agdb.extract, and is curried to form srdb.extract.
#'
#' @param db list: an agdb (see `?agdb`)
#' @param ... attributes to match
#' @param inherit: bool: whether to inherit attributes from parents
#' @param inherit.gene char: which gene to inherit along
#' @param simplify bool: whether to simplify list of length 1 to `list[[1]]`
#'
#' @return list
#' @export
#'
#' @examples
acdb.extract <- function(
  db,
  ...,
  inherit = TRUE,
  inherit.gene = "HA"
){

  result <- acdb.search(
    db = db,
    ...,
    inherit = inherit,
    inherit.gene = inherit.gene
  )

  if(length(result) != 1){
    stop(sprintf("Search returned %s matches", length(result)))
  }

  result[[1]]

}





# Function to find name matches in an antigen or serum database
#' @export
acdb.matchNamesIndices <- function(strain_names, acdb, include_aliases = TRUE, multiple_match_error = TRUE){

  # Do some common replacements
  strain_names <- standardise_db_name(strain_names)

  # Get names from database
  dbnames <- standardise_db_name(collate(acdb%$%long))
  if(include_aliases) dbaliases <- lapply(acdb, function(x){
    standardise_db_name(x[["aliases"]])
  })

  # Check there are no duplicates
  db_strains_matches <- dbnames[dbnames %in% strain_names]
  if(sum(duplicated(db_strains_matches)) > 0){
    if(multiple_match_error){
      stop(sprintf(
        "\n\nMultiple matches found for:\n\n%s\n\n",
        paste(db_strains_matches[duplicated(db_strains_matches)], collapse = "\n")
      ))
    } else {
      warning(sprintf(
        "\n\nMultiple matches found for:\n\n%s\n\n",
        paste(db_strains_matches[duplicated(db_strains_matches)], collapse = "\n")
      ))
    }
  }

  # Find matches
  db_matches <- match(strain_names, dbnames)

  # Look for alias matches in any failures
  if(include_aliases){
    for(x in which(is.na(db_matches))){
      strain_match <- which(sapply(dbaliases, function(aliases){
        strain_names[x] %in% aliases
      }))
      if(length(strain_match) == 0) next
      if(length(strain_match) > 1) {
        if(multiple_match_error){
          stop("Multiple matches found for ", strain_names[x])
        } else {
          warning("Multiple matches found for ", strain_names[x])
        }
      }
      db_matches[x] <- strain_match
    }
  }

  # Return dbmatches
  attr(db_matches, "parsed_names") <- strain_names
  db_matches

}

#' Find name matches in an antigen or serum database
#'
#' @param strain_names names to be matched
#' @param acdb
#' @param include_aliases should aliases be searched?
#' @param multiple_match_error should a name matchingmultiple db entries throw an error?
#'
#' @return Returns a list of standardised names and extracted information
#' @export
#'
acdb.matchNames <- function(strain_names, acdb, include_aliases = TRUE, multiple_match_error = TRUE){

  # Get matching indices
  db_matches <- acdb.matchNamesIndices(strain_names, acdb, include_aliases, multiple_match_error)

  # Stop if some names not matched
  if(sum(is.na(db_matches)) > 0){
    stop("No matches found for:\n'", paste(strain_names[is.na(db_matches)], collapse = "',\n'"), "'")
    # stop("No matches found for:\n'", paste(attr(db_matches, "parsed_names")[is.na(db_matches)], collapse = "',\n'"), "'")
  }

  # Return the matching entries
  acdb[db_matches]

}




#' @export
acdb.matchSubstitutionsIndices <- function(ags, subs, gene = "HA", no_extras = TRUE){

  ag_subs <- agdb.substitutions(ags)
  which(vapply(ag_subs, function(ag_subs){
    if(no_extras){
      sum(!ag_subs %in% subs) == 0 && sum(!subs %in% ag_subs) == 0
    } else {
      sum(!subs %in% ag_subs) == 0
    }
  }, logical(1)))

}

#' @export
acdb.matchSubstitutions <- function(ags, subs, gene = "HA", no_extras = TRUE){

  ags[
    acdb.matchSubstitutionsIndices(
      ags = ags,
      subs = subs,
      gene = gene,
      no_extras = no_extras
    )
  ]

}

#' @export
srdb.matchAnimalIDs <- function(srdb, animal_id){
  matches <- acdb.search(srdb, animal_id = animal_id)
  match_lengths <- vapply(matches, length, numeric(1))
  if(sum(match_lengths > 1) > 0) stop(sprintf("\n\nMultiple matches found for:\n\n'%s'\n\n", paste(animal_id[match_lengths > 1])))
  if(sum(match_lengths < 1) > 0) stop(sprintf("\n\nNo matches found for:\n\n'%s'\n\n", paste(animal_id[match_lengths < 1])))
  unlist(matches, recursive = FALSE)
}

#' @export
srdb.matchAnimalIDsIndices <- function(srdb, animal_id){
  matches <- acdb.find(srdb, animal_id = animal_id)
  match_lengths <- vapply(matches, length, numeric(1))
  if(sum(match_lengths > 1) > 0) stop(sprintf("\n\nMultiple matches found for:\n\n'%s'\n\n", paste(animal_id[match_lengths > 1])))
  vapply(matches, function(x){
    if(length(x) == 0) return(NA)
    x
  }, numeric(1))
}

##########################
#
#    OPERATIONS
#
##########################

#' Apply a function over an acdb
#'
#' Efficiently applies a function over acdb.
#'
#' @param acdb list
#' @param fn function
#'
#' @return list
#' @export
#'
#' @examples
acdb.applyFunction <- function(acdb, fn){

  # Fetch ids and convert to factors
  ids        <- collate(lapply(acdb, function(x){ x$id }))
  ids_factor <- as.factor(ids)

  # Apply the function over the factor levels
  results <- lapply(
    acdb[match(levels(ids_factor), ids)],
    fn
  )

  # Return the matched results
  results[as.numeric(ids_factor)]

}







