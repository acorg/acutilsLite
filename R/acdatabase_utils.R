##########################
#
#    HELPER
#
##########################

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


#####
# general attributes
#####



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



#####
# names
#####





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
  return(ag[['long']])
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
ag.short <- function(ag){

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
  if(ag.passageIncludesEgg(ag)) year_short <- paste0(year_short, "E")

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

}



#####
# parents and inheritance
#####


#' Get parent
#'
#' @param ag environment: an antigen
#'
#' @return list
#' @export
#'
#' @examples
ag.parent <- function(ag, gene = 'BB'){

  ha_alteration = match(gene, sapply(ag[["alterations"]], function(x) x$gene))

  # search for gene parent in alterations
  if (!is.na(ha_alteration) & !is.null(ag$alterations[[ha_alteration]]$parent_id)) {
    return( ag[['alterations']][[ha_alteration]]$.parent )
  }
  # check for ag$.parent
  else if (!is.null(ag$parent_id)) return( ag$.parent )
  else return(NULL)
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
#' Passed antigen is last element of list
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

  # Add antigen to inheritance
  inheritance <- c(list(ag), inheritance)
  parent = ag.parent(ag, gene)
  if (length(parent) == 0) return(inheritance)
  else return( ag.inheritance(parent, gene, inheritance) )

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
ag.children <- function(ag, agdb, gene = 'BB', parents = NULL){
  if (is.null(parents)) parents = acdb.applyFunction(agdb, ag.parent, gene)
  parent_ids = parents%$%id

  parent_ids = (lapply(parent_ids, function(x) {
    if (is.null(x)) return (NA)
    else return(x)} ))


  agdb[which(parent_ids == ag$id)]
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
ag.descendentsTree <- function(ag, agdb, gene = 'BB', .parents = NULL){
  if (is.null(.parents)) .parents <- acdb.applyFunction(agdb, ag.parent, gene)
  ag_children <- ag.children(ag, agdb, gene, .parents)
  list(ag = ag, children = acdb.applyFunction(ag_children, ag.descendentsTree, agdb = agdb, gene = gene, .parents = .parents ))
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
ag.descendents <- function(ag, agdb, gene = 'BB'){
  ag_tree <- ag.descendentsTree(ag, agdb, gene)
  descendents <- list()
  tree_children <- function(ag_tree){
    descendents <<- c(descendents, ag_tree$ag)
    lapply(ag_tree$children, tree_children)
  }
  tree_children(ag_tree)
  descendents[-1]
}


#####
# sequences
#####


#' Get sequence for a gene from a record
ag.sequence.pluck <- function(ag, gene = 'HA'){

  # Get the ha alterations
  ha_sequence <- match(gene, sapply(ag[["genes"]], function(x) x$gene))

  # Add any substitutions listed
  if(!is.na(ha_sequence)){
    sequence <- ag[["genes"]][[ha_sequence]][["sequence"]]
  } else {
    sequence <- NA
  }

  # Return the substitutions
  sequence

}

#'@export
ag.sequence.pluckFromID <- function(id, agdb = get_agdb(), gene = 'HA'){
  acutilsLite:::ag.sequence.pluck(acdb.extract(agdb, id = id), gene)
}

ag.sequence.applySubs <- function(sequence, subs){

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

  sequence

}

#' Calculate sequence for a gene from an inheritance (by finding a sequence and applying substitutions)
ag.sequence.calculate <- function(inheritance, gene = 'HA'){
  # Go through the inheritance tree
  subs <- list()
  for(ag in rev(inheritance)){
    sequence <- ag.sequence.pluck(ag, gene)
    if(!is.na(sequence)) break
    subs     <- c(subs, list(ag.substitutions(ag, gene)))

  }
  id = ag$id
  if(is.na(sequence)) return(NA)

  # Apply any modifications to the sequence
  sequence = acutilsLite:::ag.sequence.applySubs(sequence, subs)

  # Return the sequence
  seq = sequence
  names(seq) = paste0('',id, '')
  return(seq)
}

format_passage_history <- function(passage_history){

  passage_history[sapply(passage_history, is.null)] = '.'

  passage_history=sapply(passage_history, combineTiters)

  passage_history=paste0(passage_history, collapse = ' | ')

  return(passage_history)
}

#' Calculate gene sequence
#'
#' Returns tibble of all sequences found for an antigen, and their relationship to the passed antigen
#'
#' @param ag environment: an antigen
#' @param gene char: which gene to get sequence for
#'
#' @return list
#' @export
#' @examples
ag.sequence.allVariants <- function(
  ag,
  gene = 'HA',
  agdb = acutils:::get_agdb(),
  inherit = T,
  full_search = F
){
  sequences_found = tibble(id = character(), ag_record = list(), passage_history = character(), relationship = character(), sequence = character())
  # is the ag sequenced?
  sequence_self = ag.sequence.pluck(ag, gene)
  sequences_found = sequences_found = dplyr::add_row(sequences_found,
                                                     id=ag$id,
                                                     ag_record = list(ag),
                                                     passage_history = format_passage_history(list(unlist(ag.inheritedPassageHistory(ag)))),
                                                     relationship='self',
                                                     sequence=sequence_self)

  # are passage variants sequenced?
  passage_sibs = ag.passageSiblings(ag, agdb)
  passage_sibs = passage_sibs[passage_sibs%$%id != ag$id]
  sequence_passage_sibs = sapply(passage_sibs, ag.sequence.pluck, gene)
  names(sequence_passage_sibs) = passage_sibs%$%id
  sequence_passage_sibs = sequence_passage_sibs[!is.na(sequence_passage_sibs)]

  for (i in seq_along(sequence_passage_sibs)){
    sequences_found = dplyr::add_row(sequences_found,
                                     id = names(sequence_passage_sibs)[[i]],
                                     ag_record = list(passage_sibs[[i]]),
                                     passage_history = format_passage_history(list(unlist(ag.inheritedPassageHistory(passage_sibs[[i]])))),
                                     relationship = 'passage variant',
                                     sequence = sequence_passage_sibs[[i]])
  }



  # are parents sequenced?
  inheritance = ag.inheritance(ag, gene)
  sequence_parents = ag.sequence.calculate(inheritance, gene)

  if (!is.na(sequence_parents) & isFALSE(names(sequence_parents) %in% sequences_found$id)){
    id_parent = names(sequence_parents)
    sequences_found = dplyr::add_row(sequences_found,
                                     id = id_parent,
                                     ag_record = unlist(acdb.getIDs(id_parent, agdb)),
                                     passage_history = format_passage_history(list(unlist(ag.inheritedPassageHistory(acdb.getIDs(id_parent, agdb)[[1]])))),
                                     relationship='parent',
                                     sequence=sequence_parents)
  }


  # are parent passage variants sequenced?
  sequences_passaged_parents = list()
  for (i in seq_along(inheritance[-length(inheritance)])){
    parent_i = inheritance[[i]]

    # list of parent roots
    sequences_passaged_parents[[parent_i$id]] = character()
    parent_i_sibs = ag.passageSiblings(parent_i, agdb)



    # try calculating sequence with parent_i replaced by one of its siblings in the inheritance
    for (sib in parent_i_sibs){
      sequence_original = acutilsLite:::ag.sequence.pluck(sib)
      if (!is.na(sequence_original) & !(sib$id %in% sequences_found$id)){
        subs = unlist(agdb.substitutions(inheritance[(i+1):length(inheritance)], inherit = F))
        if (!is.na(subs)) sequence_original = ag.sequence.applySubs(sequence_original, subs)

        sequences_found = dplyr::add_row(sequences_found,
                                         id=sib$id,
                                         ag_record = list(unlist(sib)),
                                         passage_history = format_passage_history(list(unlist(ag.inheritedPassageHistory(unlist(sib))))),
                                         relationship='parent passage variant',
                                         sequence=sequence_original)
      }

    }

  }


  return(sequences_found)
}


#'@export
ag.inheritance.variants <- function(ag, agdb = get_agdb()){
  inheritance = ag.inheritance(ag, 'HA')
  inheritance_sibs = lapply(inheritance, ag.passageSiblings, agdb)
  inheritance_sibs = lapply(inheritance_sibs, function(sibs){sibs = sibs[sibs%$%id != ag$id]})
  inheritance_sibs = lapply(inheritance_sibs, function(sibs){names(sibs) = sibs%$%id})

  for (i in seq_along(inheritance_sibs)){
    sibs = inheritance_sibs[[i]]
    for(j in seq_along(sibs)){
      sib = sibs[[j]]
      if (!is.na(acutilsLite:::ag.sequence.pluck(sib))){
        names(inheritance_sibs)[[i]][[j]] = paste0('*', names(inheritance_sibs[[i]][[j]]), '*')
      }
    }
  }

  return(names(inheritance_sibs))
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


#'@export
ag.sequence.auto <- function(){}


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




#####
# clades and clusters
#####


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

  if (stringr::str_length(identifier) == 6 & from != 'id') print('ARE YOU SURE THIS IS NOT AN id?')
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

  homologous_ag = acdb.getIDs(sr$strain_id, agdb)
  return(ag.cluster(homologous_ag))

}




#' Get antigen clade
#'
#'
#' @param ag
#' @param gene
#'
#' @return char
#' @export
#'
#' @examples
ag.clade.self <- function(ag){
  clades = intersect(ag$groups, h3_clade_grouping_heirarchy)
  if (length(clades) == 0) return(NA)
  first_matching_clade = h3_clade_grouping_heirarchy[[min(which(h3_clade_grouping_heirarchy %in% clades))]]
  return(first_matching_clade)
}

#' Get antigen clade
#' You can choose whether clade is searched for only in ag (how = 'self), only in the root virus (how = 'root), or anywhere in the inheritance (how = 'any')
#'
#' @param ag
#' @param how = c('self', 'root', 'any')
#'
#' @return char
#' @export
#'
#' @examples
ag.clade <- function(ag, how = 'self', agdb = acutilsLite:::get_agdb()){
  inh = ag.inheritance(ag, gene = 'HA')
  inh.clades = lapply(inh, ag.clade.self)
  if (how == 'root') return(inh.clades[[1]])
  if (how == 'self') return(rev(inh.clades)[[1]])
  if (how == 'any') {
    # check self:
    if (!is.na(ag.clade.self(ag))) return(ag.clade.self(ag))

    # check inheritance
    if (length(unique(  inh.clades[!is.na(inh.clades)]  )) > 1) stop('Multiple clade matches:', unique(  inh.clades[!is.na(inh.clades)]  ))
    else if (length(unique(  inh.clades[!is.na(inh.clades)]  )) == 1) return(unlist(unique(  inh.clades[!is.na(inh.clades)]  )))



    # check descendants
    desc = ag.descendents(ag, agdb = agdb, gene = 'HA')
    desc.clades = lapply(inh, ag.clade.self)

    if (length(unique(  desc.clades[!is.na(desc.clades)]  )) > 1) stop('Multiple clade matches:', unique(  desc.clades[!is.na(desc.clades)]  ))
    else if (length(unique(  desc.clades[!is.na(desc.clades)]  )) == 1) return(unlist(unique(  desc.clades[!is.na(desc.clades)]  )))

    return('Clade missing')
  }
  return(ag.clade(ag.root(ag)))
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
ac.inGroups <- function(ac, ..., all = T){
  groups <- c(...)

  if (all) return(sum(!groups %in% ac[["groups"]]) == 0)

  return(sum(groups %in% ac[["groups"]]) > 0 )
}


#####
# passaging
#####

#' Get passage root
#'@export
ag.passageParent <- function(ag){
  ag.passageInheritance(ag)[[1]]
}


#' Get passage inheritance
#'@export
ag.passageInheritance <- function(ag){
  inheritance <- list()
  while(1==1){
    inheritance <- c(inheritance, list(ag))
    if(is.null(ag$parent_id) || !is.null(ag[["alterations"]])) break
    ag <- ag$.parent
  }
  rev(inheritance)

}


#' Get passage children
#'@export
ag.passageChildren <- function(ag, agdb){
  children = ag.children(ag, agdb)
  children = Filter(function(ag){is.null(ag[['alterations']])}, children)
  return(children)
}


#' Get all antigens which differ only by passage
#'@export
ag.passageSiblings = function(ag, agdb){
  passage_sibs = c(ag.passageChildrenFlat(ag.passageParent(ag), agdb), ag.passageParent(ag))
  passage_sibs = passage_sibs[passage_sibs%$%id != ag$id]
  return(passage_sibs)
}

#' Get all passage descendants
#' @export
ag.passageChildrenTree <- function(ag, agdb){
  ag_children <- ag.passageChildren(ag, agdb)
  list(ag = ag, children = lapply(ag_children, ag.passageChildrenTree, agdb = agdb))
}


#' Get all passage descendants
#' @export
ag.passageChildrenFlat <- function(ag, agdb){
  ag_tree <- ag.passageChildrenTree(ag, agdb)
  descendents <- list()
  tree_children <- function(ag_tree){
    descendents <<- c(descendents, ag_tree$ag)
    lapply(ag_tree$children, tree_children)
  }
  tree_children(ag_tree)
  descendents[-1]
}



#'@export
ag.passageIncludesEgg <- function(ag){

  passage_history <- unlist(ag.inheritedPassageHistory(ag))
  sum(grepl("E", passage_history, fixed = T)) > 0

}

#'@export
ag.inheritedPassageHistory <- function(ag){

  ag_passage_inheritance <- ag.passageInheritance(ag)
  lapply(ag_passage_inheritance, function(x){ x[["passage"]][["history"]] })

}

#####
# experiments
#####


#' Get experiments that an antigen is in
#' @export
ag.experiments <- function(
  ag,
  expdb,
  include_passage_variants = FALSE,
  agdb
){
  ag_id <- ag$id


  # Keep a list of matching experiments
  ag_exps <- list()

  # Work out the experiment ids to search on
  for(exp in expdb){
    exp_ag_ids <- expdb.ags(list(exp))
    if(include_passage_variants){
      exp_ags <- acdb.getIDs(exp_ag_ids, agdb)
      exp_ag_ids <- vapply(
        exp_ags,
        function(ag){
          ag.passageSiblings(ag, agdb)[[1]]$id
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

#' Get experiments that a serum is in
#' @export
sr.experiments <- function(
  sr,
  expdb,
  include_passage_variants = FALSE, #not implemented
  srdb
){
  sr_id <- sr$id


  # Keep a list of matching experiments
  sr_exps <- list()

  # Work out the experiment ids to search on
  for(exp in expdb){
    exp_sr_ids <- expdb.sr(list(exp))
    if(include_passage_variants){
      # pass
    }

    # Add experiment to list if antigen is found
    if(sr_id %in% exp_sr_ids){
      sr_exps <- c(sr_exps, exp)
    }
  }

  # Return the list of matching experiments
  sr_exps

}


# Check if antigens are present in certain experimental results
#' @export
ag.inResults <- function(agdb, exp){
  agids <- collate(agdb%$%id)
  lapply(agids, function(id){
    which(vapply(exp$results, function(result){
      id %in% result$antigen_ids
    }, logical(1)))
  })
}

# Check if sera are present in certain experimental results
#' @export
sr.inResults <- function(srdb, exp){
  srids <- collate(srdb%$%id)
  lapply(srids, function(id){
    which(vapply(exp$results, function(result){
      id %in% result$serum_ids
    }, logical(1)))
  })
}



#######
# From acdb
#######


#####
# general attributes
#####


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


#' @export
agdb.year <- function(agdb, gene = 'HA'){
  acdb.applyFunction(agdb, function(x){
    as.numeric(ag.attribute(x, "isolation.date", inherit.gene = gene),
               unlist=T)
  })
}


#####
# names
######


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
  if (name == 'short') return(agdb.short(acdb.getIDs(ids, acdb)))
  else stop('name must be "long" or "short"')
  }


#' Generate short name from antigen entries
#'
#' @param ags
#'
#' @return Returns a vector of short names
#' @export
#'
agdb.short <- function(agdb){
  acdb.applyFunction(agdb, ag.short, unlist = T)}



#' Select long name from antigen entries
#'
#' @param ags
#'
#' @return Returns a vector of long names
#' @export
#'
acdb.long <- function(acdb){
  acdb.applyFunction(acdb, function(x){ x[["long"]] }, unlist = T)
}



#####
# clades and clusters
#####



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


#' @export
agdb.root <- function(agdb, gene){
  acdb.applyFunction(agdb, function(ag){
    ag.root(ag, gene)
  })
}

#####
# ag-sr homology
#####


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



#####
# sequences
#####



#' Get substitutions from agdb
#'
#' @param agdb
#'
#' @return
#' @export
agdb.substitutions <- function(agdb, inherit = T){
  o = acdb.applyFunction(agdb, function(ag){ unlist(ag.substitutions(ag, gene = "HA", inherit)) })
  lapply(o, function(x){
    if (is.null(x)) return(NA)
    else return(x)
    })
}



#' Are there subsitutions in the inheritance of HA?
#'
#' @param agdb
#'
#' @return
#' @export
agdb.hasSubstitutions <- function(agdb, inherit = T){
  acdb.applyFunction(agdb, function(ag){ length(unlist(ag.substitutions(ag, gene = "HA", inherit))) > 0 }, unlist = T)
}


#' Get sequences from agdb
#'
#' @param agdb
#'
#' @return
#' @export
agdb.sequencesFromMap <- function(ags, map, gene = 'HA', inherit = T, agdb = get_agdb()){

  acdb.applyFunction(ags,
                     function(ag){
                       inh = ag.inheritance(acdb.extract(agdb, id = map[ag$id]), gene)
                       ag.sequence.calculate(inh, gene)
                       },
                     unlist = T)
}






#####
# ordering
#####



#'@export
agdb.getOrder <- function(agdb){
  agdb = unique(agdb)

  clades = acdb.applyFunction(agdb, ag.clade, 'any')
  years = agdb.year(agdb)
  w = lapply(h3_clade_order, function(clade){
    which(clades == clade)
  })

  out = unlist(sapply(w, function(x){
    x[order(years[x])]
  }))
  return(unlist(agdb%$%id)[out])
}

#'@export
srdb.getOrder <- function(sr, agdb){
  sr = unique(sr)
  srags = srdb.homologousAntigens(sr, agdb)
  srags_order = agdb.getOrder(srags)
  sr_order = unlist(lapply(srags_order, function(srag_id)return(which(unlist(sr%$%strain_id) == srag_id))))
  return(unlist(sr%$%id)[sr_order])
}

#'@export
agdb.getRank <- function(agdb){
  agdb = unique(agdb)
  order = agdb.getOrder(agdb)
  rank = sapply(agdb%$%id, function(n)which(order == n))
  names(rank) = agdb%$%id

  return(rank)
}

#'@export
srdb.getRank <- function(sr, agdb){
  sr = unique(sr)
  sr_order = srdb.getOrder(sr, agdb)
  sr_rank = sapply(sr%$%id, function(n)which(sr_order == n))
  names(sr_rank) = sr%$%id
  return(sr_rank)

}

#####
# experiments
#####

#' Get experiments that antigens are in
#' @export
agdb.experiments <- function(ags, expdb, include_passage_variants, agdb){
  acdb.applyFunction(ags, function(x){ ag.experiments(x, expdb, include_passage_variants, agdb) })
}

# Get experiment names from expdb
#' @export
expdb.expName <- function(expdb){
  acdb.applyFunction(expdb, function(x){ x[["name"]] }, unlist = T)
}


#' Get antigens in expdb
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

#' Get sera in expdb
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


##########################
#
#    misc
#
##########################



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






##########################
#
#    SEARCHING
#
##########################

#####
# by attribute
#####


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
acdb.getIDs <- function(ids, db){

  dbids <- unlist(db%$%id)
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
  acdb.getIDs(ids, agdb)
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
  acdb.getIDs(ids, srdb)
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

  matches = do.call(cbind, lapply(seq_along(db_criteria), function(i){db_criteria[[i]] == values[[i]]}))
  result = which(rowSums(!matches) == 0)


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
      db[[x]]
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



#####
# match by name
#####


#' Function to find db entries with names matching passed strain names
#' Useful for finding antigen match when adding new tables to database
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
#' @param multiple_match_error should a name matching multiple db entries throw an error?
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



#' Find db entries with same subsititutions as subs
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

#' Find db entries with same subsititutions as subs
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

#' Find sera matching an animal id
#' @export
srdb.matchAnimalIDs <- function(srdb, animal_id){
  matches <- acdb.search(srdb, animal_id = animal_id)
  match_lengths <- vapply(matches, length, numeric(1))
  if(sum(match_lengths > 1) > 0) stop(sprintf("\n\nMultiple matches found for:\n\n'%s'\n\n", paste(animal_id[match_lengths > 1])))
  if(sum(match_lengths < 1) > 0) stop(sprintf("\n\nNo matches found for:\n\n'%s'\n\n", paste(animal_id[match_lengths < 1])))
  unlist(matches, recursive = FALSE)
}

#' Find sera matching an animal id
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

#'@export
unlist_safe <- function(l){
    if (length(l) == 0) return(l)
    long_entries <- sum(sapply(l, function(y) length(y) > 1 )) != 0
    null_entries <- sum(sapply(l, function(y) length(y) == 0 )) != 0
    if (long_entries){
      message('some list entries with length > 1')
      return(l)
    }
    if (null_entries) {
      message('NULL coersed to NA')
      l[lengths(l) == 0] = NA
      return(unlist(l))
    }
    return(unlist(l))
}

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
acdb.applyFunction <- function(acdb, fn, ..., unlist = F){

  # Fetch ids and convert to factors
  ids        <- collate(lapply(acdb, function(x){ x$id }))
  ids_factor <- as.factor(ids)

  # Apply the function over the factor levels
  results <- lapply(
    acdb[match(levels(ids_factor), ids)],
    fn,
    ...
  )

  # Return the matched results
  matched_res = results[as.numeric(ids_factor)]

  if (unlist) return(unlist_safe(matched_res))

  return(matched_res)

}




#' Apply a function over an acdb subset, found by id
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
acids.applyFunction <- function(ids, acdb, fn){
  acdb <- acdb.getIDs(ids, acdb)
  return(acdb.applyFunction(acdb, fn))
}

##########################
#
#    OTHER
#
##########################

#'@export
toTable = function(list_of_vectors){
  longest = max(sapply(list_of_vectors, length))

  list_of_vectors = lapply(list_of_vectors, function(v){c(v, rep('', longest - length(v)))})

  mat = do.call(cbind, list_of_vectors)

  colnames(mat) = names(list_of_vectors)
  mat
}


#'@export
simpleTable = function(text_entries, ...){

  extra_properties = list(...)

  out_table = text_entries


  for (col_i in seq_len(dim(text_entries)[[2]]) ){


      arg_list = list(x = out_table[,col_i])

      for (i in seq_along(extra_properties)){
        property_name = names(extra_properties)[[i]]
        property_value = unname(extra_properties)[[i]]

        arg_list[property_name] = list(property_value[,col_i])
      }

      out_table[,col_i] = do.call(kableExtra::cell_spec, arg_list)
  }

  out_table_html = (kableExtra::kbl(out_table, format = 'html',  escape = F)
                                )


  kableExtra::kable_material_dark(out_table_html)
}


