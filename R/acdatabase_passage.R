
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
  return(c(ag.passageChildrenFlat(ag.passageParent(ag), agdb), ag.passageParent(ag)))
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

#'@export
merge_on_passage <- function(titertable, agdb){

  # Fetch records
  ag_records  <- acdb.getIDs(agdb, rownames(titertable))

  # Get passage parent IDs but indicate if it included egg passage
  parent_ids  <- sapply(ag_records, function(ag){ passage_inheritance(ag)[[1]]$id })
  passage_egg <- sapply(ag_records, passage_includes_egg)
  parent_ids[passage_egg] <- paste(parent_ids[passage_egg], "EGG")

  # Determine duplicate rows
  duplicate_rows <- which(duplicated(parent_ids))

  # Merge together the passage duplicate rows
  for(x in duplicate_rows){

    matched_x <- match(parent_ids[x], parent_ids)
    titers1 <- titertable[matched_x,]
    titers2 <- titertable[x,]
    if(sum(titers1 != "*" & titers2 != "*") > 0){
      stop("Would have to merge titers")
    }
    na_titers1 <- titers1 == "*"
    titertable[matched_x,na_titers1] <- titers2[na_titers1]

  }

  # Remove duplicate rows
  titertable[-duplicate_rows,,drop=F]

}
