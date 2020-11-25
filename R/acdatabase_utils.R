

##########################
#
#    SEARCHING
#
##########################

#' Match attributes
#'
#' Find positions in database which have attributes specified in '...'.  This function is same as agdb.find, and is curried to form srdb.find.
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

  matches <- agdb.find(
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
gene_inheritance <- function(
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
    gene_inheritance(
      ag            = ag[["alterations"]][[ha_alteration]][[".parent"]],
      gene          = gene,
      inheritance   = inheritance
    )
  } else if(!is.null(ag$parent_id)){
    gene_inheritance(
      ag            = ag[['.parent']],
      gene          = gene,
      inheritance   = inheritance
    )
  } else{
    inheritance
  }

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


