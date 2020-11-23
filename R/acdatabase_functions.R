
# Get the merged results
#' @export
resultsTable <- function(exp,
                         agids = NULL,
                         srids = NULL){

  # Get antigen ids
  allagids <- unique(unlist(lapply(exp$results, function(results){
    results$antigen_ids
  })))

  # Get serum ids
  allsrids <- unique(unlist(lapply(exp$results, function(results){
    results$serum_ids
  })))

  # Assemble a table
  titers <- matrix("*", length(allagids), length(allsrids))
  rownames(titers) <- allagids
  colnames(titers) <- allsrids

  # Fill in the table
  for(result in exp$results){
    titers[match(result$antigen_ids, allagids),
           match(result$serum_ids, allsrids)] <- result$titers
  }

  # Set antigen and sera defaults
  if(is.null(agids)) agids <- allagids
  if(is.null(srids)) srids <- allsrids

  # Return the titer table
  titers[agids, srids]

}




# Convert to a long table format
#' @export
gatherResultsTable <- function(exp,
                               agdb = NULL,
                               srdb = NULL,
                               ag_attributes = NULL,
                               sr_attributes = NULL,
                               inherit.lists = TRUE){

  # Read the results table
  results_table <- resultsTable(exp)

  # Convert to the long format
  results_long <- gather(
    data  = as_tibble(results_table, rownames = "ag.id"),
    key   = "sr.id",
    value = "titer",
    -ag.id
  )

  # Get the antigen attributes table
  if(!is.null(ag_attributes)){

    # Get the attributes
    ag_attributes <- attributesTable(
      db              = agdb,
      attributes      = ag_attributes,
      ids             = results_long$ag.id,
      colnames.prefix = "ag.",
      inherit.lists   = inherit.lists
    )

    # Bind to the long table
    results_long <- dplyr::bind_cols(
      ag_attributes,
      results_long
    )

  }

  # Get the serum attributes table
  if(!is.null(sr_attributes)){

    # Get the attributes
    sr_attributes <- attributesTable(
      db              = srdb,
      attributes      = sr_attributes,
      ids             = results_long$sr.id,
      colnames.prefix = "sr.",
      inherit.lists   = inherit.lists
    )

    # Bind to the long table
    results_long <- dplyr::bind_cols(
      sr_attributes,
      results_long
    )

  }

  # Return the long table
  results_long

}



# Fetch attributes for an antigen
#' @export
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
    ag_inheritance <- gene_inheritance(
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

#' @export
sr.attribute <- ag.attribute


# Assemble a table of antigen attributes
#' @export
attributesTable <- function(
  db,
  attributes,
  ids = NULL,
  inherit.lists = TRUE,
  colnames.prefix = NULL
){

  if(is.null(ids)){

    # If ids are null just get all the attributes from every antigen
    values <- .attributesTable(
      db              = db,
      attributes      = attributes,
      ids             = ids,
      inherit.lists   = inherit.lists,
      colnames.prefix = colnames.prefix
    )

  } else {

    # Otherwise try and be smart and get only attributes for unique ids
    unique_ids <- unique(ids)
    unique_values <- .attributesTable(
      db              = db,
      attributes      = attributes,
      ids             = unique_ids,
      inherit.lists   = inherit.lists,
      colnames.prefix = colnames.prefix
    )

    # Match up the values
    values <- unique_values[match(ids, unique_ids),]

  }

  # Return the values
  values

}

# base function for the antigen attributes table
.attributesTable <- function(
  db,
  attributes,
  ids = NULL,
  inherit.lists = TRUE,
  colnames.prefix = NULL
){


  # Get database subset
  if(!is.null(ids)){

    dbids <- db%$%id
    if(sum(!ids %in% dbids) > 0) stop("ID not found in database")
    db <- db[match(ids, dbids)]

  }

  # Get attributes
  values <- lapply(attributes, function(attribute){

    vals <- lapply(db, ag.attribute, attribute = attribute, inherit.lists = TRUE)
    vals[sapply(vals, is.null)] <- NA
    unlist(vals)

  })

  # Gather the values as a table
  names(values) <- attributes
  values_tibble <- as_tibble(values)

  # Prefix column names
  if(!is.null(colnames.prefix)){
    colnames(values_tibble) <- paste0(colnames.prefix, colnames(values_tibble))
  }

  # Return the tibble
  values_tibble

}

