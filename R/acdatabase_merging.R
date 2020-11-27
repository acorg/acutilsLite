
#' Merge titers
#'
#' Merge titers (formatted like "240/480/240") by taking mean of logs
#'
#' @param titers char: like "240/480/240"
#'
#' @return char: merged titer
#' @export
#'
#' @examples
mergeTiters <- function(titers){
  titers <- unlist(lapply(titers, strsplit, split = "/", fixed = T))
  titers[titers == "*"] <- NA
  titers       <- titers[!is.na(titers)]
  if(length(titers) == 0) return("*")
  if(length(titers) == 1) return(titers)
  logtiters    <- Racmacs::titer_to_logtiter(titers)
  meanlogtiter <- mean(logtiters, na.rm = TRUE)
  meantiter    <- round(2^meanlogtiter*10)
  if(meantiter < 10) meantiter <- "<10"
  as.character(meantiter)
}

#' cCombine titers
#'
#' Combine vector of titers to "240/480/240"-like
#'
#' @param titers char or int
#'
#' @return char: combined titers
#' @export
#'
#' @examples
combineTiters <- function(titers){
  if(length(titers) == 0) return("*")
  if(length(titers) == 1) return(titers)
  titers <- titers[!is.na(titers) & titers != "*"]
  paste(titers, collapse = "/")
}

#' Merge from expdb
#'
#' Get merged titer table from experimental database.
#'
#' @param expdb char: like "240/480/240"
#' @param duplicate_titers char: "mean" or "combine"
#'
#' @return char: merged titer table
#' @export
#'
#' @examples
expdb.merge <- function(
  expdb,
  duplicate_titers = "mean"
){

  all_agids <- unique(unlist(lapply(expdb, function(exp){
    lapply(exp$results, function(result) {
      result$antigen_ids
    })
  })))

  all_srids <- unique(unlist(lapply(expdb, function(exp){
    lapply(exp$results, function(result) {
      result$serum_ids
    })
  })))

  all_results <- unlist(lapply(expdb, function(exp){
    exp$results
  }), recursive = FALSE)

  table_array <- array(dim = c(length(all_agids), length(all_srids),
                               length(all_results)))

  for (i in seq_along(all_results)) {
    result <- all_results[[i]]
    ag_matches <- match(result$antigen_ids, all_agids)
    sr_matches <- match(result$serum_ids, all_srids)
    table_array[ag_matches, sr_matches, i] <- result$titers
  }

  if(duplicate_titers == "mean")         { mergeFn <- mergeTiters   }
  else if(duplicate_titers == "combine") { mergeFn <- combineTiters }
  else { stop("'duplicate_titers' must be one of 'mean' or 'combine'") }

  table_merge <- apply(table_array, c(1, 2), mergeFn)
  rownames(table_merge) <- all_agids
  colnames(table_merge) <- all_srids
  table_merge
}

#' Merge from exp
#'
#' Get merged titer table from single experiment.
#'
#' @param exp char: like "240/480/240"
#' @param expect_repeats bool: whether repeat titers are expected
#' @param merge_passage bool: not implemented
#'
#' @return char: merged titer table
#'
#' @rdname
#' @export
#'
#' @examples
exper.merge <- function(
  exp,
  expect_repeats = FALSE,
  merge_passage = FALSE
){

  # Check input
  if(!inherits(exp, "acdatabase.exp")){
    stop("Input must be of class acdatabase.exp")
  }

  # Get the unique ids and setup the titer array
  all_agids <- unique(unlist(lapply(exp$results, function(result){ result$antigen_ids })))
  all_srids <- unique(unlist(lapply(exp$results, function(result){ result$serum_ids })))
  table_array <- array(dim = c(length(all_agids), length(all_srids), length(exp$results)))
  hau_table   <- matrix(ncol = length(all_agids), nrow = length(exp$results))

  # Populate the titer array
  for(i in seq_along(exp$results)){
    result <- exp$results[[i]]
    ag_matches <- match(result$antigen_ids, all_agids)
    sr_matches <- match(result$serum_ids, all_srids)
    table_array[ag_matches, sr_matches, i] <- result$titers
    if(!is.null(result$meta$hau)){
      hau_table[i, ag_matches] <- result$meta$hau
    }
  }

  # Check for repeats if not expected
  if(!expect_repeats){
    if(sum(apply(!is.na(table_array), c(1,2), sum) > 1) > 0){
      stop("There are repeat measurements, that were not specified as expected")
    }
  }

  # Merge HAU data
  hau <- apply(hau_table, 2, function(x){
    paste0(unique(x[!is.na(x)]), collapse = ";")
  })
  names(hau) <- all_agids

  # Perform the merge
  table_merge <- apply(table_array, c(1,2), mergeTiters)
  rownames(table_merge) <- all_agids
  colnames(table_merge) <- all_srids

  # Add HAU data as an attribute
  if(sum(hau != "") > 0){
    attr(table_merge, "hau") <- hau
  }

  # Return the result
  table_merge

}


#' @export
titertable.toLong <- function(
  titertable,
  agdb,
  srdb
){

  tibble::as_tibble(
    titertable,
    rownames = "ag"
  ) %>%
    tidyr::pivot_longer(
      cols      = -ag,
      names_to  = "sr",
      values_to = "titer"
    ) %>%
    dplyr::mutate(
      ag_records = acdb.getIDs(ag, agdb),
      sr_records = acdb.getIDs(sr, srdb),
      srag_records = srdb.homologousAntigens(sr_records, agdb)
    )

}


