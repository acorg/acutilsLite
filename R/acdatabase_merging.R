##########################
#
#    converting between types
#
##########################

#' titer to logtiter
#'@export
titer.toLog = function(titers){
  titer_types <- Racmacs::titer_types(titers)
  threshold_titers <- titer_types == "lessthan" | titer_types ==
    "morethan"
  log_titers <- titers
  log_titers[titer_types == "omitted"] <- NA
  log_titers[threshold_titers] <- substr(log_titers[threshold_titers],
                                         2, nchar(log_titers[threshold_titers]))
  mode(log_titers) <- "numeric"
  log_titers <- log2(log_titers/10)
  mode(log_titers) <- "character"
  log_titers[titer_types == "lessthan"] <- paste0('<', log_titers[titer_types == "lessthan"])
  log_titers[titer_types == "morethan"] <- paste0('>' , log_titers[titer_types == "morethan"])
  #if ('morethan' %in% titer_types) #message('More than titers detected!')

  log_titers
}




#' this converts logtiters with less thans to logtiters where <3 -> 2
#'@export
logtiter.toPlot <- function (titers)
{

  titers.checkThresholds(titers)

  titer_types <- Racmacs::titer_types(titers)
  threshold_titers <- titer_types == "lessthan" | titer_types ==
    "morethan"
  log_titers <- titers
  log_titers[titer_types == "omitted"] <- NA
  log_titers[threshold_titers] <- substr(log_titers[threshold_titers],
                                         2, nchar(log_titers[threshold_titers]))
  mode(log_titers) <- "numeric"
  if (!all(titer_types == 'measured') ) {
    threshold_location = max(log_titers[threshold_titers])
    #message(paste0('<', as.character(threshold_location), ' converted to ', as.character(threshold_location-1)))
  }

  log_titers[titer_types == "lessthan"] <- log_titers[titer_types == "lessthan"] - 1
  log_titers[titer_types == "morethan"] <- log_titers[titer_types == "morethan"] + 1
  #if ('morethan' %in% titer_types) #message('More than titers detected!')
  return(log_titers)
}

#' ?
#'@export
titer.toPrint = function(titers){

  numeric_titers = sapply(titers, function(t){
    if ( stringr::str_sub(t, 1,1) %in% c('<', '>')) return(stringr::str_sub(t, 2))
    else return(t)
  })

  ns = unique(as.character(sort(as.numeric(numeric_titers))))
  ns=as.character(1:10000)
  lvls = c(paste0('<', ns), ns, paste0('>', ns))

  titers_fct = factor(as.character(titers), levels = lvls)

  return(titers_fct)

}

#' remove less thans or greater thans
#'@export
titers.toNumeric = function(titers){
  if(length(titers) == 0) return(c())
  titers = unlist(titers)
  return(sapply(titers, function(t){
    if ( stringr::str_sub(t, 1,1) %in% c('<', '>')) return(as.numeric(stringr::str_sub(t, 2)))
    else return(as.numeric(t))
  }))

}

#'@export
titers.checkThresholds <- function(titers){
  titers = unlist(titers)
  thresholds = titers.getThresholds(titers)

  if (length(thresholds$less) > 1 | length(thresholds$more) > 1) stop('Multiple conflicting thresholds present. Please change the threshold manually with titers.raiseThreshold, or consider splitting up your plot.')
  if (isTRUE(max(titers.toNumeric(titers)) != titers.toNumeric(thresholds$more)) | isTRUE(min(titers.toNumeric(titers)) != titers.toNumeric(thresholds$less))){
    stop('Some titers fall ouside the threshold range. Please change threshold manually with titers.raiseThreshold')
  }
}

#' find out what thresholded titer values are present
#'@export
titers.getThresholds <- function(titers){
  lesstiters = unique(titers[Racmacs::titer_types(titers) == 'lessthan'])
  moretiters = unique(titers[Racmacs::titer_types(titers) == 'morethan'])

  if (length(lesstiters) > 1) message('Multiple lessthan threshold values. This is most common when merging titers from multiple experiments')
  if (length(moretiters) > 1) message('Multiple lessthan threshold values. This is most common when merging titers from multiple experiments')

  list(less = titers.toNumeric(lesstiters), more = titers.toNumeric(moretiters))

}

#'@export
titers.changeThreshold = function(titers, newthreshold = 'infer', which = 'lower'){

  if (which == 'lower'){
    if (newthreshold == 'infer'){
      if (length(suppressMessages(titers.getThresholds(titers)$less)) == 0) return(titers)
      newthreshold = suppressMessages(max(titers.getThresholds(titers)$less))
    }
    numerictiters = titers.toNumeric(titers)
    titers[numerictiters < newthreshold] = paste0('<', newthreshold)
  }


  if (which == 'upper'){
    if (newthreshold == 'infer'){
      if (length(suppressMessages(titers.getThresholds(titers)$more)) == 0) return(titers)
      newthreshold = suppressMessages(min(titers.getThresholds(titers)$more))
    }
    numerictiters = titers.toNumeric(titers)
    titers[numerictiters > newthreshold] = paste0('>', newthreshold)
  }
  titers
}






##########################
#
#   merging tables
#
##########################


#####
# handling titers
#####



#' Merge titers
#'
#' Merge titers (formatted like "240/480/240") by taking mean of logs (and exponentiating back)
#'
#' When logs are taken <40 is converted to 20 for averaging
#'
#' @param titers char: like "240/480/240"
#'
#' @return char: merged titer
#' @export
#'
#' @examples
mergeTiters <- function(titers, lowerthreshold = 10, upperthreshold = Inf, threshold_masking = F){

  if (any(Racmacs::titer_types(titers) == 'lessthan'))
    if (max( titers.getThresholds(titers)$less ) > lowerthreshold ) lowerthreshold = max(titers.getThresholds(titers)$less)

    if (any(Racmacs::titer_types(titers) == 'morethan')) upperthreshold = min(titers.getThresholds(titers)$more)


  titers <- unlist(lapply(titers, strsplit, split = "/", fixed = T))
  titers[titers == "*"] <- NA
  titers       <- titers[!is.na(titers)]
  if(length(titers) == 0) return("*")
  if(length(titers) == 1 & !threshold_masking) return(titers)
  logtiters    <- Racmacs::titer_to_logtiter(titers)
  meanlogtiter <- mean(logtiters, na.rm = TRUE)
  meantiter    <- round(2^meanlogtiter*10)
  if(isTRUE(meantiter < lowerthreshold)){
    meantiter <- paste0('<', as.character(lowerthreshold))
  }


  if(isTRUE(meantiter > as.numeric(upperthreshold))){
    meantiter <- paste0('>', as.character(upperthreshold))
    #warning("titers ", combineTiters(titers), " merged to ", meantiter, ' (upper thresh = ', upperthreshold, ')')

  }

  return(as.character(meantiter))

}

#' Combine titers
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

  if(length(titers) == 0 | isTRUE(all(titers %in% c('', '*'))) |isTRUE(all(is.na(titers)))) return("*")
  if(length(titers) == 1) return(titers)
  titers <- titers[!is.na(titers) & titers != "*"]
  paste(titers, collapse = "/")
}

#####
# handling exps / expdbs
#####


#' Merge from exp
#'
#' Get merged titer table from single experiment.
#'
#' @param exp char: like "240/480/240"
#' @param expect_repeats bool: whether repeat titers are expected
#' @param merge_passage bool: not implemented
#' @param threshold
#'
#' @return array of char: merged titer table
#'
#' @rdname
#' @export
#'
#' @examples
exper.merge <- function(
  exp,
  expect_repeats = FALSE,
  duplicate_titers = "mean",
  merge_passage = FALSE,
  lowerthreshold = 10,
  upperthreshold = Inf,
  threshold_masking = F
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
  if (duplicate_titers == "mean") table_merge <- apply(table_array, c(1,2), mergeTiters, lowerthreshold, upperthreshold, threshold_masking)
  if (duplicate_titers == "combine") table_merge <- apply(table_array, c(1,2), combineTiters)

  rownames(table_merge) <- all_agids
  colnames(table_merge) <- all_srids

  # Add HAU data as an attribute
  if(sum(hau != "") > 0){
    attr(table_merge, "hau") <- hau
  }

  # Return the result
  table_merge

}

#' Merge from expdb
#'
#' Get merged titer table from experimental database.
#'
#' @param expdb char
#' @param duplicate_titers char: "mean" or "combine"
#' @param threshold
#'
#' @return array of char: merged titer table
#' @export
#'
#' @examples
expdb.merge <- function(
  expdb,
  duplicate_titers = "mean",
  lowerthreshold = 10,
  upperthreshold = Inf
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

  if(duplicate_titers == "mean")         { table_merge <- apply(table_array, c(1, 2), mergeTiters, lowerthreshold, upperthreshold)   }
  else if(duplicate_titers == "combine") { table_merge <- apply(table_array, c(1, 2), combineTiters) }
  else { stop("'duplicate_titers' must be one of 'mean' or 'combine'") }

  rownames(table_merge) <- all_agids
  colnames(table_merge) <- all_srids
  table_merge
}









