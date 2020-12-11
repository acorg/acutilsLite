##########################
#
#    Square
#
##########################



#' @export
squaretiters.toLong <- function(
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

#' @export
squaretiters.addNames <- function(
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

##########################
#
#    Long
#
##########################



#' Parse an experiment to tibble
#'
#' Parse an experiment record into a long tibble containing all results.
#'
#' @param exp
#'
#' @return tibble
#'
#' @rdname
#' @export
#'
#' @examples
exper.toLongTibble <- function(exp){

  # Fetch list of long-format results tibbles
  results_list <- lapply(seq_along(exp$results), function(x){
    result <- exp$results[[x]]
    titers <- result$titers
    rownames(titers) <- result$antigen_ids
    colnames(titers) <- result$serum_ids
    tibble::as_tibble(
      titers,
      rownames = "ag"
    ) %>%
      tidyr::pivot_longer(
        -ag,
        names_to = "sr",
        values_to = "titer"
      ) %>% dplyr::mutate(
        exp.result = x,
        assay = result$assay
      )
  })

  # Merge into one long tibble
  do.call(dplyr::bind_rows, results_list) %>%
    dplyr::mutate(
      exp = exp$id
    )

}

#' Parse an expdb to tibble
#'
#' Parse an experiment record into a long tibble containing all results.
#'
#' @param expdb
#'
#' @return tibble
#' @export
#'
#' @examples
expdb.toLongTibble <- function(expdb){

  # Fetch list of long-format experiment tibbles
  exps_list <- lapply(expdb, exper.toLongTibble)

  # Merge into one long tibble
  exps_tibble <- do.call(dplyr::bind_rows, exps_list)

  # Return tibble alongside ag, sr and experiment references
  exps_tibble %>%
    dplyr::mutate(
      exp_records = acdb.getIDs(exp, expdb),
      ag_records = acdb.getIDs(ag, agdb),
      sr_records = acdb.getIDs(sr, srdb),
      srag_records = srdb.homologousAntigens(sr_records, agdb),
    )

}






#' Add tibble plotting columns
#'
#' Add columns useful for plotting titer data
#'
#' @param titertable_long tibble
#'
#' @return tibble
#' @export
#'
#' @examples
longtiters.plotdata <- function(
  titertable_long
){

  # Fill the additional data needed for plotting
  titertable_long %>%
    dplyr::mutate(
      ag_cluster  = agdb.clusters(ag_records),
      ag_year     = agdb.year(ag_records),
      ag_long     = acdb.long(ag_records),
      ag_short    = agdb.short(ag_records),
      ag_plotid   = factor(
        ag,
        levels = unique(ag[order(ag_cluster, ag_year)])
      )
    ) %>%
    dplyr::mutate(
      sr_cluster  = agdb.clusters(srag_records),
      sr_year     = agdb.year(srag_records),
      sr_long     = acdb.long(sr_records),
      sr_short    = agdb.short(srag_records),
      sr_plotid   = factor(
        sr,
        levels = unique(sr[order(sr_cluster, sr_year)])
      ),
      logtiter = Racmacs::titer_to_logtiter(titer)
    ) -> plotdata

  if("titer" %in% colnames(plotdata)){
    plotdata <- dplyr::mutate(
      plotdata,
      logtiter = Racmacs::titer_to_logtiter(titer)
    )

  }

  # Return plotdata
  plotdata

}

#'@export
longtiters.splitSubstitutions <- function(titerlong){
  if (! ('ag_substitutions') %in% names(titerlong)) stop("Need a column called 'ag_substitutions'")

  splitsubs = purrr::transpose(subs.split.list(titerlong$ag_substitutions))

  titerlong$ag_subs_from = unlist_safe(splitsubs$from)
  titerlong$ag_subs_at = unlist_safe(splitsubs$at)
  titerlong$ag_subs_to = unlist_safe(splitsubs$to)

  return(titerlong)
}


#'@export
longtiters.orderSera = function(longTiters, antigens = F){
  sr_rank = srdb.getRank(longTiters$sr_records, agdb)
  longTiters %>% mutate(sr_rank = unlist(acdb.applyFunction(sr_records, function(sr)return(sr_rank[[sr$id]])))) -> longTiters

  if (antigens){
    ag_rank_list = agdb.getRank(longTiters$ag_records)

    print(length(unlist(acdb.applyFunction(longTiters$ag_records, function(ag)return(ag_rank_list[[ag$id]])))))
    print(length((acdb.applyFunction(longTiters$ag_records, function(ag)return(ag_rank_list[[ag$id]])))))

    longTiters %>% mutate(ag_rank = unlist(acdb.applyFunction(ag_records, function(ag)return(ag_rank_list[[ag$id]])))) -> longTiters
    #longTiters = group_by(longTiters, ag_rank, sr_rank)
    longTiters = arrange(longTiters, ag_rank, sr_rank, .by_group = F)
    longTiters = tibble.factorize(longTiters, c('ag','ag_short','ag_long','sr','sr_short','sr_long'))
    return(select(longTiters, -c(sr_rank, ag_rank)))
  }
  longTiters = arrange(longTiters, sr_rank, .by_group = T)
  longTiters = tibble.factorize(longTiters, c('ag','ag_short','ag_long','sr','sr_short','sr_long'))

  return(select(longTiters, -c(sr_rank)))

}

#'@export
tibble.factorize <- function(tib, columns){
  for (col in columns){
    tib[[col]] = factor(tib[[col]], levels = unique(tib[[col]]))
  }

  return(tib)
}



#' Summarise by serum cluster
#'
#' @param plotdata tibble
#'
#' @return tibble
#' @export
#'
#' @examples
summarise_plotdata_by_sr_cluster <- function(plotdata){

  mean_titer <- function(titers){

    logtiters <- Racmacs::titer_to_logtiter(titers)
    num_measured <- sum(!is.na(logtiters))
    logtiter_mean <- mean(logtiters, na.rm = T)
    sd_logtiters <- 1
    conf.level <- 0.678

    list(
      estimate = logtiter_mean,
      ci = logtiter_mean + (stats::qnorm(c(0.5-conf.level/2, 0.5+conf.level/2), sd = sd_logtiters) / sqrt(num_measured))
    )

  }

  mean_titertype <- function(titers){
    "measured"
  }

  # Get antigen meta data
  plotdata %>%
    dplyr::ungroup() %>%
    dplyr::select("ag", dplyr::starts_with("ag_")) %>%
    dplyr::distinct() -> ag_meta

  if(sum(duplicated(ag_meta$ag)) > 0){
    stop("There was a problem extracting antigen meta data")
  }

  # Create a titer list if not already there
  if(!"titer.list" %in% plotdata){
    plotdata[["titer.list"]] <- as.list(plotdata$titer)
  }

  # Summarise by sr cluster
  plotdata %>% dplyr::group_by(
    ag, sr_cluster, add = TRUE
  ) %>%
    dplyr::summarise(
      logtiter = mean_titer(titer.list)$estimate,
      logtiter_ci_lower = mean_titer(titer.list)$ci[1],
      logtiter_ci_upper = mean_titer(titer.list)$ci[2],
      titer.list = list(unlist(titer.list)),
      titer_type = mean_titertype(titer.list),
      sr_color = unique(sr_color)
    ) %>%
    dplyr::ungroup(
      ag
    ) %>%
    dplyr::mutate(
      sr = sr_cluster,
      sr_plotid = sr_cluster
    ) -> plotdata

  # Replace antigen meta data
  dplyr::left_join(
    plotdata, ag_meta
  )

}


#' Summarise by antigen cluster
#'
#' @param plotdata tibble
#'
#' @return tibble
#' @export
#'
#' @examples
summarise_plotdata_by_ag_cluster <- function(plotdata){

  mean_titertype <- function(titers){
    "measured"
  }

  # Get antigen meta data
  plotdata %>%
    dplyr::ungroup() %>%
    dplyr::select("sr", dplyr::starts_with("sr_")) %>%
    dplyr::distinct() -> sr_meta

  if(sum(duplicated(sr_meta$sr)) > 0){
    stop("There was a problem extracting serum meta data")
  }

  # Summarise by sr cluster
  plotdata %>% dplyr::group_by(
    sr, ag_cluster, add = TRUE
  ) %>%
    dplyr::summarise(
      logtiter = mean(logtiter, na.rm = T),
      titer_type = mean_titertype(logtiter),
      ag_color = unique(ag_color)
    ) %>%
    dplyr::ungroup(
      sr
    ) %>%
    dplyr::mutate(
      ag = ag_cluster,
      ag_plotid = ag_cluster
    ) -> plotdata

  # Replace antigen meta data
  dplyr::left_join(
    plotdata, sr_meta
  )

}


#' this converts logtiters with less thans to logtiters where <3 -> 2
#'@export
logtiter.toPlot <- function (titers)
{
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
  if ('morethan' %in% titer_types) print(titers)
  log_titers
}

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
  if ('morethan' %in% titer_types) print(titers)

  log_titers
}
