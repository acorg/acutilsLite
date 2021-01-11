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

#####
# make
#####


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



#####
# modify
#####

#' Add ag, sr, exp records from IDs to tibble
#'
#'
#' @param longtiters
#'
#' @return tibble
#' @export
#'
#' @examples
longtiters.addRecords <- function(longtiters, agdb = get_agdb(), srdb = get_srdb(), expdb = get_expdb()){
  longtiters %>% mutate(
      ag_records = acdb.getIDs(ag, agdb),
      sr_records = acdb.getIDs(sr, srdb),
      srag_records = srdb.homologousAntigens(sr_records, agdb),
      exp_records = acdb.getIDs(exp, expdb)
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
      ag_rootclade = acdb.applyFunction(.$ag_records, ag.clade, 'any')
    ) %>%
    dplyr::mutate(
      sr_cluster  = agdb.clusters(srag_records),
      sr_year     = agdb.year(srag_records),
      sr_long     = acdb.long(sr_records),
      sr_short    = agdb.short(srag_records),
      logtiter = Racmacs::titer_to_logtiter(titer),
      sr_rootclade = acdb.applyFunction(.$srag_records, ag.clade, 'any')
    ) -> plotdata

  if("titer" %in% colnames(plotdata)){
    plotdata <- dplyr::mutate(
      plotdata,
      logtiter = titer.toLog(titer)
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

  titerlong$ag_subs_from.at = paste(titerlong$ag_subs_from, titerlong$ag_subs_at)
  titerlong$ag_subs_from.at = factor(titerlong$ag_subs_from.at,
                                     levels = unique(titerlong$ag_subs_from.at[order(titerlong$ag_subs_at)]))


  titerlong$ag_subs_to = factor(titerlong$ag_subs_to, levels = aa.values())




  return(titerlong)
}

#'@export
longtiters.order <- function(longtiters, antigens = T, sera = T, by = 'clade'){
  if (by == 'clade'){
    if (antigens){
      if (!all(c('ag_rootclade', 'ag_year') %in% colnames(longtiters))) stop('Some required columns are missing - run longtiters.plotdata()')

      mutate(longtiters,
             ag_rootclade = factor(ag_rootclade, levels = h3_clade_order),
             ag_year = factor(ag_year, levels = unique(ag_year[order(as.numeric(ag_year))])),
      ) -> longtiters

      longtiters = arrange(longtiters, ag_rootclade, ag_year)
      longtiters$ag = factor(longtiters$ag, unique(longtiters$ag))


    }

    if (sera){
      if (!all(c('sr_rootclade', 'sr_year') %in% colnames(longtiters))) stop('Some required columns are missing - run longtiters.plotdata()')

      mutate(longtiters,
             sr_rootclade = factor(sr_rootclade, levels = h3_clade_order),
             sr_year = factor(sr_year, levels = unique(sr_year[order(as.numeric(sr_year))])),
      ) -> longtiters

      longtiters = arrange(longtiters, sr_rootclade, sr_year)
      longtiters$sr = factor(longtiters$sr, unique(longtiters$sr))

    }
  }
  if (by == 'cluster') stop('Not implemented')


longtiters

}



#'@export
longtiters.labelRepeats <- function(longtiters, columns = c('ag', 'sr')){
  if (missing(columns)) message("columns not specified, using default c('ag', 'sr')")
  longtiters['rep'] = 1
  diff_cols = dplyr::select(longtiters, columns)
  duplicated_uniques = dplyr::distinct(diff_cols[duplicated(diff_cols),])

  for (i in seq_len(dim(duplicated_uniques)[[1]])){
    unique = duplicated_uniques[i,]
    longtiters[apply(longtiters[,names(unique)],1 , function(r) all(r==unique)),]['rep'] = 1:sum(apply(longtiters[,names(unique)],1 , function(r) all(r==unique)))

  }

  longtiters

}



#dplyr::group_by(x, ag, sr) %>% summarise( meantiter = mean(titer))
# join
# drop varying columns
# unique

#'@export
longtiters.merge <- function(longtiters, columns = c('ag', 'sr'), repeats = 'merge', ignore_cols = c('exp.result')){
  if (missing(columns)) message("columns not specified, using default c('ag', 'sr')")


  #group_by(x, ag, sr) %>% summarise( meantiter = mean(titer))

  remove = c()
  agsr_info = select(longtiters, -dplyr::any_of(c('titer', 'logtiter', 'plotttiter')))
  for (col in colnames(agsr_info)){
    ns = agsr_info %>%
    group_by_at(columns) %>%
    select_at(c(col, columns)) %>%
    distinct() %>%
    summarise(n(), .groups = 'drop')

  if (!( all(ns[['n()']] == 1)) ) remove = c(remove, col)
  }

  if (length(remove) > 0)  message( paste(c('Removing columns:', remove), collapse = ' '))



  agsr_info = distinct(select(agsr_info, -any_of(remove)))

  group_by_at(longtiters, columns) %>%
    summarise( original_titers = combineTiters(titer),
               titer = mergeTiters(titer),
               .groups = 'drop') -> longtiters.merged


  if (any(stringr::str_detect(remove, 'colbase'))) message('Colbase titers must be recalculated')
  longtiters.merged = ungroup(longtiters.merged)
  longtiters.merged$logtiter = titer.toLog(longtiters.merged$titer)
  #longtiters.merged$plottiter = logtiter.toPlot(longtiters.merged$logtiter)
  longtiters.merged = dplyr::left_join(longtiters.merged, agsr_info, columns)

longtiters.merged
}

#'@export
print_and_capture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}

#' Colbase titers
#'@export
longtiters.colbase <- function(longtiters, longtiters.maxSet, by = 'sr'){

  if (by == 'sr'){
    sr_maxes = c()
    for(sri in unique(longtiters$sr)){
      if (sri %in% longtiters.maxSet$sr){
        sr_maxes[sri] = max(logtiter.toPlot(titer.toLog(filter(longtiters.maxSet, sr == sri)$titer)))
      }
      else{
        warning('Serum ', sr, ' not present in maxSet titer table - max titer from target titer table being used.')
        sr_maxes[sri] = max(logtiter.toPlot(titer.toLog(filter(longtiters, sr == sri)$titer)))
      }
    }

    longtiters %>%
      mutate(colbase_titer = plottiter - sr_maxes[as.character(sr)]) %>%
      mutate(titertype = Racmacs::titer_types(titer)) -> longtiters

  return(longtiters)
  }

  if (by == 'ag'){

    ag_maxes = c()
    for(agi in unique(longtiters$ag)){
      if (agi %in% longtiters.maxSet$ag){
        ag_maxes[agi] = max(logtiter.toPlot(titer.toLog(filter(longtiters.maxSet, ag == agi)$titer)))
      }
      else{
        warning('Antigen ', ag, ' not present in maxSet titer table - max titer from target titer table being used.')
        ag_maxes[agi] = max(logtiter.toPlot(titer.toLog(filter(longtiters, ag == agi)$titer)))
      }
    }


    longtiters %>%
      mutate(colbase_titer_ag = plottiter - ag_maxes[as.character(ag)]) %>%
      mutate(titertype = Racmacs::titer_types(titer)) -> longtiters

    return(longtiters)
  }
}

#' convert column(s) into factor for plotting
#'@export
tibble.factorize <- function(tib, columns){
  for (col in columns){
    tib[[col]] = factor(tib[[col]], levels = unique(tib[[col]]))
  }

  return(tib)
}

#'@export
longtiters.getAbsentMuts = function(longTiters.muts){
  if (!all(c('ag_subs_from', 'ag_subs_at', 'ag_subs_to') %in% colnames(longTiters.muts) )) stop("'ag_subs_from', 'ag_subs_at', 'ag_subs_to' columns required")
  longTiters.muts %>% select(ag_subs_from, ag_subs_at, ag_subs_to) %>% distinct() %>% transmute(from.at.to = mapply(paste , ag_subs_from, ag_subs_at, ag_subs_to)) -> subs_present
  all_possible_subs = expand.grid( ag_subs_from.at = unique(longTiters.muts$ag_subs_from.at), ag_subs_to = unique(longTiters.muts$ag_subs_to))
  all_possible_subs %>%
    mutate(present = lapply(paste(ag_subs_from.at, ag_subs_to), function(x){x %in% subs_present$from.at.to} )) %>% filter(present == F) -> subs_absent #%>% filter(present == F) -> included_subs
  subs_absent
}

#'@export
longtiters.addComparisonTiter <- function(longtiters.target, longtiters.comparator, comparator_id, agsr, newcol_name = 'comparator_titer', titer_add='plottiter'){
  longTiters.add = filter(longtiters.comparator,  (!!as.symbol(agsr)) == comparator_id)
  longtiters.target[,newcol_name] = longTiters.add[match(longtiters.target$sr, longTiters.add$sr), ][,titer_add]

  longtiters.target
}

#'@export
longtiters.subsTrafficLight <- function(longtiters.muts, srag_seq_map){

  longTiters.muts %>% mutate(srag_sequence = as.list(agdb.sequencesFromMap(srag_records, srag_seq_map))) -> longTiters.muts
  longTiters.muts %>% mutate(serum_mut = mapply(function(subs, srag_sequence){
    sub = subs.split(subs)[[1]];
    sr_aa = str_sub(srag_sequence, sub['at'], sub['at'] )
    if (sr_aa == sub['from']) return('from')
    if (sr_aa == sub['to']) return('to')
    return('neither')
  },
  longTiters.muts$ag_substitutions, longTiters.muts$srag_sequence)
  ) -> longTiters.muts

  longTiters.muts

}

#'@export
longtiters.commonColumns = function(df_list){
  cols = lapply(df_list, colnames)
  common_cols = Reduce(intersect, cols)
  lapply(df_list, function(df)df[,common_cols])
}

#####
# summarising
#####


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
