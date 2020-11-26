

#' Parse an experiment to tibble
#'
#' Parse an experiment record into a long tibble containing all results.
#'
#' @param exp
#'
#' @return tibble
#' @export
#'
#' @examples
exp.longTibble <- function(exp){

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
expdb.longTibble <- function(expdb){

  # Fetch list of long-format experiment tibbles
  exps_list <- lapply(expdb, exp.longTibble)

  # Merge into one long tibble
  exps_tibble <- do.call(dplyr::bind_rows, exps_list)

  # Return tibble alongside ag, sr and experiment references
  exps_tibble %>%
    dplyr::mutate(
      exp_records = acdb.getIDs(expdb, exp),
      ag_records = acdb.getIDs(agdb, ag),
      sr_records = acdb.getIDs(srdb, sr),
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
titerlong.plotdata <- function(
  titertable_long
){

  # Fill the additional data needed for plotting
  titertable_long %>%
    dplyr::mutate(
      ag_cluster  = agdb.clusters(ag_records),
      ag_year     = agdb.year(ag_records),
      ag_long     = acdb.Long(ag_records),
      ag_short    = agdb.Short(ag_records),
      ag_plotid   = factor(
        ag,
        levels = unique(ag[order(ag_cluster, ag_year)])
      )
    ) %>%
    dplyr::mutate(
      sr_cluster  = agdb.clusters(srag_records),
      sr_year     = agdb.year(srag_records),
      sr_long     = acdb.Long(sr_records),
      sr_short    = agdb.Short(srag_records),
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




