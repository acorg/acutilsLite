


#' Experiment database
#'
#'@details
#'
#' A list of experiments: see ?ac.exp
#'
#' @name expdb
NULL



#' Experiment
#'
#'@details
#' An experiment entry (exp) is an environment containing some of:
#'
#'
#'\tabular{ll}{
#' name   \tab   \cr
#' id \tab  6 random capitals and numerics \cr
#' description \tab  \cr
#' comments \tab  \cr
#' results \tab See below \cr
#' }
#'
#' The results field contains a list, each with the assay results from a single original titer table from a collaborating group. Each of these entries is itself a list, containing:
#'#'\tabular{ll}{
#' lab   \tab  Where the experiment was carried out \cr
#' assay \tab  HI or FRA \cr
#' conducted_by \tab Who carried out the assay \cr
#' date \tab \cr
#' comments \tab  \cr
#' file \tab Name of the original titer table file in databases repo  \cr
#' antigen_ids \tab  \cr
#' serum_ids \tab  \cr
#' titers \tab Titer table wihth '*' denoting a missing value \cr
#' meta \tab  \cr
#' }
#' Notably, the parent environment is always the empty environment; $.parent is the backbone-wise parent antigen entry; alteration$.parent is the origin antigen for a transplanted gene
#'
#' It has classes c("acdatabase.ag", "acdatabase.entry", "environment")
#' @name ac.exp
NULL




exp.nameorder <- c(
  "name",
  "id",
  "description",
  "comments",
  "results"
)

results.nameorder <- c(
  "lab",
  "assay",
  "conducted_by",
  "date",
  "comments",
  "file",
  "antigen_ids",
  "serum_ids",
  "titers",
  "meta"
)



##########################
#
#    READ AND WRITE
#
##########################



# Read a database
#' @export
read.expdb <- function(file){

  expdb.new(
    jsonlite::read_json(
      path = file,
      simplifyVector    = TRUE,
      simplifyDataFrame = FALSE,
      simplifyMatrix    = FALSE
    )
  )

}

# Read all databases in directory
#' @export
read.expdbs <- function(root_dir){
  exp_dirs <- list.files(file.path(root_dir, "experiments"), full.names = T)
  expdbs <- lapply(exp_dirs, function(exp_dir){
    if(file.exists(file.path(exp_dir, "results.json"))){
      read.expdb(file.path(exp_dir, "results.json"))
    }
  })
  unlist(expdbs, recursive = FALSE)
}



# Write a database
#' @export
write.expdb <- function(db, file){

  # Check for unique ids
  if(sum(duplicated(db%$%id)) > 0) {
    stop("Database contains duplicate ids")
  }

  # Check results matches titer dims

  # Write out the database
  dbtext <- jsonlite::toJSON(
    x          = outbox.expdb(acdb.as.list(db)),
    path       = file,
    pretty     = 4,
    auto_unbox = TRUE
  )

  # # Tidy up ag ids, sr ids and titers
  # agids <- stringr::str_extract_all(dbtext, stringr::regex('"antigen_ids": (\\[.*?\\])', dotall = TRUE))
  # for(agid in agids[[1]]) dbtext <- gsub(agid, gsub("\n[ ]*", " ", agid), dbtext, fixed = TRUE)
  #
  # srids <- stringr::str_extract_all(dbtext, stringr::regex('"serum_ids": (\\[.*?\\])', dotall = TRUE))
  # for(srid in srids[[1]]) dbtext <- gsub(srid, gsub("\n[ ]*", " ", srid), dbtext, fixed = TRUE)

  write(dbtext, file)

}

# A function to mark what properties should not be unboxed in json
outbox.expdb <- function(db){

  for(n in seq_along(db)){

    if(!is.null(db[[n]]$aliases)) {
      db[[n]]$aliases <- I(db[[n]]$aliases)
    }

    # Order the names
    if(sum(!names(db[[n]]) %in% exp.nameorder) > 0){
      stop(sprintf("Unmatched names, '%s'", names(db[[n]])[!names(db[[n]]) %in% exp.nameorder]))
    }
    db[[n]] <- db[[n]][
      order(
        match(
          names(db[[n]]),
          exp.nameorder
        )
      )
    ]

    for(m in seq_along(db[[n]]$results)){

      # Order the names
      if(sum(!names(db[[n]]$results[[m]]) %in% results.nameorder) > 0){
        stop(sprintf("Unmatched names, '%s'", names(db[[n]]$results[[m]])[!names(db[[n]]$results[[m]]) %in% results.nameorder]))
      }
      db[[n]]$results[[m]] <- db[[n]]$results[[m]][
        order(
          match(
            names(db[[n]]$results[[m]]),
            results.nameorder
          )
        )
      ]

    }

  }

  db

}






##########################
#
#    CONSTRUCTORS
#
##########################



# Function to generate a new database
#' @export
expdb.new <- function(db = list()){

  # Convert titers to matrices
  db <- lapply(db, function(exp){
    exp$results <- lapply(exp$results, function(result){
      result$titers <- do.call(rbind, result$titers)
      result
    })
    exp
  })

  # Convert any object to environments
  db <- lapply(db, function(x){
    x <- as.environment(x)
    class(x) <- c("acdatabase.exp", "acdatabase.entry", "environment")
    x
  })


  # Return the database
  db

}


# Function to generate an id
#' @export
expdb.id <- agdb.id


# Functions to generate records for the database
#' @export
expdb.exp <- function(name,
                      id,
                      description,
                      comments,
                      results = list()){

  if(!is.null(names(results)) || class(results) != "list"){
    stop("results must be an unnamed list")
  }

  exp <- new.env(parent = emptyenv())
  class(exp) <- c("acdatabase.exp", "acdatabase.entry", "environment")

  if(!missing(name))        { exp$name        <- name         }
  if(!missing(id))          { exp$id          <- id           } else { stop("ID must be provided") }
  if(!missing(description)) { exp$description <- description  }
  if(!missing(comments))    { exp$comments    <- comments     }
  exp$results <- results

  exp

}


# Function to generate a new result
#' @export
expdb.result <- function(lab,
                         assay,
                         conducted_by,
                         date,
                         comments,
                         file,
                         antigen_ids,
                         serum_ids,
                         titers,
                         meta){

  result <- list()
  if(!missing(assay))        result$assay        <- assay
  if(!missing(conducted_by)) result$conducted_by <- conducted_by
  if(!missing(date))         result$date         <- date
  if(!missing(comments))     result$comments     <- comments
  if(!missing(file))         result$file         <- file
  if(!missing(meta))         result$meta         <- meta

  if(!missing(antigen_ids) ||
     !missing(serum_ids) ||
     !missing(titers)) {

    if(sum(!missing(antigen_ids),
           !missing(serum_ids),
           !missing(titers)) != 3){
      stop("Antigen ids, serum ids and titers must be provided together")
    }

    if(nrow(titers) != length(antigen_ids)) {
      stop(sprintf("The number of antigen ids (%s) does not equal the number of titer rows (%s)", length(antigen_ids), nrow(titers)))
    }

    if(ncol(titers) != length(serum_ids)) {
      stop(sprintf("The number of serum ids (%s) does not equal the number of titer columns (%s)", length(serum_ids), ncol(titers)))
    }

    result$antigen_ids  <- antigen_ids
    result$serum_ids    <- serum_ids
    mode(titers)  <- "character"
    result$titers <- titers

  }

  result

}


# Function to append a result
#' @export
exper.result.append <- function(exp,
                              result){

  exp$results <- append(
    exp$results,
    list(result)
  )
  exp

}



####################
#
#      MODIFY
#
####################

#################
#
# MISC
#
################


#' Get titer table from result entry
#' @export
result.titerTable <- function(results){
  titertable <- results$titers
  colnames(titertable) <- results$serum_ids
  rownames(titertable) <- results$antigen_ids
  titertable
}
