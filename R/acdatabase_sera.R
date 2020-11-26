
#' Serum database
#'
#'@details
#' The serum database is ultimately stored as a json file with fields (see `?sr` for more details):
#'
#' 1. id
#' 2. long
#' 3. strain_id
#' 4. parent_id
#' 5. species
#' 6. animal_id
#' 7. sample_date
#' 8. meta
#' 9. comments
#' 10. aliases
#'
#' Not all fields are filled for all sera
#'
#' In R, this is converted to a list of environments (see `?sr`) with the above fields, plus:
#' 1. A ar$.parent entry, which is the database entry for the parent serum
#'
#' Supports subsetting with `%$%` (see `?%$%`)
#'
#' @name srdb
NULL



#' Serum
#'
#'@details
#' A serum entry (`sr`) is an environment containing some of:
#'
#'
#'\tabular{ll}{
#' id   \tab  serum id: 6 random capital letters and numbers \cr
#' long \tab a full name for the serum. \cr
#' strain_id \tab  id of the antigen used to raise serum (for ferret sera). \cr
#' parent_id \tab id of another serum taken from the same individual. Generally used for pre/post vaccination studies. \cr
#' .parent \tab a previous serum taken from the same individual for pre/post vaccination studies \cr
#' species \tab human or ferret. \cr
#' animal_id \tab a ferret id \cr
#' sample_date \tab  \cr
#' meta \tab includes: vaccine_trial_year, age_group, age_at_vaccination etc \cr
#' comments \tab currently unused \cr
#' aliases \tab generally empty, sometimes contains strain name \cr

#' }
#'
#' The parent environment is always the empty environment; `$.parent` is a previous serum taken from the same individual for pre/post vaccination studies.
#'
#' It has classes c("acdatabase.sr", "acdatabase.entry", "environment")
#' @name sr
NULL



##########################
#
#    WORKING WITH SRDBs
#
##########################


#' Read database
#'
#' @param file char
#'
#' @return list
#' @export
#'
#' @examples
read.srdb <- function(file){

  srdb.new(
    jsonlite::read_json(
      path = file,
      simplifyVector    = TRUE,
      simplifyDataFrame = FALSE,
      simplifyMatrix    = FALSE
    )
  )

}





#' Write database to file
#'
#' Writes an agdb to specified location.
#'
#' @param agdb list
#' @param file char
#'
#' @return list
#' @export
#'
#' @examples
write.srdb <- function(db, file){

  # Check for unique ids
  if(sum(duplicated(db%$%id)) > 0) {
    stop("Database contains duplicate ids")
  }

  # Write out the database
  jsonlite::write_json(
    x          = outbox.srdb(acdb.as.list(db)),
    path       = file,
    pretty     = 4,
    auto_unbox = TRUE
  )

}

sr.nameorder <- c(
  "id",
  "long",
  "strain_id",
  "parent_id",
  "species",
  "animal_id",
  "sample_date",
  "meta",
  "comments",
  "aliases"
)

#' Mark which properties should be unboxed
#'
#' Helper for write.srdb
#'
#' @param db.database list
#'
#' @return list
#'
#' @examples
outbox.srdb <- function(db){

  for(n in seq_along(db)){

    if(!is.null(db[[n]]$aliases)) {
      db[[n]]$aliases <- I(db[[n]]$aliases)
    }

    # Order the names
    db[[n]] <- db[[n]][
      order(
        match(
          names(db[[n]]),
          sr.nameorder
        )
      )
    ]

  }

  db

}




##########################
#
#    CONSTRUCTORS
#
##########################




#' Generate an id
#'
#' Generates a random id (6 capital letters and numbers) which is not seen in the provided srdb or excluded_ids
#'
#' @param srdb list a serum databses (?srdb)
#' @param id should be left NULL
#' @param excluded_ids char (optional): a vector of ids which should not be assigned
#'
#' @return character
#' @export
#'
#' @examples
srdb.id <- agdb.id



#' Produce a serum entry for the database
#'
#' Produces an entry for serum databases (see ?srdb) from the passed data.
#'
#'
#' @return character
#' @export
#'
#' @examples
srdb.sr <- function(id,
                    long,
                    parent_id,
                    strain_id,
                    species,
                    animal_id,
                    sample_date,
                    meta,
                    aliases,
                    srdb,
                    .parent = NULL){

  sr <- new.env(parent = emptyenv())
  class(sr) <- c("acdatabase.sr", "acdatabase.entry", "environment")

  if(!missing(id))          { sr$id          <- id           } else { stop("ID must be provided") }
  if(!missing(long))        { sr$long        <- long         }
  if(!missing(strain_id))   { sr$strain_id   <- strain_id    }
  if(!missing(species))     { sr$species     <- species      }
  if(!missing(animal_id))   { sr$animal_id   <- animal_id    }
  if(!missing(sample_date)) { sr$sample_date <- sample_date  }
  if(!missing(meta))        { sr$meta        <- meta         }
  if(!missing(aliases))     { sr$aliases     <- aliases      }

  if(!missing(parent_id))   {
    if(missing(srdb)) {
      stop("You must provide the srdb object when specifying a parent_id so that the correct parent can be linked to the antigen.")
    } else {
      sr$parent_id   <- parent_id
      srids          <- lapply(srdb, function(x){ x$id })
      if(sum(srids == parent_id) != 1) stop("No corresponding parent serum record found")
      parent.env(sr) <- srdb[[which(srids == parent_id)]]
    }
  }

  # Link a parent serum if specified
  if(!is.null(.parent)){
    sr$parent_id <- .parent$id
    parent.env(sr) <- .parent
  }

  sr

}



#' Create serum database
#'
#' Creates a serum database from a list of antigen entries. For information on database structure, see ?srdb.
#'
#' @param db list
#'
#' @return environment
#' @export
#'
#' @examples
srdb.new <- function(db = list()){

  # Fetch serum ids
  srids <- lapply(db, function(x){ x$id })

  # Convert any object to environments
  db <- lapply(db, function(x){
    x <- as.environment(x)
    class(x) <- c("acdatabase.sr", "acdatabase.entry", "environment")
    x
  })

  # Link parent environments
  lapply(db, function(db.sr){
    if(!is.null(db.sr$parent_id)){
      db.sr[['.parent']] <- db[[which(srids == db.sr$parent_id)]]
    }
  })

  # Return the database
  db

}



