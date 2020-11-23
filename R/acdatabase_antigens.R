# Function to generate a new database
#' @export
agdb.new <- function(db = list()){

  # Fetch antigen ids
  agids <- lapply(db, function(x){ x$id })

  # Convert any object to environments
  db <- lapply(db, function(x){
    x <- as.environment(x)
    class(x) <- c("acdatabase.ag", "acdatabase.entry", "environment")
    x
  })

  # Link parent environments
  lapply(db, function(db.ag){
    parent.env(db.ag) <- emptyenv()
    if(!is.null(db.ag$parent_id)){
      matched_parent <- which(agids == db.ag$parent_id)
      if(length(matched_parent) >  1) stop(sprintf("Multiple parents '%s' found for '%s'", db.ag$parent_id, db.ag$id))
      if(length(matched_parent) == 0) stop(sprintf("Parents '%s' not found for '%s'", db.ag$parent_id, db.ag$id))
      db.ag[[".parent"]] <- db[[matched_parent]]
    }
  })

  # Link alteration environments
  lapply(db, function(db.ag){
    if(!is.null(db.ag[["alterations"]])){
      db.ag[["alterations"]] <- lapply(db.ag[["alterations"]], function(alteration){
        if(!is.null(alteration[["parent_id"]])){
          matched_parent <- which(agids == alteration[["parent_id"]])
          alteration[[".parent"]] <- db[[matched_parent]]
        }
        alteration
      })
    }
  })

  # Return the database
  db

}




# Append an antigen to an antigen db
#' @export
agdb.append <- function(agdb, ag){

  # Check for duplicate antigens
  duplicates <- do.call(
    agdb.find,
    c(list(db = agdb), as.list(ag))
  )

  # Stop if duplicates are found
  if(sum(duplicates) > 0) stop("A duplicate antigen was found")

  # Otherwise append as normal
  append(agdb, ag)

}






# Function to generate an id
#' @export
agdb.id <- function(db, id = NULL, excluded_ids = NULL){

  # Add database ids to excluded ids
  if(!missing(db)) excluded_ids <- c(unlist(db%$%id), excluded_ids)

  # Keep cycling until a unique id is found
  while(is.null(id) || id %in% excluded_ids){
    id <- paste(
      c(LETTERS, 0:9)[sample(1:36, 6, TRUE)],
      collapse = ""
    )
  }

  # Return the id
  id

}




# Functions to generate records for the database
#' @export
agdb.ag <- function(id,
                    strain,
                    long,
                    aliases,
                    wildtype,
                    type,
                    subtype,
                    lineage,
                    isolation,
                    genes,
                    parent_id,
                    alterations,
                    passage,
                    comments,
                    groups,
                    meta,
                    agdb,
                    .parent = NULL){

  ag <- new.env(parent = emptyenv())
  class(ag) <- c("acdatabase.ag", "acdatabase.entry", "environment")

  if(!missing(id))        { ag$id       <- id         } else { stop("ID must be provided") }
  if(!missing(strain))    { ag$strain   <- strain     }
  if(!missing(long))      { ag$long     <- long       }

  if(!missing(aliases)) { ag$aliases <- aliases }

  if(!missing(wildtype))  { ag$wildtype <- wildtype            }
  if(!missing(type))      { ag$type     <- type                }
  if(!missing(subtype))   { ag$subtype  <- subtype             }
  if(!missing(lineage))   { ag$lineage  <- tolower(lineage)    }

  if(!missing(isolation))   { ag$isolation <- isolation }
  if(!missing(genes))       { ag$genes     <- genes     }

  if(!missing(alterations)) { ag$alterations <- alterations }
  if(!missing(passage))     { ag$passage     <- passage     }

  if(!missing(comments))    { ag$comments  <- comments  }
  if(!missing(groups))      { ag$groups    <- groups  }

  if(!missing(meta))        { ag$meta      <- meta }

  if(!missing(parent_id))   {
    if(missing(agdb)) {
      stop("You must provide the agdb object when specifying a parent_id so that the correct parent can be linked to the antigen.")
    } else {
      ag$parent_id   <- parent_id
      agids          <- lapply(agdb, function(x){ x$id })
      if(sum(agids == parent_id) != 1) stop("No corresponding parent antigen found")
      parent.env(ag) <- agdb[[which(agids == parent_id)]]
    }
  }

  # Link a parent antigen if specified
  if(!is.null(.parent)){
    ag$parent_id <- .parent$id
    parent.env(ag) <- .parent
  }

  ag

}



# Isolation data
#' @export
agdb.isolation <- function(id,
                           date,
                           cell,
                           location,
                           continent){

  isolation <- list()
  if(!missing(id))        { isolation$id        <- id                 }
  if(!missing(date))      { isolation$date      <- date               }
  if(!missing(cell))      { isolation$cell      <- cell               }
  if(!missing(location))  { isolation$location  <- tolower(location)  }
  if(!missing(continent)) { isolation$continent <- tolower(continent) }
  isolation

}





# Genetic data
#' @export
agdb.genes <- function(gene,
                       sequence,
                       clade,
                       cluster){

  genes <- list()
  if(!missing(gene))     { genes$gene     <- gene     }
  if(!missing(sequence)) { genes$sequence <- sequence }
  if(!missing(clade))    { genes$clade    <- clade    }
  if(!missing(cluster))  { genes$cluster  <- cluster  }
  genes

}





# Alteration data
#' @export
agdb.alterations <- function(gene,
                             parent_id,
                             substitutions){

  alterations <- list()
  if(!missing(gene))          { alterations$gene          <- gene      }
  if(!missing(parent_id))     { alterations$parent_id     <- parent_id }
  if(!missing(substitutions)) { alterations$substitutions <- substitutions              }
  alterations

}





# Passage data
#' @export
agdb.passage <- function(history,
                         egg,
                         cell,
                         details,
                         comments){

  passage <- list()
  if(!missing(history))  { passage$history  <- history  }
  if(!missing(egg))      { passage$egg      <- egg      }
  if(!missing(cell))     { passage$cell     <- cell     }
  if(!missing(details))  { passage$details  <- details  }
  if(!missing(comments)) { passage$comments <- comments }
  passage

}




# Write a database
#' @export
write.agdb <- function(db, file){

  # Check for unique ids
  if(sum(duplicated(collate(db%$%id))) > 0) {
    stop("Database contains duplicate ids")
  }

  # Set variables
  nadditions <- length(db)
  ndeletions <- 0
  nedits     <- 0

  # Convert to list
  db <- acdb.as.list(db)

  # Read any original file
  if(file.exists(file)){

    dbo <- read.agdb(file)
    dbo <- acdb.as.list(dbo)
    message("")
    message("Existing database updated")

    newids  <- vapply(db, function(x){ x$id }, character(1))
    origids <- vapply(dbo, function(x){ x$id }, character(1))

    additions <- newids[!newids %in% origids]
    deletions <- origids[!origids %in% newids]

    edits <- vapply(newids[newids %in% origids], function(id){

      # Get antigens to compare
      ag1 <- db[[which(newids == id)]]
      ag2 <- dbo[[which(origids == id)]]

      # Compare alphabetically ordered values
      !identical(
        ag1[order(names(ag1))],
        ag2[order(names(ag2))]
      )
    }, logical(1))

    nadditions <- length(additions)
    ndeletions <- length(deletions)
    nedits     <- sum(edits)

    message(sprintf("%s additions", nadditions))
    message(sprintf("%s deletions", ndeletions))
    message(sprintf("%s edits",     nedits))
    message("")

  }

  # Write out the database
  jsonlite::write_json(
    x          = outbox.agdb.database(db),
    path       = file,
    pretty     = 4,
    auto_unbox = TRUE
  )

  # Return data invisibly
  invisible(
    list(
      additions = nadditions,
      deletions = ndeletions,
      edits     = nedits
    )
  )

}





# Read a database
#' @export
read.agdb <- function(file){

  agdb.new(
    jsonlite::read_json(
      path = file,
      simplifyVector    = TRUE,
      simplifyDataFrame = FALSE,
      simplifyMatrix    = FALSE
    )
  )

}





# Match an antigen
#' @export
agdb.find <- function(
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



# Match and extract antigens
#' @export
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



#' @export
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
