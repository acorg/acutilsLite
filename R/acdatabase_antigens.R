
#' Antigen database
#'
#'@details
#' The antigen database is ultimately stored as a json file with fields:
#' 1. id
#' 2. strain
#' 3. long
#' 4. aliases
#' 5. wildtype
#' 6. type
#' 7. subtype
#' 8. lineage
#' 9. isolation
#' 10. genes
#' 11. parent_id
#' 12. alterations
#' 13. passage
#' 14. comments
#' 15. groups
#' 16. meta
#'
#' Not all fields are filled for all antigens.
#'
#' In R, this is converted to a list of environments (see ?ag) with the above fields, plus:
#' 1. A ag$.parent entry, which is the database entry for the parent antigen (defined by the backbone).
#' 2. ag$alteration$.parent entries for each gene that is not inherited from ag$.parent (for recombinant viruses).
#'
#' @name agdb
NULL



#' Antigen
#'
#'@details
#' An antigen entry (ag) is an environment containing some of:
#'
#'
#'\tabular{ll}{
#' id   \tab  Antigen id: 6 random capital letters and numbers \cr
#' strain \tab  Rarely used. \cr
#' long \tab Full strain name. \cr
#' aliases \tab Other names known to be used to refer to this antigen. \cr
#' wildtype \tab Boolean: whether the strain is wildtype. \cr
#' type \tab Influenza virus type (A, B...) \cr
#' subtype \tab Influenza virus subtype (H3N2, H1N1...) \cr
#' lineage \tab Victoria etc \cr
#' isolation \tab list containing isolation id, location, date, cell, continent \cr
#' genes \tab gene sequence for HA/NA \cr
#' parent_id \tab id of parent (by backbone) \cr
#' .parent \tab the database entry for the parent antigen (by backbone) \cr
#' alterations \tab alterations to parent (substitutions or gene transplants), organised by gene. Substitutions formatted as X123Y, transgenes identified by parent id. alterattion$.parent is the database entry for the gene origin when the alteration is a transplant.  \cr
#' passage \tab passaging history \cr
#' comments \tab Extra comments - eg which study the antigen comes from \cr
#' groups \tab Useful annotations. For example: wildtype, mutatant, gen 1 root, gen 2 mutant, 3C.3A etc.  \cr
#' meta \tab Some extra data about the antigen: cluster, mutant generation number, whether it is a reference strain \cr
#' }
#'
#' Notably, the parent environment is always the empty environment; $.parent is the backbone-wise parent antigen entry; alteration$.parent is the origin antigen for a transplanted gene
#'
#' It has classes c("acdatabase.ag", "acdatabase.entry", "environment")
#' @name ag
NULL


# load test database
agdb.loadTest <- function(){
  files = list.files(recursive = T)
  split_files = lapply(files, function(f){strsplit(f, split = '/')} )
  mask = unlist(lapply(split_files, function(file)file[[1]][[length(file[[1]])]] == 'agdb_h3_small.json'))
  db_path = files[which.max(mask)]
  agdb.test <<- read.agdb(db_path)

  invisible(read.agdb(db_path))
}

##########################
#
#    READ AND WRITE
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




#' Write database to file
#'
#' Writes an agdb to specified location. If a database with the same name already exists at that location, provides summary of additions, deletions and edits to database
#'
#' @param agdb list
#' @param file char
#'
#' @return list
#' @export
#'
#' @examples
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

# Set the name order
ag.nameorder <- c(
  "id",
  "parent_id",
  "long",
  "strain",
  "wildtype",
  "type",
  "subtype",
  "lineage",
  "isolation",
  "genes",
  "alterations",
  "passage",
  "comments",
  "groups",
  "meta"
)

#' Mark which properties should be unboxed
#'
#' Helper for write.agdb
#'
#' @param db.database list
#'
#' @return list
#'
#' @examples
outbox.agdb.database <- function(db.database){

  for(n in seq_along(db.database)){

    if(!is.null(db.database[[n]][["aliases"]])) {
      db.database[[n]][["aliases"]] <- I(db.database[[n]][["aliases"]])
    }

    if(!is.null(db.database[[n]][["groups"]])) {
      db.database[[n]][["groups"]] <- I(db.database[[n]][["groups"]])
    }

    if(!is.null(db.database[[n]][["passage"]][["history"]])){
      db.database[[n]][["passage"]][["history"]] <- I(db.database[[n]][["passage"]][["history"]])
    }

    if(!is.null(db.database[[n]][["alterations"]])) {

      for(m in seq_along(db.database[[n]][["alterations"]])){
        if(!is.null(db.database[[n]][["alterations"]][[m]][["substitutions"]])){
          db.database[[n]][["alterations"]][[m]][["substitutions"]] <- I(db.database[[n]][["alterations"]][[m]][["substitutions"]])
        }
      }

    }

    if(!is.null(db.database[[n]][["meta"]][["egg-passage-mutations"]][["HA"]])) {
      db.database[[n]][["meta"]][["egg-passage-mutations"]][["HA"]] <- I(db.database[[n]][["meta"]][["egg-passage-mutations"]][["HA"]])
    }
    if(!is.null(db.database[[n]][["meta"]][["egg-passage-mutations"]][["NA"]])) {
      db.database[[n]][["meta"]][["egg-passage-mutations"]][["NA"]] <- I(db.database[[n]][["meta"]][["egg-passage-mutations"]][["NA"]])
    }

    # Order the names
    db.database[[n]] <- db.database[[n]][
      order(
        match(
          names(db.database[[n]]),
          ag.nameorder
        )
      )
    ]

  }

  db.database

}




##########################
#
#    CONSTRUCTORS
#
##########################




#' Produce an antigen entry for the database
#'
#' Produces an entry for antigen databases (see ?agdb) from the passed data.
#'
#'
#' @return character
#' @export
#'
#' @examples
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
      ag[[".parent"]] <- agdb[[which(agids == parent_id)]]
    }
  }

  # Link a parent antigen if specified
  if(!is.null(.parent)){
    ag$parent_id <- .parent$id
    ag[[".parent"]] <- .parent
  }

  ag

}








#' Create antigen database
#'
#' Creates an antigen database from a list of antigen entries. For information on database structure, see ?agdb.
#'
#' @param db list
#'
#' @return environment
#' @export
#'
#' @examples
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








#' Generate an id
#'
#' Generates a random id (6 capital letters and numbers) which is not seen in the provided agdb or excluded_ids
#'
#' @param agdb list
#' @param id should be left NULL
#' @param excluded_ids char (optional): a vector of ids which should not be assigned
#'
#' @return character
#' @export
#'
#' @examples
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



#' Constructs  an antigen isolation entry
#'
#' Constructs an antigen isolation entry. All params are optional
#'
#' @param id char: the isolation id
#' @param date char: the year of the isolation
#' @param cell char: cell or eg passaging?
#' @param location char#
#' @param continent char
#'
#' @return list
#' @export
#'
#' @examples
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



#' Constructs  an antigen gene entry
#'
#' Constructs an antigen gene entry. All params are optional
#'
#' @param gene char: gene name (HA or NA)
#' @param sequence char
#' @param clade char
#' @param cluster char
#'
#' @return list
#' @export
#'
#' @examples
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






#' Constructs  an antigen alterations entry
#'
#' Constructs an antigen alterations entry. All params are optional
#'
#' @param gene char: gene name (HA or NA)
#' @param parent_id char: the id of the origin virus for a transgene
#' @param substitutions char: vector containing substitutions of form 'X123Y'
#' #'
#' @return list
#' @export
#'
#' @examples
agdb.alterations <- function(gene,
                             parent_id,
                             substitutions){

  alterations <- list()
  if(!missing(gene))          { alterations$gene          <- gene      }
  if(!missing(parent_id))     { alterations$parent_id     <- parent_id }
  if(!missing(substitutions)) { alterations$substitutions <- substitutions              }
  alterations

}






#' Constructs  an antigen passage entry
#'
#' Constructs an antigen passage entry. All params are optional
#'
#' @param history char: vector of known passage history
#' @param egg char: not used
#' @param cell char: cell or egg passaging
#' @param details
#' @param comments
#' @return list
#' @export
#'
#' @examples
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















##########################
#
#    MODIFY
#
##########################

#
# modify single ag entry
#

#
# modify agdb
#

#' Append an antigen to an agdb
#'
#'
#' @param agdb list
#' @param ag environment
#'
#' @return list
#' @export
#'
#' @examples
agdb.append <- function(agdb, ag){

  # Check for duplicate antigens
  duplicates <- do.call(
    acdb.find,
    c(list(db = agdb), as.list(ag))
  )

  # Stop if duplicates are found
  if(sum(duplicates) > 0) stop("A duplicate antigen was found")

  # Otherwise append as normal
  append(agdb, ag)
}
