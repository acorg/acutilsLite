
#' Printing an acdb entry
#'
#' Formatted print for acdatabase entries
#'
#' @param db.ag list
#'
#' @return character
#' @export
#'
#' @examples
print.acdatabase.entry <- function(db.ag, ...){

  # Print short descriptor
  cat(crayon::green(sprintf("<acdatabase.ag: %s>\n", db.ag$id)))


  cat("\n")
  env      <- db.ag
  envnames <- c()
  envnum   <- 1
  colfn    <- crayon::reset

  while(length(env) > 0){

    listvals <- as.list(env)
    listvals <- listvals[!names(listvals) %in% envnames]
    envnames <- c(envnames, names(listvals))

    output <- utils::capture.output(print(listvals))
    output[output != ""] <- paste0(rep("..", envnum - 1), output[output != ""])
    output <- paste(output, collapse = "\n")

    cat(colfn(output))
    cat("\n")

    env    <- parent.env(env)
    envnum <- envnum + 1
    colfn  <- crayon::silver

    # Only print once unless called directly as part of print statement
    if(!identical(parent.frame(), .GlobalEnv)) break

  }

  # Return the antigen invisibly
  invisible(db.ag)

}


#' Subsetting for acdb entries
#'
#' @export
#'
#' @examples
`[.acdatabase.entry` <- function(db.ag, n){
  val        <- lapply(n, get0, envir = db.ag)
  names(val) <- n
  val
}

#' Subsetting for acdb entries
#'
#' @export
#'
#' @examples
`$.acdatabase.entry` <- function(ag, selector){
  get0(selector, envir = ag)
}




#' Subsetting for databases
#'
#' Returns field for all database entries. Multilayer subsetting is also supported with `db%$%field$subfield`
#'
#' @export
#'
#' @examples
`%$%` <- function(db, val){
  val <- as.character(match.call()[3])
  val <- strsplit(val, "$", fixed = TRUE)[[1]]
  lapply(db, function(ag){
    if(is.environment(ag)){
      result <- get0(val, envir = ag)
    } else {
      result <- ag[[val]]
    }
    for(v in val[-1]){
      result <- result[[v]]
    }
    # result <- ag
    # for(v in val){
    #   result <- result[[v]]
    # }
    result
  })
}


#' Subsetting for databases
#'
#' Returns field for all database entries. Like `%$%`, except variable val is supported.
#'
#' @export
#'
#' @examples
acdb.slice <- function(db, val){
  val <- strsplit(val, "$", fixed = TRUE)[[1]]
  lapply(db, function(ag){
    if(is.environment(ag)){
      result <- get0(val, envir = ag)
    } else {
      result <- ag[[val]]
    }
    for(v in val[-1]){
      result <- result[[v]]
    }
    # result <- ag
    # for(v in val){
    #   result <- result[[v]]
    # }
    result
  })
}


#' Assignment for databases
#'
#' Assigns value to field for all database entries
#'
#' @export
#'
#' @examples
`%$%<-` <- function(tmp, ..., value){

  attribute <- as.character(match.call()[3])
  lapply(tmp, function(t){
    eval(
      parse(
        text = paste0(
          "t$",
          gsub(".", "$", attribute, fixed = TRUE),
          " <- value"
        )
      )
    )
    t
  })

}

#' Convert database entries to lists
#'
#' Antigen entries become lists of properties rather than environments `.parent` entries are removed.
#'
#' @export
#'
#' @examples
acdb.as.list <- function(db){
  lapply(db, as.list)
}

#' Convert database entry to list
#'
#' Antigen entries become lists of properties rather than environments. `.parent` entries are removed.
#'
#' @export
#'
#' @examples
as.list.acdatabase.ag <- function(record, ...){
  recordlist <- as.list.environment(record)

  # remove .parent entries under alterations (as.list.environment automatically removes top level .parent)
  if(!is.null(recordlist[["alterations"]])){
    recordlist[["alterations"]] <- lapply(recordlist[["alterations"]], function(alteration){
      parent_link <- which(names(alteration) == ".parent")
      print(parent_link)
      if(length(parent_link) > 0){
        alteration <- alteration[-parent_link]
      }
      alteration
    })
  }
  recordlist
}



