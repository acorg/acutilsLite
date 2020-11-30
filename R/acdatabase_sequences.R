
#' Split substitutions to get from, at, to
#' @export
subs.split <- function(subs){
  if (is.null(subs)) return( list(c(from = NA, at =  NA, to = NA) ))
  lapply(subs, function(sub){

    if (is.list(subs)) sub = unlist(sub)


    first = NA
    last = NA

    if (is.null(sub)) return(c(from = first, at =  NA, to = last))


    if (stringr::str_sub(sub, 1, 1) %in% aa.values()){
      first = stringr::str_sub(sub, 1, 1)
      sub = stringr::str_sub(sub, 2, -1)
    }

    if (stringr::str_sub(sub, -1, -1) %in% aa.values()){
      last = stringr::str_sub(sub, -1, -1)
      sub = stringr::str_sub(sub, 1, -2)
    }

    return(c(from = first, at = sub, to = last))
  }
  )

}


#'@export
subs.split.list <- function(subs.list){
  l = lapply(subs.list, subs.split)
  l = lapply(l, purrr::transpose)
  l = lapply(l, function(x)lapply(x,unlist))
  return(l)
}

