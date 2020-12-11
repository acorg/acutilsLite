
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

#'@export
geneseq.diff <- function(seq1, seq2){
  min_len = min(stringr::str_length(seq1), stringr::str_length(seq2))
  if (stringr::str_length(seq1) != stringr::str_length(seq2) ) warning('Sequence lengths differ: comapring only first ', as.character(min_len), ' loci' )

  original_split <- str_split(seq1, '')[[1]][1:min_len]
  derived_split <- str_split(seq2, '')[[1]][1:min_len]


  return(unlist(mapply(function(s1, s2, n){if (s1!=s2){return(paste0(s1,n,s2))}}, original_split, derived_split, 1:length(original_split) )))
}
