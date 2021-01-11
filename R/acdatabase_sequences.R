
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




#'@export
nucs_to_protein = function(nucs){
  proteins = paste(seqinr::translate(stringr::str_split(nucs,'')[[1]], frame = 0, sens = 'F'), collapse = '')
  names(proteins) = names(nucs)
  return(proteins)
}

#'@export
get_protein_sequence <- function(fname, seqtype = 'unspecified', newlines = T){
  extension = rev(str_split(fname, fixed('.'))[[1]])[[1]]

  if (seqtype == 'unspecified'){

    if (extension == 'fas'){
      seqtype = 'protein'
    } else if (extension == 'fasta'|extension == 'seq') {
      seqtype = 'nucleotide'
    }
  }


  sequences = c()

  if (extension == 'fasta' | extension == 'fas'){
    seq_strings = readLines(fname)

    if (!newlines){
      for (i in 2*(1:(length(seq_strings)/2))-1 ){
        sequences[str_split(seq_strings[[i]],'>')[[1]][[2]]] = seq_strings[[i+1]]
      }
    }else{
      sequences[str_split(seq_strings[1],'>')[[1]][[2]]] = paste(seq_strings[-1], collapse = '')
    }

  } else if( extension == 'seq'){
    all_lines = readLines(fname)
    gene_seqs = Filter(function(line){stringr::str_length(line) > 5 & all(stringr::str_split(line, '')[[1]] %in% c('A', 'T' , 'C', 'G'))}, all_lines)
    if (length(gene_seqs) > 1) stop('Ambiguous .seq file')
    sequences = gene_seqs[[1]]
    names(sequences)[[1]] = 'ag'
  } else stop('Unrecognised file extension.')


  if (seqtype == 'nucleotide') sequences =  sapply(sequences, nucs_to_protein)


  return(sequences)

}

#'@export
trim_protein_sequence <- function(sequence, pattern){
  locs = stringr::str_locate_all(sequence, pattern)[[1]]
  if (dim(locs)[[1]] > 1){
    warning('Pattern detected more than once. Using first appearance')
    locs = locs[1,]
  }

  return(stringr::str_sub(sequence, locs[[1]]))

}






