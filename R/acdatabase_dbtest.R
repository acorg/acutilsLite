# database and object checking
## ideally this or soemthing liek this should define the structure of the database?


#' Check antigen entry formatting
#'
#' Checks that a proposed antigen entry follows the formatting rules. See ?ag for more details on antigen objects.
#'
#' @details
#' The following are currently checked:
#' 1. class is c("acdatabase.ag", "acdatabase.entry", "environment")
#' 2. only permitted fields are entered
#' 3. all names (strain, long. type, subtype, lineage) are characters
#' 4. construct checkers are fulfilled for : id, parent_id, alterations, genes, isolation, passage, .parent.
#' @param ag environment: a proposed antigen entry
#' @param null_permitted bool: whether NULL should return T or F
#' @return bool
#' @export
#'
#' @examples
agdb.checkAG <- function(ag, null_permitted = F){

  if (is.null(ag)) return(null_permitted)

  # class:
  if (!(setequal( class(ag), c("acdatabase.ag", "acdatabase.entry", "environment") ))) stop('class must be  c("acdatabase.ag", "acdatabase.entry", "environment")')

  # fields
  permitted_fields = c('id', 'strain', 'long', 'aliases', 'wildtype', 'type', 'subtype', 'lineage', 'isolation', 'genes', 'parent_id', 'alterations', 'passage', 'comments', 'groups', 'meta', '.parent')
  if (!(all(names(ag) %in% permitted_fields))) stop('Unpermitted fields in ', names(ag))


  # strain, long, wildtype, type, subtype, lineage
  if (!(all(unlist(lapply(ag[c('strain', 'long', 'type', 'subtype', 'lineage')],
                          function(x){(typeof(x) == "character" | is.null(x))})))) ) {
    stop('all names must be char type')
  }

  # aliases
  if (!(all(unlist(lapply(ag$aliases, is.character))))) stop('all aliases must be char type')

  # id
  check_condition(agdb.checkid, ag$id)

  # isolation
  check_condition(agdb.checkisolation, ag$isolation)

  # genes
  check_condition(agdb.checkgenes, ag$genes)

  # parent_id
  check_condition(agdb.checkid, ag$parent_id, null_permitted = T)

  # alterations
  check_condition(agdb.checkalterations, ag$alterations)

  # passage
  check_condition(agdb.checkpassage, ag$passage)

  # .parent
  check_condition(agdb.checkAG, ag$.parent, null_permitted = T)

  return(TRUE)
}


#' Checks Antigen Database formatting
#'
#' Checks that an antigen database follows the formatting rules. For information on database structure, see ?agdb.
#'
#' @details
#' The following rules are currently checked:
#' 1. All entries are antigens (see ?agddb.checkAG)
#' @param agdb list
#'
#' @return bool
#' @export
#'
#' @examples
agdb.checkagdb <- function(agdb){
  return(all(unlist(lapply(agdb, agdb.checkAG))))
}



#' Check an id
#'
#' Checks if a proposed id follows the formatting rules.
#'
#' @details
#' The following rules are currently checked:
#' 1. id is be character type
#' 2. id has length 6
#' 3. all characters are capital letters or numerals 0:9
#'
#' @param id char: A proposed id
#' @param null_permitted bool: whether NULL should return T or F
#'
#' @return character
#' @export
#'
#' @examples
agdb.checkid <- function(id, null_permitted = F){
  if (is.null(id)) return(null_permitted)
  if (!is.character(id)) return(F)
  split_id = strsplit(id, split = '')[[1]]
  A = length(split_id) == 6
  B = all(unlist(lapply( split_id, function(x){x %in% c(LETTERS, 0:9)})))
  return(A & B)
}



#' Check an antigen isolation entry
#'
#' Checks if a proposed isolation entry follows the formatting rules.
#'
#' @details
#' The following rules are currently checked:
#' 1. isolation is a list
#' 2. only permitted fields are present
#'
#' @param isolation char: A proposed isolation entry
#' @param null_permitted bool: whether NULL should return T or F
#'
#' @return character
#' @export
#'
#' @examples
agdb.checkisolation <- function(isolation){
  if (is.null(isolation)) return (T)
  A = typeof(isolation) == 'list'
  B = all(names(isolation) %in% c('id', 'location', 'date', 'cell', 'continent'))

  return(A & B)
}


#' Check an antigen genes entry
#'
#' Checks if a proposed genes entry follows the formatting rules.
#'
#' @details
#' The following rules are currently checked:
#' 1. genes is a list
#' 2. only permitted fields are present for each gene in genes
#' 2. all genes are HA or NA
#'
#' @param genes char: A proposed genes entry
#' @param null_permitted bool: whether NULL should return T or F
#'
#' @return character
#' @export
#'
#' @examples
agdb.checkgenes <- function(genes){
  if (is.null(genes)) return(T)

  if (!(typeof(genes) == 'list')) return (F)


  A = all(unlist(lapply(genes, function(gene){
    all(names(gene) %in% c('gene', 'sequence', 'clade', 'cluster'))
  })))

  B = all(unlist(lapply(genes, function(gene){
    gene$gene %in% c('HA', 'NA')
  })))


  return(A & B)
}


#' Check an antigen alterations entry
#'
#' Checks if a proposed alterations entry follows the formatting rules.
#'
#' @details
#' The following rules are currently checked:
#' 1. alterations is a list
#' 2. only permitted fields are present for each alteration in alterations
#' 3. all substitutions are valid
#' 4. all parent_id's are valid
#' 5. all .parent's are valid antigens
#'
#' @param alterations char: A proposed alterations entry
#' @param null_permitted bool: whether NULL should return T or F
#'
#' @return character
#' @export
#'
#' @examples
agdb.checkalterations <- function(alterations){
  if (is.null(alterations)) return(T)

  A = (typeof(alterations) == 'list')

  B = all(unlist(lapply(alterations, function(alteration){
    A = (all(names(alteration) %in% c('gene', 'substitutions', 'parent_id')));
    B = (alteration$gene %in% c('HA', 'NA'));
    C = all(unlist(lapply(alteration$substitutions, is.substitution)))
    D = is.null(alteration$parent_id) || agdb.checkid(alteration$parent_id)
    E = is.null(alteration$.parent) || agdb.checkAG(alteration$.parent)
  })))

  return(all(A, B))
}


#' Check a substitutions entry
#'
#' Checks if a proposed substitution follows the correct formatting
#'
#' @details
#' The following rules are currently checked:
#' 1. substitution is of character type
#' 2. starts and ends with valid single letter amino acid code
#' 3. central characters form a numerb in 1:484
#'
#' @param subs char: A proposed substitution entry
#'
#' @return character
#' @export
#'
#' @examples
is.substitution <- function(subs){
  if(is.null(subs)) return(T)
  A = stringr::str_sub(subs,1,1) %in% aavalues()
  B = stringr::str_sub(subs,-1,-1) %in% aavalues()
  C = as.numeric(stringr::str_sub(subs,2,-2), ''[[1]]) %in% 1:484
  all(A, B, C)
}



#' Check a passage entry
#'
#' Checks if a proposed passage entry follows the correct formatting
#'
#' @details
#' The following rules are currently checked:
#' 1. passage is of list type
#' 2. only permitted names are present
#'
#' @param passage char: A proposed passage entry
#'
#' @return character
#' @export
#'
#' @examples
agdb.checkpassage <- function(passage){
  if(is.null(passage)) return(T)

  A = all(names(passage) %in% c('history', 'egg', 'cell', 'details', 'comments'))

  return(all(A))
}




#' Check if a condition is fulfilled
#'
#' Checks if checker_fn returns TRUE when passed arguments '...'. Otherwise throws error.
#'
#' @param checker_fn function: a boolean valued function checking if a condition is fulfilled by the arguments ...
#' @param ... arguments to checker_fn
#'
#' @return bool
#' @export
check_condition <- function(checker_fn, ...){
  fn_char = as.character(enexpr(checker_fn))
  if (!((checker_fn(...)))) stop('Error in ', fn_char)
  return(T)
}
