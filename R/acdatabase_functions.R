#' Get gene attribute
#'
#' Find an attribute value for an antigen. By default this is set to traverse the inheritance list for 'HA' to find the attribute value but the gene can be set to 'NA' or 'backbone', or this behaviour turned off.
#'
#' @param ag environment: an antigen
#' @param attribute char: the name of the desired attiribute
#' @param inherit: bool: whether to inherit attributes from parents
#' @param inherit.gene char: which gene to inherit along
#'
#' @return attribute value
#' @export
#'
#' @examples
ag.attribute <- function(
  ag,
  attribute,
  inherit = TRUE,
  inherit.gene = "HA"
){

  # Assemble antigen query
  query <- paste0(
    "[['",
    paste0(
      strsplit(attribute, ".", TRUE)[[1]],
      collapse = "']][['"
    ),
    "']]"
  )

  # Get antigen inheritance
  if(inherit){
    ag_inheritance <- gene_inheritance(
      ag   = ag,
      gene = inherit.gene
    )
  } else {
    ag_inheritance <- list(ag)
  }

  # Get attribute
  for(ag in rev(ag_inheritance)){
    value <- eval(parse(text = paste0("ag", query)))
    if(!is.null(value)) break
  }

  # Return the value
  value

}


#' Get serum attribute
#'
#' Find an attribute value for an serum By default this is set to check parent sera if attribute not found at lower level.
#'
#' @param sr environment: a serum
#' @param attribute char: the name of the desired attiribute
#' @param inherit: bool: whether to inherit attributes from parents
#'
#' @return attribute value
#' @export
#'
#' @examples
sr.attribute <- function(sr, attribute, inherit = T){
  return(ag.attribute(sr, attribute, inherit, inherit.gene = 'serum'))
}
