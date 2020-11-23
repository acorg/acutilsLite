
# Fetch attributes for an antigen
#' @export
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
