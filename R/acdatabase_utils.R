#' @export
gene_inheritance <- function(
  ag,
  gene,
  inheritance = list()){

  # Check input
  if(!inherits(ag, "acdatabase.entry")){
    stop(sprintf("Class of input must inherit from 'acdatabase.entry' but was '%s'", paste(class(ag), collapse = ", ")))
  }

  # Get the ha alterations
  ha_alteration <- match(gene, sapply(ag[["alterations"]], function(x) x$gene))

  # Add any substitutions listed
  inheritance <- c(list(ag), inheritance)

  # Look for inheritance
  if(!is.na(ha_alteration) && !is.null(ag[["alterations"]][[ha_alteration]][["parent_id"]])){
    gene_inheritance(
      ag            = ag[["alterations"]][[ha_alteration]][[".parent"]],
      gene          = gene,
      inheritance   = inheritance
    )
  } else if(!is.null(ag$parent_id)){
    gene_inheritance(
      ag            = parent.env(ag),
      gene          = gene,
      inheritance   = inheritance
    )
  } else{
    inheritance
  }

}
