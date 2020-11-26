

ag.passageParent <- function(ags){
  lapply(ags, function(ag) {
    passage_inheritance(ag)[[1]]
  })
}


passage_inheritance <- function(ag){

  inheritance <- list()
  while(1==1){
    inheritance <- c(inheritance, list(ag))
    if(is.null(ag$parent_id) || !is.null(ag[["alterations"]])) break
    ag <- ag$.parent
  }
  rev(inheritance)

}


passage_includes_egg <- function(ag){

  passage_history <- unlist(inherited_passage_history(ag))
  sum(grepl("E", passage_history, fixed = T)) > 0

}


inherited_passage_history <- function(ag){

  ag_passage_inheritance <- passage_inheritance(ag)
  lapply(ag_passage_inheritance, function(x){ x[["passage"]][["history"]] })

}


merge_on_passage <- function(titertable, agdb){

  # Fetch records
  ag_records  <- acdb.getIDs(agdb, rownames(titertable))

  # Get passage parent IDs but indicate if it included egg passage
  parent_ids  <- sapply(ag_records, function(ag){ passage_inheritance(ag)[[1]]$id })
  passage_egg <- sapply(ag_records, passage_includes_egg)
  parent_ids[passage_egg] <- paste(parent_ids[passage_egg], "EGG")

  # Determine duplicate rows
  duplicate_rows <- which(duplicated(parent_ids))

  # Merge together the passage duplicate rows
  for(x in duplicate_rows){

    matched_x <- match(parent_ids[x], parent_ids)
    titers1 <- titertable[matched_x,]
    titers2 <- titertable[x,]
    if(sum(titers1 != "*" & titers2 != "*") > 0){
      stop("Would have to merge titers")
    }
    na_titers1 <- titers1 == "*"
    titertable[matched_x,na_titers1] <- titers2[na_titers1]

  }

  # Remove duplicate rows
  titertable[-duplicate_rows,,drop=F]

}
