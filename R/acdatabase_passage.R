

#'@export
merge_on_passage <- function(titertable, agdb){

  # Fetch records
  ag_records  <- acdb.getIDs(agdb, rownames(titertable))

  # Get passage parent IDs but indicate if it included egg passage
  parent_ids  <- sapply(ag_records, function(ag){ rev(ag.passageInheritance(ag))[[1]]$id })
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
