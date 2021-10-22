#' Add a column of annotations to a data frame
#'
#'
#'
#' @param x a data frame to which to add annotations
#' @param keyColumnName column name in x that matches the names in the annotation vector
#' @param annotationVector annotation lookup vector with names that correspond to keyColumnName
#' @param annotationName name for the new column that will be added to x containing annotations
#' @param parseProteinGroups flag indicating if the key column is protein groups that need to be parsed to match lookup vector (default = TRUE)
#' @param parsePTMs flag indicating if the key column contains PTM accessions that need to be parsed to match lookup vector (default = FALSE)
#'
#' @return an annotated data frame
#' @export
#'
#' @examples
annotateData <- function(x, keyColumnName, annotationVector, annotationName, parseProteinGroups = TRUE, parsePTMs = FALSE) {

  # If the key column is protein groups then split by semicolon character, otherwise don't split
  if (parseProteinGroups) {
    split_acc <- stringr::str_split(string = x[,keyColumnName], pattern = ";")
  } else {
    split_acc <- x[,keyColumnName]
  }

  # Loop through rows of split_acc
  for (i in 1:length(split_acc)) {

    # If accessions are PTM accessions, then parse out the protein accession from the PTM accession
    # E.g., parsing P01234_ph123 will give P01234 as the accession
    if (parsePTMs) {

      # ISSUE: if an accession contains non-alphanumeric characters (like underscores) this regex will not work correctly
      protein_accession <- stringr::str_match(split_acc[[i]], "([A-Z0-9]+)_[A-Za-z][A-Za-z]\\d+")[,2]
    } else {
      protein_accession <- split_acc[i]
    }

    # Create a lookup function to lookup annotations
    lookup_func <- function (x) unname(annotationVector[x])

    # Apply to accession(s) for this row and concatenate annotations together
    annotation <- paste(unlist(lapply(X = protein_accession, FUN = lookup_func)), collapse = ";")

    # Add to the annotation column in x
    x[i, annotationName] <- annotation
  }

  # Return the annotated data frame
  return(x)
}
