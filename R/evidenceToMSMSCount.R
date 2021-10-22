#' Convert an evidence file to an msmscount.txt file
#'
#' @param evidence MaxQuant evidence data frame
#' @param keys keys data frame with columns for Raw.file, Condition, BioReplicate, Experiment, and IsotopeLabelType
#' @param removeCON flag indicating whether to remove CON__ entries (default = TRUE)
#' @param removeREV flag indicating whether to remove REV__ entries (default = TRUE)
#' @param removeEmpty flag indicating whether to remove entries with no accession (default = TRUE)
#'
#' @return a tab-delimited txt file os MSMS counts in wide format
#' @export
#'
#' @examples
#' evidence <- read.delim("evidence.txt", quote = "")
#' keys <- read.delim("keys.txt", quote = "")
#' msmscount_wide <- evidenceToMSMSCount(evidence = evidence, keys = keys)
evidenceToMSMSCount <- function (evidence, keys, removeCON = TRUE, removeREV = TRUE, removeEmpty = TRUE) {

  # If removeCON flas is true then remove CON__ entries
  if (removeCON) {
    evidence <- evidence[which(!(stringr::str_detect(string = evidence$Proteins, pattern = "CON__"))),]
  }

  # If removeREV flag is true them remove REV__ entries
  if (removeREV) {
    evidence <- evidence[which(!(stringr::str_detect(string = evidence$Proteins, pattern = "REV__"))),]
  }

  # If removeEmpty flag is true them remove entries with an empty accession
  if (removeEmpty) {
    evidence <- evidence[which(stringr::str_length(evidence$Proteins) > 0),]
  }

  # Merge evidence with keys data frame
  # Note: this merge will remove any raw files not present in the keys file
  # To do: warn user when raw file names don't match between keys and evidence
  evidence_keys <- merge(x = evidence, by.x = "Raw.file", y = keys, by.y = "Raw.file")

  # Cast to wide format to sum spectral counts by protein by Experiment
  msmscount_wide <- reshape2::dcast(data = evidence_keys, formula = Proteins ~ Experiment + Raw.file, value.var = "MS.MS.count", fun.aggregate = sum)

  # Return msmscount_wide data frame
  return(msmscount_wide)

}
