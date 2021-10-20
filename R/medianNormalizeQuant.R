#' Median normalize an MSstats quant data frame
#'
#' @param quant MSstats-generated quant data frame
#'
#' @return a quant dataframe that is median normalized
#' @export
#'
#' @examples
#' evidence <- read.delim(file = "evidence.txt")
#' keys <- read.delim(file = "keys.txt")
#' proteinGroups <- read.delim(file = "proteinGroups.txt")
#' quant <- MaxQtoMSstatsFormat(evidence = evidence, annotation = keys, proteinGroups = proteinGroups)
#' quantNorm <- medianNormalizeQuant(quant)
medianNormalizeQuant <- function (quant) {

  # Filter out zero intensities for median intensity calculation
  quant_filtered <- quant[which(quant$Intensity > 0),]

  # Log 2 transform intensities
  quant_filtered$Log2Intensity <- log(quant_filtered$Intensity) / log(2)

  # Create a function to calculate median log2intensity by raw file
  median_func <- function (x) median(x = quant_filtered[which(quant_filtered$Run == x), "Log2Intensity"])

  # Get list of raw file names
  rawfiles <- as.character(unique(quant_filtered$Run))

  # Apply median function to list of raw file names
  medians <- unlist(lapply(X = rawfiles, FUN = median_func))

  # Combine raw files and medians in a data frame
  median_df <- cbind.data.frame(rawfiles, medians)

  # Get maximum median log2intensity for all raw files
  max_median <- max(median_df$medians)

  # Add a column to median_df and fill with the max median log2intensity
  median_df$Max <- max_median

  # Calculate the log2intensity needed to add to each to bring the median up
  median_df$Diff <- 2 ^ as.numeric(median_df$Max) - 2 ^ as.numeric(median_df$medians)

  for (i in 1:length(rawfiles)) {
    diff <- median_df[which(median_df$rawfiles == rawfiles[i]), "Diff"]
    quant[which((quant$Run == rawfiles[i]) & (quant$Intensity > 0)),"Intensity"] <- quant[which((quant$Run == rawfiles[i]) & (quant$Intensity > 0)),"Intensity"] + diff
  }

  # Return quant data frame
  return (quant)
}
