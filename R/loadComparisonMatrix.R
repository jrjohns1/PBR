#' Title
#'
#' @param filename comparisons file
#'
#' @return a matrix of comparisons compatible with MSstats
#' @export
#'
#' @examples
#' comparisonMatrix <- loadComparisonMatrix("comparisons.txt")
loadComparisonMatrix <- function (filename) {

  # Read in comparisons file
  comparisons <- read.delim(file = filename, quote = "", stringsAsFactors = FALSE)

  # Grab condition columns from comparisons
  comparisonmatrix <- comparisons[,2:ncol(comparisons)]

  # Reorder conditions in alphabetical order (MSstast will fail if you don't do this)
  comparisonmatrix <- as.matrix(comparisonmatrix[,order(names(comparisonmatrix))])

  # Add names of comparisons as row names
  row.names(comparisonmatrix) <- comparisons[,1]

  # Return comparisonmatrix
  return(comparisonmatrix)
}
