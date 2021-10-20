


#' Load a fasta-formatted UniProt protein sequence database
#'
#' @param dbFilenames list of database file names
#' @param getSequences returns a data frame named dbSequences of protein sequences and protein lengths (default = FALSE)
#'
#' @return List of lookup vectors of protein names (getProteinNames), protein descriptions (getProteinDescriptions) and/or data frame of sequences (dbSequences)
#' @export
#'
#' @examples
#' db <- loadUniProtFASTA(dbFilenames = c("SwissProt.Hsapiens.fasta", "Uniprot.Ggallus.fasta"))
#' getProteinNames <- db[["getProteinNames"]]
#' getProteinDescriptions <- db[["getProteinDescriptions"]]
#' proteinName <- getProteinNames[["P11388"]]
#' proteinDescription <- getProteinDescriptions[["P11388"]]
#'
#' db <- loadUniProtFASTA("SwissProt.Human.fasta", getSequences = TRUE)
#' dbSequences <- db[["dbSequences"]]
#' sequence <- dbSequences[which(dbSequences$Accession == "P11388), "Sequence"]
#' sequenceLength <- dbSequences[which(dbSequences$Accession == "P11388), "Length"]
loadUniProtFASTA <- function(dbFilenames, getSequences = FALSE) {

  # Get the number of databases
  numdb <- length(dbFilenames)

  # Loop through each filename provided as input
  for (i in 1:numdb) {

    # If it's the first database then read lines into dblines
    if (i == 1) {
      dblines <- readLines(con = dbFilenames[i])

      # If it's not the first db then read lines into tmp
      # then concatenate tmp onto the end of dblines
    } else {
      tmp <- readLines(con = dbFilenames[i])
      dblines <- c(dblines, tmp)
    }
  }

  ########## Read in protein names and descriptions ##########


  # Grab header lines starting with > character
  headerlines <- dblines[grep(">", dblines)]

  # Extract protein accession and protein name from header lines
  protein.accession <- stringr::str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,2]
  protein.description <- stringr::str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+.*)")[,3]
  protein.name <- stringr::str_match(headerlines, ">\\S\\S\\|(\\S+)\\|(\\S+)_.*")[,3]

  # Create named vectors to lookup protein names and descriptions by accession
  getProteinDescriptions <- protein.description
  names(getProteinDescriptions) <- protein.accession
  getProteinNames <- protein.name
  names(getProteinNames) <- protein.accession




  ########## Read in sequences ##########

  # If getSequences flag is true
  if (getSequences == TRUE) {

    # Set a flag for when to read sequences to zero
    grabsequences <- 0

    # Declare lists to store sequences and accessions
    sequencelist <- list()
    accessionlist <- list()

    # Loop through dbs line by line
    for (i in 1:length(dblines)) {

      # If the line contains a > then this is a description line
      if (grepl(pattern = ">", x = dblines[i])) {

        if (i > 1) {
          # If this is not the first line then append the sequences and accession to their respective lists
          sequencelist <- append(sequencelist, sequence)
          accessionlist <- append(accessionlist, accession)
        }

        # Extract accession based on the pattern: ">sp|P12345|ASDF_HUMAN blah blah", where P12345 is accession
        accession <- stringr::str_match(string = dblines[i], pattern = ">\\S\\S\\|([^\\|]*)\\|.*")[2]

        # Next line is the first line of sequence so set grab sequence flag to 1
        grabsequence <- 1

      } else {
        if (grabsequence == 1) {
          # If this is not a header line and grabsequence flag is set then read line into sequence
          sequence <- dblines[i]

          # And set the grabsequence flag back to zero
          grabsequence <- 0
        } else {
          # Otherwise append the sequence to the existing sequqence
          sequence <- paste(sequence, dblines[i], sep = '', collapse = '')
        }
      }
    }

    # Grab the last ones in the file
    sequencelist <- append(sequencelist, sequence)
    accessionlist <- append(accessionlist, accession)

    # Convert accession and sequence lists to data frame
    accessiondf <- as.data.frame(unlist(accessionlist))
    sequencedf <- as.data.frame(unlist(sequencelist))

    # Bind them together
    accseqdf <- cbind(accessiondf, sequencedf)

    # Rename columns
    colnames(accseqdf) <- c("Accession", "Sequence")

    # Add a column for sequence length
    accseqdf$Length <- stringr::str_length(accseqdf$Sequence)

  }

  if (getSequences == TRUE) {
    return(list("getProteinNames" = getProteinNames, "getProteinDescriptions" = getProteinDescriptions, "dbSequences" = accseqdf))
  } else {
    return(list("getProteinNames" = getProteinNames, "getProteinDescriptions" = getProteinDescriptions))
  }


}
