#' Create files for SAINT scoring from an evidence file
#'
#' @param evidence MaxQuant evidence data frame
#' @param keys keys data frame where Condition indicates the baits (columns must include Raw.file, Condition, BioReplicate, Experiment, IsotopeLabelType, TestControl, and an optional SAINTGroups column)
#' @param db db object created by loadUniProtFASTA
#' @param prefix optional prefix for filenames created by this function (default = NULL)
#'
#' @return returns 1 if all files were written successfully
#' @export
#'
#' @examples db <- loadUniProtFASTA("SwissProt.Hsapiens.fasta")
#' evidence <- read.delim("evidence.txt")
#' keys <- read.delim("keys.txt")
#' evidenceToSAINT(evidence = evidence, keys = keys, db = db, prefix = "010121-JJ")
evidenceToSAINT <- function(evidence, keys, db, prefix = NULL) {

  # Get dbSequences from db object
  dbSequences <- db[["dbSequences"]]

  # Get names from db object
  dbNames <- as.data.frame(db[["getProteinNames"]])

  # Add a column for accessions
  dbNames$Accession <- names(db[["getProteinNames"]])

  # Change column names
  colnames(dbNames) <- c("Name", "Accession")


  # Need to cast evidence to wide format then melt to long format to fill in zeros across all samples
  evidence_wide <- reshape2::dcast(data = evidence, formula = Proteins ~ Raw.file, value.var = "MS.MS.count", fun.aggregate = sum)
  evidence_long <- reshape2::melt(data = evidence_wide, id.vars = "Proteins", variable.name = "Raw.file", value.name = "MS.MS.count")

  # Merge evidence_long with keys
  # To do: warn the user if raw files in evidence and keys don't match
  evidence_long_baits <- merge(x = evidence_long, by.x = "Raw.file", y = keys, by.y = "Raw.file")

  # Merge with sequence lengths
  # Note: this merge will ignore protein groups with more than one protein accession; SAINT wouldn't know how to handle anyway
  # To do: warn the user if proteins don't have sequence lengths in the dbSequence data frame
  evidence_long_baits_seqlength <- merge(x = evidence_long_baits, by.x = "Proteins", y = dbSequences[,c("Accession", "Length")], by.y = "Accession")

  # Merge with protein names df
  # Note: this merge will remove any accessions that do not have protein names!
  # To do: warn the user if accessions were removed
  evidence_long_baits_seqlength_names <- merge(x = evidence_long_baits_seqlength, by.x = "Proteins", y = dbNames, by.y = "Accession")



  # Check if there are groups in the keys file
  if ("SAINTGroups" %in% colnames(keys)) {

    # Get the groups
    groups <- unique(keys$SAINTGroups)

    # Loop through the groups and write SAINT files separately for each group
    for (i in 1:length(groups)) {

      # Get info for SAINT interactions file
      interactions <- evidence_long_baits_seqlength[which(evidence_long_baits_seqlength$SAINTGroups == groups[i]),c("Raw.file", "Condition", "Proteins", "MS.MS.count")]

      # Write out the interactions file
      if (!(is.null(prefix))) {
        write.table(x = interactions, file = paste(prefix, "-", groups[i], "-interactions.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      } else {
        write.table(x = interactions, file = paste(groups[i], "-interactions.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }

      # Get info for the SAINT preys file
      preys <- unique(evidence_long_baits_seqlength_names[which(evidence_long_baits_seqlength_names$SAINTGroups == groups[i]),c("Proteins", "Length", "Name")])

      # Write out the preys file
      if (!(is.null(prefix))) {
        write.table(x = preys, file = paste(prefix, "-", groups[i], "-preys.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      } else {
        write.table(x = preys, file = paste(groups[i], "-preys.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }

      # Get info for the baits file
      baits <- keys[which(keys$SAINTgroups == groups[i]), c("Raw.file", "Condition", "TestControl")]

      # Write out the baits file
      if (!(is.null(prefix))) {
        write.table(x = baits, file = paste(prefix, "-", groups[i], "-baits.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      } else {
        write.table(x = baits, file = paste(groups[i], "-baits.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      }

    }


  } else {

    # If there aren't groups then write SAINT files

    # Get info for SAINT interactions file
    interactions <- evidence_long_baits_seqlength[,c("Raw.file", "Condition", "Proteins", "MS.MS.count")]

    # Write to a file
    if (!(is.null(prefix))) {
      write.table(x = interactions, file = paste(prefix, "-interactions.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
      write.table(x = interactions, "interactions.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }

    # Get info for SAINT preys file
    preys <- unique(evidence_long_baits_seqlength_names[,c("Proteins", "Length", "Name")])

    # Write to a file
    if (!(is.null(prefix))) {
      write.table(x = preys, file = paste(prefix, "-preys.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
      write.table(x = preys, "preys.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }

    # Get info for the baits file
    baits <- keys[, c("Raw.file", "Condition", "TestControl")]

    # Write out the baits file
    if (!(is.null(prefix))) {
      write.table(x = baits, file = paste(prefix, "-baits.txt", collapse = "", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
      write.table(x = baits, file = "baits.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  }

  return(1)
}
