% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evidenceToMSMSCount.R
\name{evidenceToMSMSCount}
\alias{evidenceToMSMSCount}
\title{Convert an evidence file to an msmscount.txt file}
\usage{
evidenceToMSMSCount(
  evidence,
  keys,
  removeCON = TRUE,
  removeREV = TRUE,
  removeEmpty = TRUE
)
}
\arguments{
\item{evidence}{MaxQuant evidence data frame}

\item{keys}{keys data frame with columns for Raw.file, Condition, BioReplicate, Experiment, and IsotopeLabelType}

\item{removeCON}{flag indicating whether to remove CON__ entries (default = TRUE)}

\item{removeREV}{flag indicating whether to remove REV__ entries (default = TRUE)}

\item{removeEmpty}{flag indicating whether to remove entries with no accession (default = TRUE)}
}
\value{
a tab-delimited txt file os MSMS counts in wide format
}
\description{
Convert an evidence file to an msmscount.txt file
}
\examples{
evidence <- read.delim("evidence.txt", quote = "")
keys <- read.delim("keys.txt", quote = "")
msmscount_wide <- evidenceToMSMSCount(evidence = evidence, keys = keys)
}
