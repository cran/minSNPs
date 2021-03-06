% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{flag_position}
\alias{flag_position}
\title{\code{flag_position}}
\usage{
flag_position(
  pro_seqc,
  dash_ignore = TRUE,
  accepted_char = c("A", "C", "T", "G"),
  ignore_case = TRUE,
  bp = SerialParam()
)
}
\arguments{
\item{pro_seqc}{Sequences after processed, i.e. all with the same length}

\item{dash_ignore}{whether to treat '-' as another type}

\item{accepted_char}{character to accept, default to c("A", "C", "T", "G")}

\item{ignore_case}{whether to be case insensitive, default to TRUE}
}
\value{
Will return a list of positions that need to be ignored.
}
\description{
\code{flag_position} is used to find out positions that will be ignored in
calculation (either not A,C,G,T or '-'), can be case sensitive
or insensitive.
}
\keyword{internal}
