% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{check_multistate}
\alias{check_multistate}
\title{\code{check_multistate}}
\usage{
check_multistate(position, sequences)
}
\arguments{
\item{position}{position to check}

\item{sequences}{sequences from group of interest}
}
\value{
return `TRUE` if the position contains multistate otherwise `FALSE`
}
\description{
\code{check_multistate} is used to remove positions where there are
more than 1 state within the group of interest.
}
