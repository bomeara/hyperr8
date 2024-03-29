% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperr8.R
\name{clean_input_data}
\alias{clean_input_data}
\title{Check and clean input data}
\usage{
clean_input_data(all_data)
}
\arguments{
\item{all_data}{A data frame with columns of time, rate, and perhaps more.}
}
\value{
A data frame with columns of time, rate, citation, numerator, denominator, and other columns in the input.
}
\description{
This will check and clean the input data. It expects column names of time, rate, and optionally citation; it can also include a numerator, denominator, and/or total_time columns. It will add missing columns.
}
