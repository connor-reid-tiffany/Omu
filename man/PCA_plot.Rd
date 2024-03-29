% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCA_plot.R
\name{PCA_plot}
\alias{PCA_plot}
\title{Create a PCA plot}
\usage{
PCA_plot(
  count_data,
  metadata,
  variable,
  color,
  response_variable = "Metabolite",
  label = FALSE,
  size = 2,
  ellipse = FALSE
)
}
\arguments{
\item{count_data}{Metabolomics count data}

\item{metadata}{Metabolomics metadata}

\item{variable}{The independent variable you wish to compare and contrast}

\item{color}{String of what you want to color by. Usually should be the same as variable.}

\item{response_variable}{String of the response_variable, usually should be "Metabolite"}

\item{label}{TRUE or FALSE, whether to add point labels or not}

\item{size}{An integer for point size.}

\item{ellipse}{TRUE or FALSE, whether to add confidence interval ellipses or not.}
}
\description{
Performs an ordination and outputs a PCA plot using a metabolomics
count data frame and metabolomics metadata
}
\examples{
PCA_plot(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
variable = "Treatment", color = "Treatment", response_variable = "Metabolite")
}
