#' read.metabo
#' wrapper for read.csv that appends the "cpd" class
#' @param filepath a file path in "" to your metabolomics count table
#' @importFrom utils read.csv

#' @export

read.metabo <- function(filepath){
  df <- read.csv(filepath)
  class(df) <- append(class(df),"cpd")
  return(df)
}
