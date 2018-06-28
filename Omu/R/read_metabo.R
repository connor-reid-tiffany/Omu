#' read_metabo
#' wrapper for read.csv that appends the "cpd" class and sets blank cells to NA
#' @param filepath a file path in "" to your metabolomics count table
#' @importFrom utils read.csv
#' @examples
#' #See vignette for information on read_metabo()
#' @export


read_metabo <- function(filepath){
  df <- read.csv(filepath, na.strings = "")
  class(df) <- append(class(df),"cpd")
  return(df)
}
