#' Import a metabolomics count data frame
#'
#' @description Wrapper for read.csv that appends the "cpd" class and sets blank cells to NA. Used to
#' import metabolomics count data into R.
#' @param filepath a file path to your metabolomics count data
#' @importFrom utils read.csv
#' @examples
#' filepath_to_yourdata = paste0(system.file(package = "omu"), "/extdata/read_metabo_test.csv")
#' count_data <- read_metabo(filepath_to_yourdata)
#' @export


read_metabo <- function(filepath){
  df <- read.csv(filepath, na.strings = "")
  class(df) <- append(class(df),"cpd")
  return(df)
}
