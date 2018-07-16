#' Import a metabolomics count data frame
#'
#' @description Wrapper for read.csv that appends the "cpd" class and sets blank cells to NA. Used to
#' import metabolomics count data into R.
#' @param filepath a file path in "" to your metabolomics count table
#' @importFrom utils read.csv
#' @examples
#' \dontrun{count_data <- read_metabo("path/to/your/data")}
#' @export


read_metabo <- function(filepath){
  df <- read.csv(filepath, na.strings = "")
  class(df) <- append(class(df),"cpd")
  return(df)
}
