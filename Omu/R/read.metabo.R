#'read.metabo
#'wrapper for read.csv that appends the "cpd" class
#'@param filepath a file path in "" to your metabolomics count table
#'@example read.metabo("~/users/you/Desktop/your_metabolomics_data")
#'@export

read.metabo <- function(filepath){
  df <- read.csv(filepath)
  class(df) <- append(class(df),"cpd")
  return(df)
}
