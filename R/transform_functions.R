#' transform_samples
#' @description A functional to transform metabolomics data by samples
#' @param count_data Metabolomics data
#' @param func a function to transform samples by. can be an anonymous function
#' @examples
#' data_ln <- transform_samples(count_data = c57_nos2KO_mouse_countDF, log)
#' @export

transform_samples <- function(count_data, func){

  count_data[,sapply(count_data, is.numeric)] <- apply(count_data[,sapply(count_data, is.numeric)], 2, func)
  return(count_data)

}

#' transform_metabolites
#' @description A functional to transform metabolomics data by metabolites
#' @param count_data Metabolomics data
#' @param func a function to transform metabolites by. can be an anonymous function
#' @examples
#' data_pareto_scaled <- transform_samples(count_data = c57_nos2KO_mouse_countDF,
#' function(x) x/sqrt(sd(x)))
#' @export
transform_metabolites <- function(count_data,func){
  #set metabolite to rownames
  rownames(count_data) <- count_data$Metabolite
  #store non-numeric data in dataframe to remerge later
  char_data_cols <- sapply(count_data, function(x) !is.numeric(x))
  char_data <- count_data[,char_data_cols]
  #remove character data and transpose
  metabo_num <- count_data[,which(char_data_cols==FALSE)]
  metabo_num <- t(metabo_num)
  metabo_num <- apply(metabo_num, 2, func)
  #transpose, rejoin with character data by metabolite values
  metabo_num <- as.data.frame(t(metabo_num))
  metabo_num$Metabolite <- rownames(metabo_num)
  metabo_merge <- merge(char_data, metabo_num, by = "Metabolite")

  return(metabo_merge)

}
