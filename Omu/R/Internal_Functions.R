#' Make_Omelette
#' Internal function for KEGG_Gather
#'@export


make_omelette <- function(countDF, column, req){

  #Create input of nested lists of length 10 and feed to KEGG API
  input <- data.frame(i = countDF[, column])
  input = as.vector(input[rowSums(is.na(input)) != ncol(input),])
  input_split  <- split(input,  ceiling(seq_along(input)/10))
  input_split = llply(input_split, as.list)

  #Send nested list of Cpd numbers to KEGG API
  output <- llply(input_split, function(x)keggGet(x), .progress = "text")
  unlist_output <- unlist(output, recursive = F)

  #Pull out Rxn and Cpd numbers
  extract_reqs <- lapply(unlist_output, '[', req)

  #Convert data into list of character vectors that can be added to
  #an R data frame object

  n_obs <- sapply(extract_reqs, length)
  seq_max <- seq_len(max(n_obs))
  matrix <- t(sapply(extract_reqs, '[', i = seq_max))

  return(matrix)

}

#' Plate_Omelette_rxnKO
#' Internal function for KEGG_Gather.rxn method
#' KEGG_Gather.rxn requires dispatch on multiple elements, so
#' There was no way to incorporate as a method
#' @export

plate_omelette_rxnko<- function(countDF, matrix){

ko_df = as.data.frame(unlist(matrix[,2], recursive = F))
ko_df = rownames_to_column(df, "KO_Number")
ko_df = with(ko_df, cbind(DF[,2],
colsplit(ko_df$KO_Number, pattern = "\\.",
names = c('Element_number', 'KO_Number'))))
colnames(ko_df)[1] <- "KO"

#Do the above but for Rxn numbers this time
rxn_df = as.data.frame(unlist(matrix[,1], recursive = F))
rxn_df = rownames_to_column(rxn_df, "Element_number")
rxn_df = with(rxn_df, cbind(rxn_df[,2],
  colsplit(rxn_df$Element_number, pattern = "\\.",
  names = c('Element_number', 'remove'))))
colnames(rxn_df)[1] <- "Rxn"
rxn_df = rxn_df[,1:2]

#Join origianl dataframe to dataframe with KO by Rxn number
ko_with_rxn <- left_join(ko_df, rxn_df, by = "Element_number")
countDF= left_join(ko_with_rxn, countDF, by = "Rxn")

#Drop the element number column as it is no longer needed
countDF = countDF[ , !(names(countDF) %in% "Element_number")]

return(countDF)
}
