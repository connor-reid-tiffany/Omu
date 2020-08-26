#' Get metadata from KEGG API
#'
#' @description Internal function for KEGG_Gather
#' @param count_data The metabolomics count data
#' @param column The name of the KEGG identifier being sent to the KEGG API
#' @param req The character vector with variables to pull out of the metadata
#' @importFrom KEGGREST keggGet
#' @importFrom plyr llply
#' @export


make_omelette <- function(count_data, column, req){

  #Create input of nested lists of length 10 and feed to KEGG API
  input <- data.frame(i = count_data[, column])
  input = as.vector(input[rowSums(is.na(input)) != ncol(input),])
  input_split  <- split(input,  ceiling(seq_along(input)/10))
  input_split = llply(input_split, as.list)

  #Send nested list of identifier values to KEGG API
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

#' Clean up orthology metadata
#'
#' @description Internal function for KEGG_Gather.rxn method
#' KEGG_Gather.rxn requires dispatch on multiple elements, so
#' There was no way to incorporate as a method
#' @param count_data Metabolomics count data
#' @param matrix  the matrix of KEGG metadata
#' @importFrom reshape2 colsplit
#' @importFrom dplyr left_join
#' @export

plate_omelette_rxnko<- function(count_data, matrix){

  #Unlist orthology data, extract KO numbers, and element numbers to join with Rxn numbers
  ko_df = as.data.frame(unlist(matrix[,2], recursive = F))
  ko_df <- cbind(rownames(ko_df), data.frame(ko_df, row.names=NULL))
  colnames(ko_df)[1] <- "KO_Number"
  ko_df = with(ko_df, cbind(ko_df[,2],
                            colsplit(ko_df$KO_Number, pattern = "\\.",
                                     names = c('Element_number', 'KO_Number'))))
  colnames(ko_df)[1] <- "KO"

  #Do the above but for Rxn numbers this time
  rxn_df = as.data.frame(unlist(matrix[,1], recursive = F))
  rxn_df <- cbind(rownames(rxn_df), data.frame(rxn_df, row.names=NULL))
  colnames(rxn_df)[1] <- "Element_number"
  rxn_df = with(rxn_df, cbind(rxn_df[,2],
                              colsplit(rxn_df$Element_number, pattern = "\\.",
                                       names = c('Element_number', 'remove'))))
  colnames(rxn_df)[1] <- "Rxn"
  rxn_df = rxn_df[,1:2]

  #Join original dataframe to dataframe with KO by Rxn number
  ko_with_rxn <- left_join(ko_df, rxn_df, by = "Element_number")
  count_data$KO <- ko_with_rxn$KO[match(count_data$Rxn, ko_with_rxn$Rxn)]
  count_data$Element_number <- ko_with_rxn$Element_number[match(count_data$Rxn, ko_with_rxn$Rxn)]
  count_data$KO_Number <- ko_with_rxn$KO_Number[match(count_data$Rxn, ko_with_rxn$Rxn)]
  #Drop the element number column as it is no longer needed
  count_data = count_data[ , !(names(count_data) %in% "Element_number")]

  return(count_data)
}
