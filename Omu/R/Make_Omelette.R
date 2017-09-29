
Make_Omelette <- function(countDF, column, Reqs){

  #Create input of nested lists of length 10 and feed to KEGG API
  Input <- data.frame(i = test_metabo[, column])
  Input = as.vector(Input[rowSums(is.na(Input)) != ncol(Input),])
  Input_Split  <- split(Input,  ceiling(seq_along(Input)/10))
  Input_Split = llply(Input_Split, as.list)

  #Send nested list of Cpd numbers to KEGG API
  Output <- llply(Input_Split, function(x)keggGet(x))
  Unlist_Output <- unlist(Output, recursive = F)

  #Pull out Rxn and Cpd numbers
  Extract_reqs <- lapply(Unlist_Output, '[', Reqs)

  #Convert data into list of character vectors that can be added to an R data frame object
  #Note that character vectors
  n_obs <- sapply(Extract_reqs, length)
  seq_max <- seq_len(max(n_obs))
  Matrix <- t(sapply(Extract_reqs, '[', i = seq_max)

  return(Matrix)

}
