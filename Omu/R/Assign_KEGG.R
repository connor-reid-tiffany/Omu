#' Assign_KEGG
#' Assigns a KEGG value based off of Metabolite name
#'@param data A metabolomics data frame with a KEGG column 
#'@param KEGG_key A data frame from the Omu git repository

Assign_KEGG <- function(data, KEGG_key){
  data$KEGG2 <- KEGG_key$KEGG[match(data$Metabolite, KEGG_key$Metabolite)]
  data$KEGG[is.na(data$KEGG)] <- as.character(data$KEGG2[is.na(data$KEGG)])
  data = subset(data, select=-c(KEGG2))
  data = as.data.frame(data)
  
}
