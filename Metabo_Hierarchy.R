#' Metabo_Path
#' Annotates a metabolomics data frame with metabolite hiearchy data
#'@param data A metabolomics data frame with a KEGG column 
#'@param KEGG_key A data frame from the Omu git page with metabolite hiearchy data in the form of Class, Subclass_1, Subclass_2, Subclass_3, and Subclass_4 columns
#'@example Metabo_Path(data = data, KEGG_Key = KEGG_Key)
Metabo_Path <- function(data, KEGG_key){
  data$Class = KEGG_key$Class[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_1 = KEGG_key$Subclass_1[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_2 = KEGG_key$Subclass_2[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_3 = KEGG_key$Subclass_3[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_4 = KEGG_key$Subclass_4[match(data$KEGG, KEGG_key$KEGG)]
  rownames(data) <- data[,1]
  data = data[, !names(data) %in% c('Metabolite')]
  
}
  