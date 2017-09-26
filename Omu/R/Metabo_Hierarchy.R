#' Metabo_Hierarchy
#' Annotates a metabolomics data frame with metabolite hiearchy data
#'@param data A metabolomics data frame with a KEGG column 
#'@param KEGG_key A data frame from the Omu git page with metabolite hiearchy data in the form of Class, Subclass_1, Subclass_2, Subclass_3, and Subclass_4 columns
#'@param keep_unknowns FALSE removes unknown metabolites, use TRUE to keep
#'@example Metabo_Path(data = data, KEGG_Key = KEGG_Key)
Metabo_Hierarchy <- function(data, Metabolite_Hierarchy, keep_unknowns){
  if (keep_unknowns ==FALSE){
    data = inner_join(data, Metabolite_Hierarchy, by = "KEGG")
    data <- data[1:(length(data)-2)]
    colnames(data)[1] <- "Metabolite"
    data = distinct(data, Metabolite,.keep_all = TRUE)
    return(data)
    }
    else if (keep_unknowns == TRUE) {
      data = left_join(data, Metabolite_Hierarchy, by = "KEGG")
      data <- data[1:(length(data)-2)]
      colnames(data)[1] <- "Metabolite"
      data = distinct(data, Metabolite,.keep_all = TRUE)
      return(data)
    }
  
}
  
