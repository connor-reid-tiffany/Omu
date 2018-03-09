#' assign_hierarchy
#' 
#' This function assigns hierarchy metadata to a metabolomics count matrix with KEGG identifier numbers.
#' It can assign KEGG compound hierarchy, orthology hierarchy, or organism hierarchy
#' @param data a metabolomics count matrix with either a KEGG compound, orthology, or gene identifier
#' @param file_path a file path to either the metabolite hierarchy, orthology hierarchy, or organism hierarchy
#' @param keep_unknowns a boolean of either TRUE or FALSE. TRUE keeps unannotated compounds, FALSE prunes them
#' @param identifier a string that is either "KEGG" for metabolite classification, "KO_Number" 
#' for orthology classification, or "Org" for organism classification
#' @example assign_hierarchy(data = yourmetabolomicsdata, file_path = "~/Desktop/Metabolite_Hierarchy.csv", keep_unknowns = TRUE, identifier = "KEGG")
#' @export 





assign_hierarchy <- function(data, file_path, keep_unknowns, identifier){
  hierarchy <- read.csv(file_path)

if (identifier == "KEGG"){
  if (keep_unknowns ==FALSE){
    data <- inner_join(data, hierarchy, by = identifier)
    data <- distinct(data, Metabolite,.keep_all = TRUE)
    return(data)
    }
    else if (keep_unknowns == TRUE) {
      data <- left_join(data, hierarchy, by = identifier)
      data <- distinct(data, Metabolite,.keep_all = TRUE)
      return(data)
    }
  } else if (identifier == "KO_Number"){

    data$KO_Class <- hierarchy$KO_Class[match(data$KO_Number,
      hierarchy$KO_Number)]
    data$KO_Subclass1 <- hierarchy$KO_Subclass1[match(data$KO_Number,
      hierarchy$KO_Number)]
    data$KO_Subclass2 <- hierarchy$KO_Subclass2[match(data$KO_Number,
        hierarchy$KO_Number)]

    return(data)
  } else if (identifier == "Org"){
    data <- inner_join(data, hierarchy, by = identifier)
    return(data)
    

    return(data)
  }
}
