#' assign_hierarchy
#'
#' This function assigns hierarchy metadata to a metabolomics count matrix with KEGG identifier numbers. It can assign KEGG compound hierarchy, orthology hierarchy, or organism hierarchy
#' @param data a metabolomics count matrix with either a KEGG compound, orthology, or gene identifier
#' @param keep_unknowns a boolean of either TRUE or FALSE. TRUE keeps unannotated compounds, FALSE prunes them
#' @param identifier a string that is either "KEGG" for metabolite, "KO_Number" for orthology,"Prokaryote" for organism, or "Eukaryote" for organism
#' @example assign_hierarchy(data = yourmetabolomicsdata, file_path = "~/Desktop/Metabolite_Hierarchy.csv", keep_unknowns = TRUE, identifier = "KEGG")
#' @export

assign_hierarchy <- function(data, keep_unknowns, identifier){

if (identifier == "KEGG"){
  if (keep_unknowns ==FALSE){
    data <- inner_join(data, Metabolite_Hierarchy_Table, by = identifier)
    data <- distinct(data, Metabolite,.keep_all = TRUE)
    return(data)
    }
    else if (keep_unknowns == TRUE) {
      data <- left_join(data, Metabolite_Hierarchy_Table, by = identifier)
      data <- distinct(data, Metabolite,.keep_all = TRUE)
      return(data)
    }
  } else if (identifier == "KO_Number"){

    data$KO_Class <- Orthology_Hierarchy_Table$KO_Class[match(data$KO_Number,
      Orthology_Hierarchy_Table$KO_Number)]
    data$KO_Subclass_1 <- Orthology_Hierarchy_Table$KO_Subclass_1[match(data$KO_Number,
      Orthology_Hierarchy_Table$KO_Number)]
    data$KO_Subclass_2 <- Orthology_Hierarchy_Table$KO_Subclass_2[match(data$KO_Number,
        Orthology_Hierarchy_Table$KO_Number)]

    return(data)
  } else if (identifier == "Prokaryote"){
    data <- inner_join(data, Prokaryote_Hierarchy_Table, by = identifier)



    return(data)
  } else if (identifier == "Eukaryote"){
    data <- inner_join(data, Eukaryote_Hierarchy_Table, by = identifier)

    return(data)
  }
}
