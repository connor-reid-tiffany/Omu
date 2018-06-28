#' assign_hierarchy
#'
#' This function assigns hierarchy metadata to a metabolomics count matrix with KEGG identifier numbers.
#' It can assign KEGG compound hierarchy, orthology hierarchy, or organism hierarchy
#' @param count_data a metabolomics count data frame with either a KEGG compound, orthology,
#' or a gene identifier column
#' @param keep_unknowns a boolean of either TRUE or FALSE. TRUE keeps unannotated compounds,
#' FALSE removes them
#' @param identifier a string that is either "KEGG" for metabolite, "KO_Number" for orthology,
#' "Prokaryote" for organism, or "Eukaryote" for organism
#' @importFrom dplyr inner_join
#' @importFrom dplyr left_join
#' @importFrom dplyr distinct
#' @examples
#' assign_hierarchy(count_data = c57_nos2KO_mouse_countDF, keep_unknowns = TRUE, identifier = "KEGG")
#' @export

assign_hierarchy <- function(count_data, keep_unknowns, identifier){
Metabolite <- NULL
identifier = match.arg(arg = identifier, choices = c("KEGG", "KO_Number", "Prokaryote", "Eukaryote"))

if (identifier == "KEGG"){
  if (keep_unknowns ==FALSE){
    count_data <- inner_join(count_data, Metabolite_Hierarchy_Table, by = identifier)
    count_data <- distinct(count_data, Metabolite,.keep_all = TRUE)
    class(count_data) = append(class(count_data), "cpd")
    return(count_data)
    }
    else if (keep_unknowns == TRUE) {
      count_data <- left_join(count_data, Metabolite_Hierarchy_Table, by = identifier)
      count_data <- distinct(count_data, Metabolite,.keep_all = TRUE)
      class(count_data) = append(class(count_data), "cpd")
      return(count_data)
    }
  } else if (identifier == "KO_Number"){

    count_data$KO_Class <- Orthology_Hierarchy_Table$KO_Class[match(count_data$KO_Number,
      Orthology_Hierarchy_Table$KO_Number)]
    count_data$KO_Subclass_1 <- Orthology_Hierarchy_Table$KO_Subclass_1[match(count_data$KO_Number,
      Orthology_Hierarchy_Table$KO_Number)]
    count_data$KO_Subclass_2 <- Orthology_Hierarchy_Table$KO_Subclass_2[match(count_data$KO_Number,
        Orthology_Hierarchy_Table$KO_Number)]

    return(count_data)
  } else if (identifier == "Prokaryote"){

    count_data <- inner_join(count_data, Prokaryote_Hierarchy_Table, by = 'Org')



    return(count_data)
  } else if (identifier == "Eukaryote"){

    count_data <- inner_join(count_data, Eukaryote_Hierarchy_Table, by = 'Org')

    return(count_data)
  }

}
