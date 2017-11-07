#' assign_hierarchy
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

    data$Kingdom <- hierarchy$Kingdom[match(data$Org, hierarchy$Org)]
    data$Supergroup <- hierarchy$Supergroup[match(data$Org, hierarchy$Org)]
    data$Supergroup <- hierarchy$Supergroup[match(data$Org, hierarchy$Org)]
    data$Genus <- hierarchy$Genus[match(data$Org, hierarchy$Org)]
    data$Species <- hierarchy$Species[match(data$Org, hierarchy$Org)]
    data$Additional_metadata <- hierarchy$Additional_metadata[match(data$Org,
      hierarchy$Org)]

    return(data)
  }
}
