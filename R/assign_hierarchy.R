#' Assign hierarchy metadata
#'
#' @description Assigns hierarchy metadata to a metabolomics count matrix using identifier values.
#' It can assign KEGG compound hierarchy, orthology hierarchy, or organism hierarchy data.
#' @param count_data a metabolomics count data frame with either a KEGG compound, orthology,
#' or a gene identifier column
#' @param keep_unknowns a boolean of either TRUE or FALSE. TRUE keeps unannotated compounds,
#' FALSE removes them
#' @param identifier a string that is either "KEGG" for metabolite, "KO" for orthology,
#' "Prokaryote" for organism, or "Eukaryote" for organism
#' @importFrom stats complete.cases
#' @examples
#' assign_hierarchy(count_data = c57_nos2KO_mouse_countDF, keep_unknowns = TRUE, identifier = "KEGG")
#' @export

assign_hierarchy <- function(count_data, keep_unknowns, identifier){


  Metabolite <- NULL
  identifier = match.arg(arg = identifier, choices = c("KEGG", "KO", "Prokaryote", "Eukaryote"))

  if (identifier == "KEGG"){

    if (any(names(test_df) %in% "KEGG"!=FALSE)){

      stop("dataframe is missing KEGG compound number column")

    }
    if (keep_unknowns ==FALSE){

      count_data <- count_data[complete.cases(count_data$KEGG),]
      count_data$Class <- Metabolite_Hierarchy_Table$Class[match(count_data$KEGG,
                                                                   Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_1 <- Metabolite_Hierarchy_Table$Subclass_1[match(count_data$KEGG,
                                                                 Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_2 <- Metabolite_Hierarchy_Table$Subclass_2[match(count_data$KEGG,
                                                                 Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_3 <- Metabolite_Hierarchy_Table$Subclass_3[match(count_data$KEGG,
                                                                 Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_4 <- Metabolite_Hierarchy_Table$Subclass_4[match(count_data$KEGG,
                                                                 Metabolite_Hierarchy_Table$KEGG)]


      class(count_data) = append(class(count_data), "cpd")
      return(count_data)
    }
    else if (keep_unknowns == TRUE) {
      count_data$Class <- Metabolite_Hierarchy_Table$Class[match(count_data$KEGG,
                                                                 Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_1 <- Metabolite_Hierarchy_Table$Subclass_1[match(count_data$KEGG,
                                                                           Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_2 <- Metabolite_Hierarchy_Table$Subclass_2[match(count_data$KEGG,
                                                                           Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_3 <- Metabolite_Hierarchy_Table$Subclass_3[match(count_data$KEGG,
                                                                           Metabolite_Hierarchy_Table$KEGG)]
      count_data$Subclass_4 <- Metabolite_Hierarchy_Table$Subclass_4[match(count_data$KEGG,
                                                                           Metabolite_Hierarchy_Table$KEGG)]

      class(count_data) = append(class(count_data), "cpd")
      return(count_data)
    }
  } else if (identifier == "KO"){

    if (any(names(test_df) %in% "KO"!=FALSE)){

      stop("dataframe is missing KO number column")

    }

    count_data$KO_Class <- Orthology_Hierarchy_Table$KO_Class[match(count_data$KO,
                                                                    Orthology_Hierarchy_Table$KO_Number)]
    count_data$KO_Subclass_1 <- Orthology_Hierarchy_Table$KO_Subclass_1[match(count_data$KO,
                                                                              Orthology_Hierarchy_Table$KO_Number)]
    count_data$KO_Subclass_2 <- Orthology_Hierarchy_Table$KO_Subclass_2[match(count_data$KO,
                                                                              Orthology_Hierarchy_Table$KO_Number)]
    count_data$KO_Subclass_3_Enzyme <- Orthology_Hierarchy_Table$KO_Subclass_3_Enzyme[match(count_data$KO,
                                                                                            Orthology_Hierarchy_Table$KO_Number)]

    return(count_data)
  } else if (identifier == "Prokaryote"){

    if (any(names(test_df) %in% "Org"!=FALSE)){

      stop("dataframe is missing Org column")

    }

    count_data$Kingdom <- Prokaryote_Hierarchy_Table$Kingdom[match(count_data$Org,
                                                               Prokaryote_Hierarchy_Table$Org)]
    count_data$Phylum.Class.Family <- Prokaryote_Hierarchy_Table$Phylum.Class.Family[match(count_data$Org,
                                                                   Prokaryote_Hierarchy_Table$Org)]
    count_data$Genus <- Prokaryote_Hierarchy_Table$Genus[match(count_data$Org,
                                                                   Prokaryote_Hierarchy_Table$Org)]
    count_data$Species.Strain.Serotype <- Prokaryote_Hierarchy_Table$Species.Strain.Serotype[match(count_data$Org,
                                                                   Prokaryote_Hierarchy_Table$Org)]


    return(count_data)
  } else if (identifier == "Eukaryote"){

    if (any(names(test_df) %in% "Org"!=FALSE)){

      stop("dataframe is missing Org column")

    }

    count_data$Kingdom <- Eukaryote_Hierarchy_Table$Kingdom[match(count_data$Org,
                                                                    Eukaryote_Hierarchy_Table$Org)]
    count_data$Phylum <- Eukaryote_Hierarchy_Table$Phylum[match(count_data$Org,
                                                                  Eukaryote_Hierarchy_Table$Org)]
    count_data$Class <- Eukaryote_Hierarchy_Table$Class[match(count_data$Org,
                                                                  Eukaryote_Hierarchy_Table$Org)]
    count_data$Genus.Species <- Eukaryote_Hierarchy_Table$Genus.Species[match(count_data$Org,
                                                                  Eukaryote_Hierarchy_Table$Org)]
    count_data$Common.Name <- Eukaryote_Hierarchy_Table$Common.Name[match(count_data$Org,
                                                                  Eukaryote_Hierarchy_Table$Org)]
    count_data$X <- Eukaryote_Hierarchy_Table$X[match(count_data$Org,
                                                                  Eukaryote_Hierarchy_Table$Org)]
    return(count_data)
  }

}
