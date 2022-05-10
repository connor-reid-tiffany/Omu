#' Gather metadata from KEGG for metabolites
#'
#' @description Method for gathering metadata from the KEGG API.
#' @param count_data A metabolmics count dataframe with a KEGG identifier columns
#' @importFrom dplyr if_else
#' @examples
#' count_data <- assign_hierarchy(count_data = c57_nos2KO_mouse_countDF,
#' keep_unknowns = TRUE, identifier = "KEGG")
#'
#' count_data <- subset(count_data, Subclass_2=="Aldoses")
#'
#' count_data <- KEGG_gather(count_data = count_data)
#' @export

KEGG_gather <- function(count_data) UseMethod("KEGG_gather", count_data)

#' @rdname KEGG_gather
#' @export
KEGG_gather.cpd <- function(count_data){

  #stop conditions
  if(any(names(count_data) %in% "KEGG")==FALSE){

    stop("dataframe is missing a KEGG compound number column")

  }
  #Set variables
  first_char <- "C"
  column <- "KEGG"


  #print patience snail into terminal =)
  text_art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
                 "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
                 "  \\ _______________/______\\................please be patient =)"
  )
  cat(text_art, sep = "\n")

  #Send identifier count_data to KEGG API
  df <- make_omelette(count_data = count_data, column = column, first_char = first_char)
  class(df)[2] <- "rxn"
  #call parser function
  df = plate_omelette(df)
  #Append acquired data
  count_data <- merge(df, count_data, "KEGG")
  #Assign rxn class to data.frame
  class(count_data)[2] <- "rxn"
  #We want Orthologies, so need to run new DF through KEGG_Gather again
  count_data <- KEGG_gather(count_data)

  return(count_data)

  }

#' @rdname KEGG_gather
#' @export
KEGG_gather.rxn <- function(count_data){

#Set variables
first_char <- "R"
column <- "Rxn"


#print patience snail into terminal =)
text_art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
               "  \\ _______________/______\\................please be patient =)"
)
cat(text_art, sep = "\n")


#Send indentifier data to KEGG API
df <- make_omelette(count_data = count_data, column = column, first_char = first_char)
class(df)[2] <- "KO"
#Call Plate_Omelette method to clean data up
df <- plate_omelette(df)

count_data <- merge(df, count_data, "Rxn")
#append KO class in case user wishes to KEGG_Gather genes
class(count_data) <- c("data.frame", "KO")

return(count_data)
}


#' @rdname KEGG_gather
#' @export
KEGG_gather.KO <- function(count_data){

#Set variables
column <- "KO"
first_char <- "K"

#print patience snail into terminal =)
text_art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
               "  \\ _______________/______\\................please be patient =)"
)
cat(text_art, sep = "\n")

#Send indentifier to KEGG API
df <- make_omelette(count_data = count_data, column = column, first_char = first_char)
#append class "genes"
class(df) <- append(class(df), "genes")
#Call Plate_Omelette to make it human readable
df <- plate_omelette(df)
#Append acquired data
count_data <- merge(df, count_data, "KO")

class(df) <- append(class(df), "seqs")

return(count_data)
}
