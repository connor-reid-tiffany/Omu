#' KEGG_gather
#' Method for gathering metadata from KEGG API
#' @param count_data A metabolmics count dataframe with a KEGG identifier columns
#' @importFrom dplyr if_else
#' @export

KEGG_gather <- function(count_data) UseMethod("KEGG_gather", count_data)

#' @rdname KEGG_gather
#' @export
KEGG_gather.cpd <- function(count_data){

  #Set variables
  req <- c('ENTRY', 'REACTION')
  column <- "KEGG"
  req_name <- "REACTION"

  #print patience snail into terminal =)
  text_art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
                 "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
                 "  \\ _______________/______\\................please be patient =)"
  )
  cat(text_art, sep = "\n")

  #Send identifier count_data to KEGG API
  matrix <- make_omelette(count_data = count_data, column = column, req = req)

  #Convert to data.frame and append acquired data
  df = as.data.frame(matrix)
  count_data$Rxn = df[,req_name][match(count_data[, column], df[, 'ENTRY'])]

  #Assign rxn class to data.frame
  class(count_data)[2] <- "rxn"

  #Call function from method Plate_Omelette to make data human readable
  count_data = plate_omelette(count_data)

  #We want Orthologies, so need to run new DF through KEGG_Gather again

  count_data <- KEGG_gather(count_data)

  return(count_data)

  }

#' @rdname KEGG_gather
#' @export
KEGG_gather.rxn <- function(count_data){

#Set variables
req <- c("ENTRY","ORTHOLOGY")
column <- "Rxn"
req_name <- "Rxn"


#print patience snail into terminal =)
text_art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
               "  \\ _______________/______\\................please be patient =)"
)
cat(text_art, sep = "\n")


#Send indentifier data to KEGG API
matrix <- make_omelette(count_data = count_data, column = column, req = req)

#append rxnKO class for calling Plate_Omelette

#Call Plate_Omelette method to clean data up
count_data = plate_omelette_rxnko(count_data, matrix)

#append KO class in case user wishes to KEGG_Gather genes
class(count_data) <- append(class(count_data), "KO")

return(count_data)
}


#' @rdname KEGG_gather
#' @export
KEGG_gather.KO <- function(count_data){

#Set variables
req <- c("NAME", "ENTRY","DEFINITION", "GENES")
column <- "KO_Number"
req_name <- "GENES"

#print patience snail into terminal =)
text_art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
               "  \\ _______________/______\\................please be patient =)"
)
cat(text_art, sep = "\n")

#Send indentifier to KEGG API
matrix <- make_omelette(count_data = count_data, column = column, req = req)
df = as.data.frame(matrix)

#append class "genes"
class(count_data) <- append(class(count_data), "genes")

#append columns
count_data$Genes = df[,req_name][match(count_data[, column], df$ENTRY)]
count_data$GeneOperon = df$NAME[match(count_data[, column], df$ENTRY)]

#Call Plate_Omelette to make it human readable
count_data = plate_omelette(count_data = count_data)

return(count_data)
}
