#Method for gathering metadata from KEGG API
#'@export

KEGG_gather <- function(countDF, sig_threshold) UseMethod("KEGG_gather", countDF)

KEGG_gather.cpd <- function(countDF, sig_threshold){

#create value column, subset data based on significance
if (missing(sig_threshold)){
    countDF$Val <- if_else(countDF$log2FoldChange > 0, 'Increase', 'Decrease')
  }else {
    countDF = countDF[which(countDF['padj'] <= sig_threshold), ]
    countDF$Val <- if_else(countDF["log2FoldChange"] > 0, 'Increase', 'Decrease')}

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

  #Send identifier countDF to KEGG API
  matrix <- Make_Omelette(countDF = countDF, column = column, req = req)

  #Convert to data.frame and append acquired data
  df = as.data.frame(matrix)
  countDF$Rxn = df[,req_name][match(countDF[, column], df[, 'ENTRY'])]

  #Assign rxn class to data.frame
  class(countDF)[2] <- "rxn"

  #Call function from method Plate_Omelette to make data human readable
  countDF = plate_omelette(countDF)

  #We want Orthologies, so need to run new DF through KEGG_Gather again

  countDF <- KEGG_gather(countDF)

  return(countDF)

  }


KEGG_gather.rxn <- function(countDF){

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
matrix <- make_omelette(countDF = countDF, column = column, req = req)

#append rxnKO class for calling Plate_Omelette

#Call Plate_Omelette method to clean data up
countDF = plate_omelette_rxnko(countDF, matrix)

#append KO class in case user wishes to KEGG_Gather genes
class(countDF) <- append(class(countDF), "KO")

return(countDF)
}



KEGG_gather.KO <- function(countDF){

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
matrix <- make_omelette(countDF = countDF, column = column, req = req)
df = as.data.frame(matrix)

#append class "genes"
class(countDF) <- append(class(countDF), "genes")

#append columns
countDF$Genes = df[,req_name][match(countDF[, column], df$ENTRY)]
countDF$GeneOperon = df$NAME[match(countDF[, column], df$ENTRY)]

#Call Plate_Omelette to make it human readable
countDF = plate_omelette(countDF = countDF)

return(countDF)
}
