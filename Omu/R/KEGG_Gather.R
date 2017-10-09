#Method for gathering metadata from KEGG API
#'@export

KEGG_Gather <- function(countDF, sig_threshold) UseMethod("KEGG_Gather", countDF)

KEGG_Gather.cpd <- function(countDF, sig_threshold){

#create value column, subset data based on significance
if (missing(sig_threshold)){
    countDF$Val <- if_else(countDF$log2FoldChange > 0, 'Increase', 'Decrease')
  }else {
    countDF = countDF[which(countDF['padj'] <= sig_threshold), ]
    countDF$Val <- if_else(countDF["log2FoldChange"] > 0, 'Increase', 'Decrease')}

  #Set variables
  Reqs = c('ENTRY', 'REACTION')
  column = "KEGG"
  req_name = "REACTION"

  #print patience snail into terminal =)
  Text_Art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
                 "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
                 "  \\ _______________/______\\................please be patient =)"
  )
  cat(Text_Art, sep = "\n")

  #Send identifier countDF to KEGG API
  Matrix <- Make_Omelette(countDF = countDF, column = column, Reqs = Reqs)

  #Convert to data.frame and append acquired data
  DF = as.data.frame(Matrix)
  countDF$Rxn = DF[,req_name][match(countDF[, column], DF[, 'ENTRY'])]

  #Assign rxn class to data.frame
  class(countDF)[2] <- "rxn"

  #Call function from method Plate_Omelette to make data human readable
  countDF = Plate_Omelette(countDF)

  #We want Orthologies, so need to run new DF through KEGG_Gather again

  New_DF <- KEGG_Gather(countDF)

  return(New_DF)

  }


KEGG_Gather.rxn <- function(countDF){

#Set variables
Reqs = c("ENTRY","ORTHOLOGY")
column = "Rxn"
req_name = "Rxn"


#print patience snail into terminal =)
Text_Art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
               "  \\ _______________/______\\................please be patient =)"
)
cat(Text_Art, sep = "\n")


#Send indentifier data to KEGG API
Matrix <- Make_Omelette(countDF = countDF, column = column, Reqs = Reqs)

#append rxnKO class for calling Plate_Omelette
class(Matrix) <- append(class(Matrix), "rxnKO")
class(countDF)[2] <- "rxnKO"
#Call Plate_Omelette method to clean data up
countDF = Plate_Omelette(countDF,Matrix)

#append KO class in case user wishes to KEGG_Gather genes
class(countDF) <- append(class(countDF), "KO")

return(countDF)
}



KEGG_Gather.KO <- function(countDF){

#Set variables
Reqs = c("NAME", "ENTRY","DEFINITION", "GENES")
column = "KO_Number"
req_name = "GENES"

#print patience snail into terminal =)
Text_Art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
               "  \\ _______________/______\\................please be patient =)"
)
cat(Text_Art, sep = "\n")

#Send indentifier to KEGG API
Matrix <- Make_Omelette(countDF = countDF, column = column, Reqs = Reqs)
DF = as.data.frame(Matrix)

#append class "genes"
class(countDF) <- append(class(countDF), "genes")

#append columns
countDF$Genes = DF[,req_name][match(countDF[, column], DF$ENTRY)]
countDF$GeneOperon = DF$NAME[match(countDF[, column], DF$ENTRY)]

#Call Plate_Omelette to make it human readable
countDF = Plate_Omelette(countDF = countDF)

return(countDF)
}
