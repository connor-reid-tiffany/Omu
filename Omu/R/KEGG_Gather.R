#Method for gathering metadata from KEGG API


KEGG_Gather <- function(countDF, sig_threshold) UseMethod("KEGG_Gather", countDF)

KEGG_Gather.cpd <- function(countDF, sig_threshold){

#create value column, subset data based on significance
if (missing(sig_threshold)){
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')
  }else {
    data = data[which(data$padj <= sig_threshold ), ]
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')}

  #Set variables
  Reqs = c('ENTRY', 'REACTION')
  column = "KEGG"
  req_name = "REACTION"
  appended_col = "Rxn"

  #print patience snail into terminal =)
  Text_Art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\",
                 "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____",
                 "  \\ _______________/______\\................please be patient =)"
  )
  cat(Text_Art, sep = "\n")

  #Send identifier data to KEGG API
  Matrix <- Make_Omelette(countDF = countDF, column = column, Reqs = Reqs)

  #Convert to data.frame and append acquired data
  DF = as.data.frame(Matrix)
  countDF[,appended_col] = DF[,req_name][match(countDF[, column], DF[, column])]

  #Assign rxn class to data.frame
  class(countDF) <- append(class(countDF), "rxn")

  #Call function from method Do_Dishes to make data human readable
  countDF = Do_Dishes(countDF)

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

#Extract Element_number and KO_Number from rownames and put into 2 columns
DF = as.data.frame(unlist(Matrix[,2], recursive = F))
DF = rownames_to_column(KO_df, "KO_Number")  DF = with(DF, cbind(DF[,2],
colsplit(DF$KO_Number, pattern = "\\.",
names = c('Element_number', 'KO_Number'))))
colnames(DF)[1] <- "KO"

#Do the above but for Rxn numbers this time
Rxn_df = as.data.frame(unlist(Matrix[,1], recursive = F))
Rxn_df = rownames_to_column(Rxn_df, "Element_number")
Rxn_df = with(Rxn_df, cbind(Rxn_df[,2],
  colsplit(Rxn_df$Element_number, pattern = "\\.",
  names = c('Element_number', 'remove'))))
colnames(Rxn_df)[1] <- "Rxn"
Rxn_df = Rxn_df[,1:2]

#Join origianl dataframe to dataframe with KO by Rxn number
KO_with_Rxn <- left_join(DF, Rxn_df, by = "Element_number")
data = left_join(KO_with_Rxn, data, by = "Rxn")

#Drop the element number column as it is no longer needed
data = data[ , !(names(data) %in% "Element_number")]

return(data)
}



KEGG_Gather.KO <- function(countDF){

#Set variables
Reqs = c("NAME", "ENTRY","DEFINITION", "GENES")
column = "KO_Number"
req_name = "GENES"
appended_col = "Genes"


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
countDF[,appended_col] = DF[,req_name][match(countDF[, column], DF[, column])]
countDF$GeneOperon = DF$NAME[match(DF[, column], DF[, column])]

#Call Do_Dishes to make it human readable
countDF = Do_Dishes(countDF = countDF)

return(countDF)
}
