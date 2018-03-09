#Internal method for KEGG_Gather
#'@export



plate_omelette <- function(countDF, Matrix) UseMethod("plate_omelette")

plate_omelette.rxn <- function(countDF){

#Clean up using regex
countDF$Rxn = gsub("\\c\\("," \\[",countDF$Rxn)
countDF$Rxn = gsub("[[:punct:]]", "", countDF$Rxn)
countDF$Rxn = gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",countDF$Rxn))

#Split multiple values per identifier by a delimiter
lst <- strsplit(as.character(countDF$Rxn), ",")
Split_DF <- transform(countDF[rep(1:nrow(countDF),
lengths(lst)),-5], Rxn= unlist(lst))
Split_DF$Rxn = as.vector(Split_DF$Rxn)
Split_DF[Split_DF =="NULL"] <- NA
Split_DF = Split_DF[which(complete.cases(Split_DF$Rxn)),]

class(Split_DF) <- append(class(Split_DF), "rxn")
return(Split_DF)
}

plate_omelette.genes <- function(countDF){

#Clean up using regex
countDF$Genes = gsub("\\c\\("," \\[",countDF$Genes)
countDF$Genes = gsub("\\[|\\]", "", countDF$Genes)

#Split multiple values per identifier by a delimiter
lst <- strsplit(as.character(countDF$Genes), ",")
DF <- transform(countDF[rep(1:nrow(countDF), lengths(lst)),-1], Genes= unlist(lst))
DF <- as.data.frame(sapply(DF, function(x) gsub("\"", "", x)))
DF$Gene_number <- DF$Genes
DF2 <- str_split_fixed(DF$Genes, ":", 2)
DF2 = as.data.frame(DF2)
colnames(DF2)[2] <- "Genes"
DF = DF[,-which(colnames(DF)=="Genes")]

#Create organism identifier
colnames(DF2)[1] <- "Org"
DF_Complete = as.data.frame(cbind(DF2, DF))

DF_Complete$Org <- tolower(DF_Complete$Org)
DF_Complete$Org <- trimws(DF_Complete$Org)
DF_Complete = DF_Complete[ , -which(names(DF_Complete) %in% c("Gene_number"))]
return(DF_Complete)
}

plate_omelette.rxnKO <- function(countDF, Matrix){

DF = as.data.frame(unlist(Matrix[,2], recursive = F))
DF = rownames_to_column(DF, "KO_Number")
DF = with(DF, cbind(DF[,2],
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
countDF= left_join(KO_with_Rxn, countDF, by = "Rxn")

#Drop the element number column as it is no longer needed
countDF = countDF[ , !(names(countDF) %in% "Element_number")]

return(countDF)
}
