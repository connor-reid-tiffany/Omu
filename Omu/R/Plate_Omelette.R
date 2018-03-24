#'plate_omelette
#Internal method for KEGG_Gather
#'@export



plate_omelette <- function(countDF, Matrix) UseMethod("plate_omelette")

#' @rdname plate_omelette
#' @export
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

#' @rdname plate_omelette
#' @export
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
