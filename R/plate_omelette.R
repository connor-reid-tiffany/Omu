#' plate_omelette
#' Internal method for KEGG_Gather
#' @param count_data The metabolomics count dataframe
#' @importFrom stringr str_split_fixed
#' @importFrom stats complete.cases
#' @export



plate_omelette <- function(count_data) UseMethod("plate_omelette")

#' @rdname plate_omelette
#' @export
plate_omelette.rxn <- function(count_data){

#Clean up using regex
count_data$Rxn = gsub("\\c\\("," \\[",count_data$Rxn)
count_data$Rxn = gsub("[[:punct:]]", "", count_data$Rxn)
count_data$Rxn = gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",count_data$Rxn))

#Split multiple values per identifier by a delimiter
lst <- strsplit(as.character(count_data$Rxn), ",")
Split_DF <- transform(count_data[rep(1:nrow(count_data),
lengths(lst)),-5], Rxn= unlist(lst))
Split_DF$Rxn = as.vector(Split_DF$Rxn)
Split_DF[Split_DF =="NULL"] <- NA
Split_DF = Split_DF[which(complete.cases(Split_DF$Rxn)),]

class(Split_DF) <- append(class(Split_DF), "rxn")
return(Split_DF)
}

#' @rdname plate_omelette
#' @export
plate_omelette.genes <- function(count_data){

#Clean up using regex
count_data$Genes = gsub("\\c\\("," \\[",count_data$Genes)
count_data$Genes = gsub("\\[|\\]", "", count_data$Genes)

#Split multiple values per identifier by a delimiter
lst <- strsplit(as.character(count_data$Genes), ",")
DF <- transform(count_data[rep(1:nrow(count_data), lengths(lst)),-1], Genes= unlist(lst))
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
