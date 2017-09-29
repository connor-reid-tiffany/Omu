#Method for cleaning output of KEGG_Gather using regex and strsplit



Do_Dishes <- function(countDF) UseMethod("Do_Dishes", countDF)

Do_Dishes.rxn <- function(countDF) "rxn"{

#Clean up using regex
countDF$Rxn = gsub("\\c\\("," \\[",countDF$Rxn)
countDF$Rxn = gsub("[[:punct:]]", "", countDF$Rxn)
countDF$Rxn = gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",countDF$Rxn))

#Split multiple values per identifier by a delimiter
lst <- strsplit(as.character(countDF$Rxn), ",")
Split_DF <- transform(countDF[rep(1:nrow(countDF),
lengths(lst)),-5], Rxn= unlist(lst))
class(DF_gathered) <- append(class(DF_gathered), "rxn")
return(Split_DF}
}

Do_Dishes.genes <- function(countDF) "genes"{

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
colnames(DF2)[1] <- "Org"
DF_Complete <- data.frame(do.call('rbind', strsplit(as.character(DF2$Org),':',fixed=TRUE)))
DF_Complete = as.data.frame(cbind(DF, DF_Complete))
colnames(DF_Complete)[17] <- "Org"
DF_Complete$Org <- tolower(DF_Complete$Org)
DF_Complete$Org <- trimws(DF_Complete$Org)
DF_Complete = DF_Complete[ , -which(names(DF_Complete) %in% c("Gene_number"))]
return(DF_Complete)
}
