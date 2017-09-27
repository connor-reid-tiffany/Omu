#Method for cleaning output of KEGG_Gather using regex and strsplit



Clean_Up <- function(countDF) UseMethod("Clean_Up", countDF)

Clean_Up.rxn <- function(countDF) "rxn"{
  data$Rxn = gsub("\\c\\("," \\[",data$Rxn)
  data$Rxn = gsub("[[:punct:]]", "", data$Rxn)
  data$Rxn = gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",data$Rxn))
  
  lst <- strsplit(as.character(data$Rxn), ",")
  DF_gathered <- transform(data[rep(1:nrow(data), lengths(lst)),-5], Rxn= unlist(lst))
  class(DF_gathered) <- append(class(DF_gathered), "rxn")
  return(DF_gathered)}
  }

Clean_Up.genes <- function(countDF) "genes"{
countDF$Genes = gsub("\\c\\("," \\[",countDF$Genes)
countDF$Genes = gsub("\\[|\\]", "", countDF$Genes)

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
