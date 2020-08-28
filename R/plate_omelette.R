#' plate_omelette
#' Internal method for KEGG_Gather which parses flat text files
#' @param output The metabolomics count dataframe
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @importFrom stringr str_trim
#' @export



plate_omelette <- function(output) UseMethod("plate_omelette")

#' @rdname plate_omelette
#' @export
plate_omelette.rxn <- function(output){

  .strip <- function(str)
   {
     gsub("^\\s+|\\s+$", "", str)
   }
  #remove newline text delimeters
    content <- lapply(output, function(x) strsplit(.strip(x), "\n", fixed=TRUE)[[1]])
    #replace delimeter elements with END_OF_ENTRY to separate entries
    content <- lapply(content, function(x) gsub(x, pattern = "///", replacement = "END_OF_ENTRY"))
    #convert to a string
    #content <- paste(content, sep = "", collapse = "")
    content <- lapply(content, function(x) paste(x, sep = "", collapse = ""))
    #split into character matrix by End of Entry
    #content <- str_split(content, "END_OF_ENTRY", simplify = TRUE)
    content <- lapply(content, function(x) str_split(x, "END_OF_ENTRY", simplify = TRUE))
    #remove elements that don't contain REACTION (broken record but needs control flow for each class)
    #content <-t(content[,str_detect(content, pattern = "REACTION")==TRUE])
    content <- lapply(content, function(x) t(x[,str_detect(x, pattern = "REACTION")==TRUE]))
    #convert each column into a vector within a list
    #content <- as.list(as.data.frame(content))
    #change element names to compound, this will need control flow in the future for each class of cpd, rxn, KO because the word after ENTRY will
    #be different!
    change_names <- function(x){

      colnames(x) <- gsub('^.*ENTRY\\s*|\\s*Compound.*$', '', x)

      return(x)

    }
    content <- lapply(content, change_names)
    content <- lapply(content, function(x) as.list(as.data.frame(x)))

    #remove everything but REACTION identifiers (again this will need control flow for each class)
    content <- lapply(content, function(x) lapply(x, function(x) gsub('^.*REACTION\\s*|\\s*PATHWAY.*$|MODULE.*$|ENZYME.*$|BRITE.*$|DBLINKS.*$|ATOM.*$|BOND.*$', '', x)))
    #split the strings into vectors of length n again
    content <- lapply(content, function(x) sapply(x, function(x) str_split(x, " ")))
    content <- lapply(content, function(x) sapply(x, function(x) x[x!=""]))
    content <- lapply(content, as.list)

    content_cpd <- lapply(content, names)
    content_cpd <- lapply(content_cpd, as.list)

    content<- lapply(rapply(content, enquote, how="unlist"), eval)
    content_cpd<- lapply(rapply(content_cpd, enquote, how="unlist"), eval)

    content <- Map("rbind", content, content_cpd)

    content <- lapply(content, as.data.frame(t))

    content <- lapply(content, function(x){ colnames(x) <- c("Rxn", "KEGG"); return(x)})

    content <- do.call("rbind", content)

    return(content)
}

#' @rdname plate_omelette
#' @export
plate_omelette.genes <- function(output){

#Clean up using regex
content <- lapply(output, function(x) gsub(x, pattern = "///", replacement = "END_OF_ENTRY"))
  #convert to a string
  #content <- paste(content, sep = "", collapse = "")
  content <- lapply(content, function(x) paste(x, sep = "", collapse = ""))
  #split into character matrix by End of Entry
  #content <- str_split(content, "END_OF_ENTRY", simplify = TRUE)
  content <- lapply(content, function(x) str_split(x, "END_OF_ENTRY", simplify = TRUE))
  #remove elements that don't contain REACTION (broken record but needs control flow for each class)
  #content <-t(content[,str_detect(content, pattern = "REACTION")==TRUE])
  content <- lapply(content, function(x) t(x[,str_detect(x, pattern = "GENES")==TRUE]))
  #convert each column into a vector within a list
  #content <- as.list(as.data.frame(content))
  #change element names to compound, this will need control flow in the future for each class of cpd, rxn, KO because the word after ENTRY will
  #be different!
  change_names <- function(x){

    colnames(x) <- gsub('^.*ENTRY\\s*|\\s*KO.*$', '', x)

    return(x)

  }
  content <- lapply(content, change_names)
  content <- lapply(content, function(x) as.list(as.data.frame(x)))

  #remove everything but REACTION identifiers (again this will need control flow for each class)
  content <- lapply(content, function(x) lapply(x, function(x) gsub('^.*GENES\\s*|\\s*JOURNAL.*$|DOI.*$|SEQUENCE.*$|REFERENCE.*$|AUTHORS.*$|TITLE.*$|JOURNAL.*$', '', x)))
  #split the strings into vectors of length n again
  content <- lapply(content, function(x) sapply(x, function(x) gsub(pattern = ":. *", replacement = ":", x = x)))
  content <- lapply(content, function(x) sapply(x, function(x) str_split(x, "\n")))
  #content <- lapply(content, function(x) sapply(x, function(x) x[x!=""]))
  content <- lapply(content, as.list)

  content_KO <- lapply(content, names)
  content_KO <- lapply(content_KO, as.list)

  content<- lapply(rapply(content, enquote, how="unlist"), eval)
  content_KO<- lapply(rapply(content_KO, enquote, how="unlist"), eval)

  #content <- lapply(content, as.data.frame)
  #content <- lapply(content, function(x){ x$first_char <- substring(x[,1], 1,1); return(x)})
  #content <- lapply(content, function(x) {x <- x[x$first_char=="K",]; return(x)})
  #content <- lapply(content, function(x) {x$char <- nchar(x[,1]); x[x$char==6,]; x <- x[,!names(x) %in% c("first_char", "char")]; return(x)})

  content <- Map("rbind", content, content_KO)

  content <- lapply(content, as.data.frame(t))

  content <- lapply(content, function(x){ colnames(x) <- c("Org", "KO"); return(x)})

  content <- do.call("rbind", content)

  content$Genes <- gsub(".*:","",content$Org)

  content$Genes <- gsub(" ",",",content$Genes)

  content$Org <- gsub(":.*\\s*", "", content$Org)

  content$Org <- str_trim(content$Org)

  content$Org <- tolower(content$Org)

  return(content)
}

#' @rdname plate_omelette
#' @export
plate_omelette.KO <- function(output){

  .strip <- function(str)
   {
     gsub("^\\s+|\\s+$", "", str)
   }
  content <- lapply(output, function(x) strsplit(.strip(x), "\n", fixed=TRUE)[[1]])
    #replace delimeter elements with END_OF_ENTRY to separate entries
    content <- lapply(content, function(x) gsub(x, pattern = "///", replacement = "END_OF_ENTRY"))
    #convert to a string
    #content <- paste(content, sep = "", collapse = "")
    content <- lapply(content, function(x) paste(x, sep = "", collapse = ""))
    #split into character matrix by End of Entry
    #content <- str_split(content, "END_OF_ENTRY", simplify = TRUE)
    content <- lapply(content, function(x) str_split(x, "END_OF_ENTRY", simplify = TRUE))
    #remove elements that don't contain REACTION (broken record but needs control flow for each class)
    #content <-t(content[,str_detect(content, pattern = "REACTION")==TRUE])
    content <- lapply(content, function(x) t(x[,str_detect(x, pattern = "ORTHOLOGY")==TRUE]))
    #convert each column into a vector within a list
    #content <- as.list(as.data.frame(content))
    #change element names to compound, this will need control flow in the future for each class of cpd, rxn, KO because the word after ENTRY will
    #be different!
    change_names <- function(x){

      colnames(x) <- gsub('^.*ENTRY\\s*|\\s*Reaction.*$', '', x)

      return(x)

    }
    content <- lapply(content, change_names)
    content <- lapply(content, function(x) as.list(as.data.frame(x)))

    #remove everything but REACTION identifiers (again this will need control flow for each class)
    content <- lapply(content, function(x) lapply(x, function(x) gsub('^.*ORTHOLOGY\\s*|\\s*DBLINKS.*$|RHEA.*$', '', x)))
    #split the strings into vectors of length n again
    content <- lapply(content, function(x) sapply(x, function(x) str_split(x, " ")))
    #content <- lapply(content, function(x) sapply(x, function(x) x[x!=""]))
    content <- lapply(content, as.list)

    content_rxn <- lapply(content, names)
    content_rxn <- lapply(content_rxn, as.list)

    content<- lapply(rapply(content, enquote, how="unlist"), eval)
    content_rxn<- lapply(rapply(content_rxn, enquote, how="unlist"), eval)

    content <- lapply(content, function(x) as.data.frame(x, stringsAsFactors = FALSE))
    content <- lapply(content, function(x){ x$first_char <- substring(x[,1], 1,1); return(x)})
    content <- lapply(content, function(x) {x <- x[x$first_char=="K",]; return(x)})
    content <- lapply(content, function(x) {x$char <- nchar(x[,1]); x[x$char==6,]; x <- x[,!names(x) %in% c("first_char", "char")]; return(x)})

    content <- Map("rbind", content, content_rxn)

    content <- lapply(content, as.data.frame(t))

    content <- lapply(content, function(x){ colnames(x) <- c("KO", "Rxn"); return(x)})

    content <- do.call("rbind", content)

  return(content)
}
