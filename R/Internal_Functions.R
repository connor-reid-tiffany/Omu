#' Get metadata from KEGG API
#'
#' @description Internal function for KEGG_Gather
#' @param count_data The metabolomics count data
#' @param column The name of the KEGG identifier being sent to the KEGG API
#' @param first_char firct character in number being fed to KEGG database
#' @importFrom httr GET
#' @importFrom httr message_for_status
#' @importFrom httr content
#' @export


make_omelette <- function(count_data, column, first_char){


 strip <- function(str)
  {
    gsub("^\\s+|\\s+$", "", str)
  }

 column <- column
 first_char <- first_char

 input <- data.frame(i = count_data[, column])
 #remove values that do not have a KEGG cpd
 input = data.frame(i = input[!is.na(input$i),])
 #remove values that are not 6 characters long
 input$char <- nchar(as.matrix(input))
 input <- input[input$char==6,]
 #remove values that do not begin with first_char
 input$first_char <- substring(as.matrix(input$i), 1,1)
 input <- input[input$first_char==first_char,]
 #convert to character vector and split into a list of length 10 character vectors to feed to KEGG API
 input <- as.vector(input$i)
 input_split  <- split(input,  ceiling(seq_along(input)/10))

 #function to pull flat text files from KEGG API. will encapsulate errors that are due to
 #bad requests, lost connection, or server issues and return a NULL instead of throwing an error
 kegg_get <- function(x){

  x <- paste(x, collapse="+")
  url <- sprintf("%s/get/%s", getOption("KEGG_REST_URL", "http://rest.kegg.jp"), x)
  response <- GET(url)
  status_report <- tryCatch(message_for_status(GET(url)),
           http_404 = function(c) "That url doesn't exist",
           http_403 = function(c) "Authentication required",
           http_400 = function(c) "Incorrect input from KEGG column",
           http_500 = function(c) "KEGG server is unavailable"
  )
  content <- .strip(content(response, "text"))
  if (nchar(content) == 0)
    return(status_report)
  return(content)
 }

 #iterate over the list of vectors
  output <- lapply(input_split, kegg_get)
  #control flow to determine if server was successfully reached.
  if (any(c("That url doesn't exist","Authentication required","Incorrect input from KEGG column","KEGG server is unavailable") %in% output[[1]])==FALSE){

    output <- output

  } else if (identical(output[[1]], "That url doesn't exist")==TRUE){

    message("Ending process because url does not exist.")

    return(NULL)

  } else if (identical(output[[1]], "Authentication required")==TRUE){

    message("Ending process because Authentication is required.")

    return(NULL)

  } else if (identical(output[[1]], "Incorrect input from KEGG column")==TRUE){

    message("Ending process because KEGG column in count_data did not contain valid KEGG cpd numbers.")

    return(NULL)

  } else if (identical(output[[1]], "KEGG server is unavailable")==TRUE){

    message("Ending process because KEGG server is unavailable.")

    return(NULL)

  }

  return(output)

}

#' Clean up orthology metadata
#'
#' @description Internal function for KEGG_Gather.rxn method
#' KEGG_Gather.rxn requires dispatch on multiple elements, so
#' There was no way to incorporate as a method
#' @param output output from plate_omelette
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @export

plate_omelette_rxnko<- function(output){

  content <- lapply(output, function(x) strsplit(strip(x), "\n", fixed=TRUE)[[1]])
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
    content <- lapply(content, function(x) sapply(x, function(x) x[x!=""]))
    content <- lapply(content, as.list)

    content_rxn <- lapply(content, names)
    content_rxn <- lapply(content_rxn, as.list)

    content<- lapply(rapply(content, enquote, how="unlist"), eval)
    content_rxn<- lapply(rapply(content_rxn, enquote, how="unlist"), eval)

    content <- lapply(content, as.data.frame)
    content <- lapply(content, function(x){ x$first_char <- substring(x[,1], 1,1); return(x)})
    content <- lapply(content, function(x) {x <- x[x$first_char=="K",]; return(x)})
    content <- lapply(content, function(x) {x$char <- nchar(x[,1]); x[x$char==6,]; x <- x[,!names(x) %in% c("first_char", "char")]; return(x)})

    content <- Map("rbind", content, content_rxn)

    content <- lapply(content, as.data.frame(t))

    content <- lapply(content, function(x){ colnames(x) <- c("KO", "Rxn"); return(x)})

    content <- do.call("rbind", content)

  return(content)
}
