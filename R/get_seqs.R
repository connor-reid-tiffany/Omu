#' Get nucleotide and amino acid sequences for genes
#'
#' @description Function that gets nt and aa seqs for gene data from KEGG_gather
#' @param gene_data A dataframe with genes from KEGG_gather, with class seqs
#' @importFrom httr GET
#' @importFrom stringr str_split
#' @importFrom stringr str_detect
#' @examples
#'
#' gene_data <- c57_nos2KO_mouse_countDF[(1:2),]
#'
#' gene_data <- KEGG_gather(gene_data)
#'
#' gene_data <- KEGG_gather(gene_data)
#' gene_data <- gene_data[1:2,]
#'
#' gene_data <- get_seqs(gene_data)
#' @export

get_seqs <- function(gene_data){


  .strip <- function(str)
  {
    gsub("^\\s+|\\s+$", "", str)
  }

  gene_data$Gene_code <- paste0(gene_data$Org, ":", gene_data$Genes)

  #remove any gene name by targetting parenthesis and characters after with regex
  gene_data$Gene_code <- gsub(pattern = "\\(.*", replacement = "", x = gene_data$Gene_code)

  input <- data.frame(i = gene_data$Gene_code)
  #remove values that do not have a KEGG cpd
  input = data.frame(i = input[!is.na(input$i),])
  #convert to character vector and split into a list of length 10 character vectors to feed to KEGG API
  input <- as.vector(input$i)
  input_split  <- split(input,  ceiling(seq_along(input)/10))

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

    message("Ending process because KEGG column in gene_data did not contain valid KEGG cpd numbers.")

    return(NULL)

  } else if (identical(output[[1]], "KEGG server is unavailable")==TRUE){

    message("Ending process because KEGG server is unavailable.")

    return(NULL)

  }





  content <- lapply(output, function(x) strsplit(.strip(x), "\n", fixed=TRUE)[[1]])
  #replace delimeter elements with END_OF_ENTRY to separate entries
  content <- lapply(content, function(x) gsub(x, pattern = "///", replacement = "END_OF_ENTRY"))
  #need to keep ENTRY, NAME, AASEQ, and NTSEQ fields
  #content <- paste(content, sep = "", collapse = "")
  content <- lapply(content, function(x) paste(x, sep = "", collapse = ""))
  #split into character matrix by End of Entry
  #content <- str_split(content, "END_OF_ENTRY", simplify = TRUE)
  content <- lapply(content, function(x) str_split(x, "END_OF_ENTRY", simplify = TRUE))
  #remove elements that don't contain NTSEQ (broken record but needs control flow for each class)
  #content <-t(content[,str_detect(content, pattern = "REACTION")==TRUE])
  content <- lapply(content, function(x) t(x[,str_detect(x, pattern = "NTSEQ")==TRUE]))
  #convert each column into a vector within a list
  #content <- as.list(as.data.frame(content))
  #change element names to compound, this will need control flow in the future for each class of cpd, rxn, KO because the word after ENTRY will
  #be different!
  change_names <- function(x){

    colnames(x) <- gsub('^.*ENTRY\\s*|\\s*ENTRY.*$', '', x)

    return(x)

  }
  content <- lapply(content, change_names)
  content <- lapply(content, function(x) as.list(as.data.frame(x)))

  make_seq_dfs <- function(x){

    df <- data.frame(Gene_code = regmatches(x, regexpr(pattern = "(ENTRY).+?(?=CDS)", text = x, perl = TRUE)),
                     NAME = regmatches(x, regexpr(pattern = "(NAME).+?(?=ORTHOLOGY|ORGANISM)", text = x, perl = TRUE)),
                     NTSEQ = regmatches(x, regexpr(pattern = "NTSEQ.*", text = x, perl = TRUE)),
                     AASEQ = regmatches(x, regexpr(pattern = "(AASEQ).+?(?=NTSEQ)", text = x, perl = TRUE)))

  }



  content <- lapply(content, make_seq_dfs)

  content <- do.call("rbind", content)

  content$Gene_code <- gsub(pattern = "ENTRY", replacement = "", x = content$Gene_code)

  content$Gene_code <- trimws(content$Gene_code, "both")

  content$NAME <- gsub(pattern = "NAME", replacement = "", x = content$NAME)

  content$NAME <- trimws(content$NAME, "both")

  gene_data$Gene_code <- gsub(pattern = ".*:", replacement = "", x = gene_data$Gene_code)

  gene_data_m <- merge(gene_data, content, by = "Gene_code")

  return(gene_data_m)

}
