#Method for gathering metadata from KEGG API


KEGG_Gather <- function(countDF, sig_threshold) UseMethod("KEGG_Gather", countDF)

KEGG_Gather.cpd <- function(countDF, sig_threshold) "cpd"{

if (missing(sig_threshold)){
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')
  }else {
    data = data[which(data$padj <= sig_threshold ), ]
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')}
  #Set conditions for classes "cpd" and "KO"
  if(class(countDF)[2] = "cpd"){
  Reqs = c('ENTRY', 'REACTION')
  column = "KEGG"
  matchval = "KEGG"
  req_name = "REACTION"
  appended_col = "Rxn"
  } else if(class(countDF)[2] = "KO"{
  Reqs = c("NAME", "ENTRY","DEFINITION", "GENES")
  column = "KO_Number"
  match_val = "KO_Number"
  req_name = "GENES"
  appended_col = "Genes"
  }
  #Create input of nested lists of length 10 and feed to KEGG API         
  Input <- data.frame(i = test_metabo[, column])
  Input = as.vector(Input[ grep("nan", Input[, column], invert = TRUE) , ])
  Input_Split  <- split(Input,  ceiling(seq_along(Input)/10))
  Input_Split = llply(Input_Split, as.list)
  #Send nested list of Cpd numbers to KEGG API
  Output <- llply(Input_Split, function(x)keggGet(x))
  Unlist_Output <- unlist(Output, recursive = F)
  #Pull out Rxn and Cpd numbers
  Extract_reqs <- lapply(Unlist_Output, '[', Reqs)
  #Convert data into list of character vectors that can be added to an R data frame object
  #Note that character vectors 
  n_obs <- sapply(Extract_reqs, length)
  seq_max <- seq_len(max(n_obs))
  Matrix <- t(sapply(Extract_reqs, '[', i = seq_max)
  DF = as.data.frame(Matrix)
  countDF[,appended_col] = DF[,req_name][match(data[match_val], DF[, match_val])]
  
  if(appended_col=="Rxn"){
    class(countDF) <- append(class(countDF), "rxn")
    } else if(appended_col=="genes"){
    class(countDF) <- append(class(countDF), "genes"
       }
  #Call function from method Clean_Up
  countDF = Clean_Up(countDF)
  
  return(countDF)
                             }


KEGG_Gather.rxn <- function(countDF) "rxn"

KEGG_Gather.genes <- function(countDF) "KO"
