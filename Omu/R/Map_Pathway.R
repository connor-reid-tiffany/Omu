#'Map_Orthology
#'Map functional orthologies, i.e. enzymes, based on compounds and their rxns in KEGG database
#'@param data is a metabolomics count data frame with a column of KEGG compound numbers
#'@param sig_threshold is an optional significance threshold to be used as a validation step
#'@param ortho_hierarchy is a path to the csv file of the KEGG brite hierarchy of orthologies
#'@export
#'@example Map_Orthology(data = yourdataframe, sig_threshold = 0.05, Orthology_Hierarchy = Orthology_Hierarchy)



Map_Orthology <- function(data, sig_threshold, Orthology_Hierarchy){
OG_data = data
if (missing(sig_threshold)){
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')
  }else {
    data = data[which(data$padj <= sig_threshold ), ]
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')}

#First need to get reaction numbers from KEGG API
Map_Rxns <- function(data){
  #create a nested list of KEGG numbers, ten elements each due to server side limitations of API
  Vector <- as.data.frame(data$KEGG)
  colnames(Vector)[1] <- "KEGG"
  Vector2 = as.vector(Vector[ grep("nan", Vector$KEGG, invert = TRUE) , ])
  Vector_Split  <- split(Vector2,  ceiling(seq_along(Vector2)/10))
  Vector_Split2 = llply(Vector_Split, as.list)
  #Send nested list of Cpd numbers to KEGG API
  Cpd_to_Rxn <- llply(Vector_Split2, function(x)keggGet(x))
  Unlist_Rxn <- unlist(Cpd_to_Rxn, recursive = F)
  #Pull out Rxn and Cpd numbers
  Rxn_DF <- lapply(Unlist_Rxn, '[', c("ENTRY","REACTION"))
  #Convert data into list of character vectors that can be added to an R data frame object
  #Note that character vectors 
  n_obs_Rxn <- sapply(Rxn_DF, length)
  seq_max_Rxn <- seq_len(max(n_obs_Rxn))
  Rxn_Matrix <- t(sapply(Rxn_DF, '[', i = seq_max_Rxn))
  colnames(Rxn_Matrix)[1] <- "KEGG"
  Rxn_Matrix = as.data.frame(Rxn_Matrix)
  data$Rxn = Rxn_Matrix$REACTION[match(data$KEGG, Rxn_Matrix$KEGG)]
  #For unidentified compounds, Null values need to be accounted for and changed to NA
  data <- data[complete.cases(data$KEGG),]
  data$Rxn = as.character(data$Rxn)
  data <- data[complete.cases(data$Rxn),]
  
  data[data == "NULL"] <- NA
  data[data == "NA"] <- NA
  data <- data[complete.cases(data$Rxn),]
  #Use REGEX to remove special characters except for commas, and using commas create individual rows for each Rxn
  #This will duplicate the data in other rows that was associated with those character vectors for each Rxn number
  data$Rxn = gsub("\\c\\("," \\[",data$Rxn)
  data$Rxn = gsub("[[:punct:]]", "", data$Rxn)
  data$Rxn = gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",data$Rxn))
  
  lst <- strsplit(as.character(data$Rxn), ",")
  DF <- transform(data[rep(1:nrow(data), lengths(lst)),-5], Rxn= unlist(lst))
  
  return(DF)
}

data_rxn <- Map_Rxns(data)


Rxn_to_Orthology <- function(data){
  
  #Convert column of Rxn numbers into a nested list. Each list needs 10 entries
  DRxn = as.list(data$Rxn)
  DRxn <- split(DRxn,  ceiling(seq_along(DRxn)/10))
  #Send list to API, pull out required data
  KO <- llply(DRxn, function(x)keggGet(x))
  RX_KO <- as.data.frame(cbind(DRxn, KO))
  KO = unlist(KO, recursive = F)
  KO = lapply(KO, '[', c("ENTRY", "ORTHOLOGY"))
   #Unlist and convert to matrix
  n_obs_KO <- sapply(KO, length)
  seq_max_KO <- seq_len(max(n_obs_KO))
  KO_Matrix <- t(sapply(KO, '[', i = seq_max_KO))
  #Convert to data frame, unlist and obtain names attribute from KO values
  #These are the KO numbers needed to assign a hierarchy. Then Create
  #a second data frame of Rxn numbers
  #finally join the two tables by the shared element numbers
  #Joining the Rxn to KO is necessary to join this dataframe with the original data frame
  KO_df = as.data.frame(unlist(KO_Matrix[,2], recursive = F))
  KO_df = rownames_to_column(KO_df, "KO_Number")
  KO_df = with(KO_df, cbind(KO_df[,2], colsplit(KO_df$KO_Number, pattern = "\\.", names = c('Element_number', 'KO_Number'))))
  colnames(KO_df)[1] <- "KO"
  
  Rxn_df = as.data.frame(unlist(KO_Matrix[,1], recursive = F))
  Rxn_df = rownames_to_column(Rxn_df, "Element_number")
  Rxn_df = with(Rxn_df, cbind(Rxn_df[,2], colsplit(Rxn_df$Element_number, pattern = "\\.", names = c('Element_number', 'remove'))))
  colnames(Rxn_df)[1] <- "Rxn"
  Rxn_df = Rxn_df[,1:2]
  #Join origianl dataframe to dataframe with KO by Rxn number
  KO_with_Rxn <- left_join(KO_df, Rxn_df, by = "Element_number")
  data = left_join(KO_with_Rxn, data, by = "Rxn")
  #Drop the element number column as it is no longer needed
  data = data[ , !(names(data) %in% "Element_number")]
  return(data)
}
data_KO <- Rxn_to_Orthology(data_rxn)

data_ortho = Orthology_hierarchy(data = data_ortho, ortho_map = ortho_map)
Order = c("Metabolite","Val","log2FoldChange", "Class", "Subclass_1", "Subclass_2", "Subclass_3", 
          "Subclass_4", "Rxn", "KEGG", "Orthology_number", "KO", "KO_Class", "KO_Sub_class1", "KO_Sub_class2")
data_ortho = data_ortho[, Order]

return(data_ortho)
}

}
#If a Count table of significant changes by Orthology Hierarchy meta is required
Count_Changes <- function(data, ..., column, alpha){
  
  data = data %>% group_by_(...) %>%
    mutate(Significant_Changes = sum(log2FoldChange>0),
           neg = sum(log2FoldChange<0))
  data = data[,c(column, "Significant_Changes", "neg")]
  data$neg <- data$neg * -1
  output <- data[,c(1,3)]
  colnames(output)[2] <- "Significant_Changes"
  data <- rbind(data, output)
  data = subset(data, select = -neg)
  
  data$colour = ifelse(data$Significant_Changes < 0, "Decrease","Increase")
  unique(data[])
  data = data[apply(data[2],1,function(z) !any(z==0)),]
}





