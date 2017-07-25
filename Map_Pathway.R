Orthology_Mapping <- function(data, sig_threshold, KEGG_key, ortho_map){
OG_data = data
if (missing(sig_threshold)){
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')
  }else {
    data = data[which(data$padj <= sig_threshold ), ]
    data$Val <- if_else(data$log2FoldChange > 0, 'Increase', 'Decrease')}


Rxn_Path <- function(data){


  Vector <- as.data.frame(data$KEGG)
  colnames(Vector)[1] <- "KEGG"
  Vector2 = as.vector(Vector[ grep("nan", Vector$KEGG, invert = TRUE) , ])
  Vector_Split  <- split(Vector2,  ceiling(seq_along(Vector2)/10))
  Vector_Split2 = llply(Vector_Split, as.list)
  
  
  
  KEGG_to_Path <- llply(Vector_Split2, function(x)keggGet(x))
  Unlist_Pathway_Map <- unlist(KEGG_to_Path, recursive = F)
  
  Pathway_DF_Rxn <- lapply(Unlist_Pathway_Map, '[', c("ENTRY","REACTION"))
  
  n_obs_Rxn <- sapply(Pathway_DF_Rxn, length)
  seq_max_Rxn <- seq_len(max(n_obs_Rxn))
  Pathway_Matrix_Rxn <- t(sapply(Pathway_DF_Rxn, '[', i = seq_max_Rxn))
  colnames(Pathway_Matrix_Rxn)[1] <- "KEGG"
  Pathway_Matrix_Rxn = as.data.frame(Pathway_Matrix_Rxn)
  data$Rxn = Pathway_Matrix_Rxn$REACTION[match(data$KEGG, Pathway_Matrix_Rxn$KEGG)]
  
  
  RXN = data
  
  
  
  RXN <- RXN[complete.cases(RXN$KEGG),]
  RXN$Rxn = as.character(RXN$Rxn)
  RXN <- RXN[complete.cases(RXN$Rxn),]
  
  RXN[RXN == "NULL"] <- NA
  RXN[RXN == "NA"] <- NA
  RXN <- RXN[complete.cases(RXN$Rxn),]
  
  
  
  Rxn_Vect = RXN[, c("KEGG", "Rxn", "Val", "log2FoldChange")]
  Rxn_Vect$Rxn = gsub("\\c\\("," \\[",Rxn_Vect$Rxn)
  Rxn_Vect$Rxn = gsub("[[:punct:]]", "", Rxn_Vect$Rxn)
  Rxn_Vect$Rxn = gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",Rxn_Vect$Rxn))
  
  lst <- strsplit(as.character(Rxn_Vect$Rxn), ",")
  DF4 <- transform(Rxn_Vect[rep(1:nrow(Rxn_Vect), lengths(lst)),-5], Rxn= unlist(lst))
  
  return(DF4)
}

data_rxn <- Rxn_Path(data)



Rxn_to_Orthology <- function(data){
  
  DKEGG = as.list(data$KEGG)
  DRxn = as.list(data$Rxn)
  DVal = as.list(data$Val)
  DKEGG <- split(DKEGG,  ceiling(seq_along(DKEGG)/10))
  DRxn <- split(DRxn,  ceiling(seq_along(DRxn)/10))
  DVal <- split(DVal,  ceiling(seq_along(DVal)/10))
  K_Key = data
  
  KO <- llply(DRxn, function(x)keggGet(x))
  RX_KO <- as.data.frame(cbind(DRxn, KO))
  KO1 = unlist(KO, recursive = F)
  KO1 <- lapply(KO1, '[', c("ENTRY", "ORTHOLOGY"))
  
  n_obs_KO <- sapply(KO1, length)
  seq_max_KO <- seq_len(max(n_obs_KO))
  Pathway_Matrix_KO <- t(sapply(KO1, '[', i = seq_max_KO))
  
  TKO = as.data.frame(unlist(Pathway_Matrix_KO[,2], recursive = F))
  TKO = rownames_to_column(df = TKO, var = "Orthology_number")
  TKO$Orthology_number = gsub("^.*?.",".",TKO$Orthology_number)
  TKO$Orthology_number = gsub("[[:punct:]]", "", TKO$Orthology_number)
  TKO$Orthology_number = gsub("^.*?K","K",TKO$Orthology_number)
  colnames(TKO)[2] <- "TKO"
  
  
  Pathway_Matrix_KO = Pathway_Matrix_KO[,c(1,2)]
  colnames(Pathway_Matrix_KO)[1] <- "Rxn" 
  colnames(Pathway_Matrix_KO)[2] <- "KO"
  Pathway_Matrix_KO = as.data.frame(Pathway_Matrix_KO)
  
  Pathway_Matrix_KO$Rxn = unlist(Pathway_Matrix_KO$Rxn)
  Unlist_PKO <- data.frame(Rxn = rep(Pathway_Matrix_KO$Rxn, sapply(Pathway_Matrix_KO$KO, length)), KO = TKO$TKO) 
  
  Unlist_PKO$KEGG <- K_Key$KEGG[match(Unlist_PKO$Rxn, K_Key$Rxn)]
  Unlist_PKO$Val <- K_Key$Val[match(Unlist_PKO$Rxn, K_Key$Rxn)]
  Unlist_PKO$log2FoldChange <- K_Key$log2FoldChange[match(Unlist_PKO$Rxn, K_Key$Rxn)]
  D = Unlist_PKO
  D$Orthology_number <- TKO$Orthology_number[match(D$KO, TKO$TKO)]
  return(D)
}

data_ortho <- Rxn_to_Orthology(data_rxn)

Metabo_Path <- function(data, KEGG_key){
  data$Class = KEGG_key$Class[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_1 = KEGG_key$Subclass_1[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_2 = KEGG_key$Subclass_2[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_3 = KEGG_key$Subclass_3[match(data$KEGG, KEGG_key$KEGG)]
  data$Subclass_4 = KEGG_key$Subclass_4[match(data$KEGG, KEGG_key$KEGG)]
  data$Metabolite = OG_data$Metabolite[match(data$KEGG, OG_data$KEGG)]
  return(data)
  
}
data_ortho = Metabo_Path(data = data_ortho, KEGG_key = KEGG_key)


Orthology_hierarchy <- function(data, ortho_map){
  
  data$KO_Class <- ortho_map$KO_Class[match(data$Orthology_number, ortho_map$Orthology_number)]
  data$KO_Sub_class1 <- ortho_map$KO_Sub_class1[match(data$Orthology_number, ortho_map$Orthology_number)]
  data$KO_Sub_class2 <- ortho_map$KO_Sub_class2[match(data$Orthology_number, ortho_map$Orthology_number)]
  return(data)
}

data_ortho = Orthology_hierarchy(data = data_ortho, ortho_map = ortho_map)
Order = c("Metabolite","Val","log2FoldChange", "Class", "Subclass_1", "Subclass_2", "Subclass_3", 
          "Subclass_4", "Rxn", "KEGG", "Orthology_number", "KO", "KO_Class", "KO_Sub_class1", "KO_Sub_class2")
data_ortho = data_ortho[, Order]

return(data_ortho)
}

Orthology_Counts <- function(data, level){
  data_pos = data[data$Val %in% "Increase",]
  data_neg = data[data$Val %in% "Decrease",]
  
  data_pos_table = as.data.frame(table(data_pos[, level]))
  colnames(data_pos_table)[1] <- level
  colnames(data_pos_table)[2] <- "Significant_Changes"
  
  data_neg_table = as.data.frame(table(data_neg[,level]))
  colnames(data_neg_table)[1] <- level
  colnames(data_neg_table)[2] <- "Significant_Changes"
  
  data_pos_table = data_pos_table[apply(data_pos_table[2],1,function(z) !any(z==0)),]
  
  data_neg_table = data_neg_table[apply(data_neg_table[2],1,function(z) !any(z==0)),]
  
  data_neg_table[2] = data_neg_table[2] * -1
  
  
  data_counts = rbind(data_pos_table, data_neg_table)
  
  data_counts$Val = ifelse(data_counts$Significant_Changes < 0, "Decrease","Increase")
  
  return(data_counts)
}

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





