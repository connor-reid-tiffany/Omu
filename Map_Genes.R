Map_Genes <- function(data, Organism_map){

Get_Genes <- function(data){
Vector <- as.data.frame(data$Orthology_number)
colnames(Vector)[1] <- "Orthology_number"
Vector2 = as.vector(Vector[ grep("nan", Vector$Orthology_number, invert = TRUE) , ])
Vector_Split  <- split(Vector2,  ceiling(seq_along(Vector2)/10))
Vector_Split2 = llply(Vector_Split, as.list)

Text_Art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\", 
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____", 
               "  \\ _______________/______\\................please be patient =)"
)
cat(Text_Art, sep = "\n")
KO_to_Genes <- llply(Vector_Split2, function(x)keggGet(x), .progress = "text")
KO_to_Genes_unlist <- unlist(KO_to_Genes, recursive = F)

KO_Genes_DF <- lapply(KO_to_Genes_unlist, '[', c("NAME", "ENTRY","DEFINITION", "GENES"))

n_obs <- sapply(KO_Genes_DF, length)
seq_max <- seq_len(max(n_obs))
KO_Genes_Mat <- t(sapply(KO_Genes_DF, '[', i = seq_max))
colnames(KO_Genes_Mat)[2] <- "Orthology_number"
KO_Genes_Mat = as.data.frame(KO_Genes_Mat)
data$Genes = KO_Genes_Mat$GENES[match(data$Orthology_number, KO_Genes_Mat$Orthology_number)]
data$GeneOperon = KO_Genes_Mat$NAME[match(data$Orthology_number, KO_Genes_Mat$Orthology_number)]

Genes_Vect = data

#Genes_subset <- c("Orthology_number", "Genes", "KEGG", "Metabolite", "Val", "log2FoldChange", "Class", "Subclass_1", "Subclass_2", 
                  #"Subclass_3", "Subclass_4", "KO", "KO_Class", "KO_Sub_class1", "KO_Sub_class2", "GeneOperon")
#Genes_Vect = Genes[, Genes_subset]
Genes_Vect$Genes = gsub("\\c\\("," \\[",Genes_Vect$Genes)
Genes_Vect$Genes = gsub("\\[|\\]", "", Genes_Vect$Genes)


lst <- strsplit(as.character(Genes_Vect$Genes), ",")
DF <- transform(Genes_Vect[rep(1:nrow(Genes_Vect), lengths(lst)),-1], Genes= unlist(lst))
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

  }

data = Get_Genes(data)
Organism_Hierarchy <- function(data, Organism_map){
  data$Org = as.factor(data$Org)
  data$Domain <- Organism_map$Domain[match(data$Org, Organism_map$Org)]
  data$Kingdom <- Organism_map$Kingdom[match(data$Org, Organism_map$Org)]
  data$Group <- Organism_map$Group[match(data$Org, Organism_map$Org)]
  data$Subgroup_1 <- Organism_map$Subgroup_1[match(data$Org, Organism_map$Org)]
  data$Genus <- Organism_map$Genus[match(data$Org, Organism_map$Org)]
  data$Species <- Organism_map$Species[match(data$Org, Organism_map$Org)]
  data$Meta1 <- Organism_map$Meta1[match(data$Org, Organism_map$Org)]
  data$Meta2 <- Organism_map$Meta2[match(data$Org, Organism_map$Org)]
  data$Meta3 <- Organism_map$Meta3[match(data$Org, Organism_map$Org)]
  data$Meta4 <- Organism_map$Meta4[match(data$Org, Organism_map$Org)]
  data$Meta5 <- Organism_map$Meta5[match(data$Org, Organism_map$Org)]
  data$Meta6 <- Organism_map$Meta6[match(data$Org, Organism_map$Org)]
  return(data)
}

data = Organism_Hierarchy(data, Organism_map)  
return(data)
}
