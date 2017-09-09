#Function for assigning genes associated with functional orthologies
#data is a data table that has undergone KO mapping
#Organism map is an organism hierarchy file that must be downloaded to 
#assign an Organism hierarchy to the genes



Map_Genes <- function(data, Organism_map){

#Create a nested list of Orthology numbers that can be sent to the KEGG API
Get_Genes <- function(data){
Vector <- as.data.frame(data$Orthology_number)
colnames(Vector)[1] <- "Orthology_number"
Vector2 = as.vector(Vector[ grep("nan", Vector$Orthology_number, invert = TRUE) , ])
Vector_Split  <- split(Vector2,  ceiling(seq_along(Vector2)/10))
Vector_Split2 = llply(Vector_Split, as.list)

#ASCII art "patience snail" to amuse user. This step can take up to 20 minutes if data wasn't subsetted previously 
Text_Art <-  c(" / /", " L_L_", "/    \\", "|00  |       _______", "|_/  |      /  ___  \\", 
               "|    |     /  /   \\  \\", "|    |_____\\  \\_  /  /", " \\          \\____/  /_____", 
               "  \\ _______________/______\\................please be patient =)"
)
cat(Text_Art, sep = "\n")
#Send nested list to KEGG API, list is in groups of ten due to server side limitations
KO_to_Genes <- llply(Vector_Split2, function(x)keggGet(x), .progress = "text")
KO_to_Genes_unlist <- unlist(KO_to_Genes, recursive = F)

#Extract desired information from new nested list received from KEGG API
KO_Genes_DF <- lapply(KO_to_Genes_unlist, '[', c("NAME", "ENTRY","DEFINITION", "GENES"))

#Reformat list into a Data frame, then add gene names to original data frame by matching by Orthology number                    
n_obs <- sapply(KO_Genes_DF, length)
seq_max <- seq_len(max(n_obs))
KO_Genes_Mat <- t(sapply(KO_Genes_DF, '[', i = seq_max))
colnames(KO_Genes_Mat)[2] <- "Orthology_number"
KO_Genes_Mat = as.data.frame(KO_Genes_Mat)
data$Genes = KO_Genes_Mat$GENES[match(data$Orthology_number, KO_Genes_Mat$Orthology_number)]
data$GeneOperon = KO_Genes_Mat$NAME[match(data$Orthology_number, KO_Genes_Mat$Orthology_number)]

#Rename for convenience
Genes_Vect = data
#Genes_subset <- c("Orthology_number", "Genes", "KEGG", "Metabolite", "Val", "log2FoldChange", "Class", "Subclass_1", "Subclass_2", 
                  #"Subclass_3", "Subclass_4", "KO", "KO_Class", "KO_Sub_class1", "KO_Sub_class2", "GeneOperon")
#Genes_Vect = Genes[, Genes_subset]
#Clean up genes column using REGEX
Genes_Vect$Genes = gsub("\\c\\("," \\[",Genes_Vect$Genes)
Genes_Vect$Genes = gsub("\\[|\\]", "", Genes_Vect$Genes)

#Genes column is still in the form of a character vector, needs to be reduced to one gene per row
#To accomplish this, all other rows will be replicated for each gene in each character vector until there is one gene per row
#Also cleans up data with another REGEX in order to easily assign organism hierarchies
lst <- strsplit(as.character(Genes_Vect$Genes), ",")
DF <- transform(Genes_Vect[rep(1:nrow(Genes_Vect), lengths(lst)),-1], Genes= unlist(lst))
DF <- as.data.frame(sapply(DF, function(x) gsub("\"", "", x)))
DF$Gene_number <- DF$Genes
DF2 <- str_split_fixed(DF$Genes, ":", 2)
DF2 = as.data.frame(DF2)
colnames(DF2)[1] <- "Org" 

#Merge split data frames after making it one gene per row, rename Column for convenience and remove duplicate Gene column
#Change gene names to lowercase for assigning organism hierarchy                           
DF_Complete <- data.frame(do.call('rbind', strsplit(as.character(DF2$Org),':',fixed=TRUE)))
DF_Complete = as.data.frame(cbind(DF, DF_Complete))
colnames(DF_Complete)[17] <- "Org"
DF_Complete$Org <- tolower(DF_Complete$Org)
DF_Complete$Org <- trimws(DF_Complete$Org)
DF_Complete = DF_Complete[ , -which(names(DF_Complete) %in% c("Gene_number"))]

  }

#Assign Organism hierarchy. KEGG BRITE db is a little complicated due to bacterial strains, serovars, and different clinical isolates                           
data = Get_Genes(data)
Organism_Hierarchy <- function(data, Organism_map){
  DF <- left_join(data, Organism_map, by = "KEGG)
  return(DF)
}

data = Organism_Hierarchy(data, Organism_map)  
return(data)
}
