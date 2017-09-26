#'Count Fold Changes
#'
#'This function takes an input data frame that has been run in a statistical modeling function 
#'and returns the number of compounds that significantly changed in each metabolite Class or Subclass.
#'@param data Output dataframe from the Clean_DESeq_results function
#'@param ... Either a Class or Subclass column listed in paretheses, i.e. "Class
#'@param column The same value entered for the ... parameter, i.e. column = "Class
#'@param alpha Significance threshold, i.e. alpha = 0.05
#'@keywords metabo
#'@export
#'@examples Count_Fold_Changes(data = data, "Class", column = "Class", alpha = 0.05)
#'Count_Fold_Changes()
Count_Fold_Changes <- function(data, ..., column, alpha){
  data = data[which(data$padj < alpha), ]
  data = data %>% group_by_(...) %>%
    mutate(Significant_Changes = sum(log2FoldChange>0),
           neg = sum(log2FoldChange<0))
  data = data[,c(column, "Significant_Changes", "neg")]
  data$neg <- data$neg * -1
  output <- subset(data, select = c(1,3))
  colnames(output)[2] <- "Significant_Changes"
  data <- rbind(data, output)
  data = subset(data, select = -neg)
  
  data$colour = ifelse(data$Significant_Changes < 0, "Decrease","Increase")
  
  data = data[apply(data[2],1,function(z) !any(z==0)),]
  unique(data[])
}



