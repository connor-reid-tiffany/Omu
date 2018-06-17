#'count_fold_changes
#'
#' This function takes an input data frame that has been run in a statistical modeling function
#' and returns the number of compounds that significantly changed in each metabolite Class or Subclass.
#' @param data Output dataframe from the Clean_DESeq_results function
#' @param ... Either a Class or Subclass column listed in paretheses, i.e. "Class
#' @param column The same value entered for the ... parameter, i.e. column = "Class
#' @param sig_threshold Significance threshold, i.e. sig_threshold = 0.05
#' @importFrom dplyr group_by_
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#' t_test_df <- t_test(data = c57_nos2KO_mouse_countDF, colData = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE)
#' fold_change_counts <- count_fold_changes(data = t_test_df, "Class",
#' column = "Class", sig_threshold = 0.05)
#' @export

count_fold_changes <- function(data, ..., column, sig_threshold){
  data <- data[which(data[,"padj"] < sig_threshold),]
  data <- data %>% group_by_(...) %>%
    mutate(Significant_Changes = sum(log2FoldChange>0),
           neg = sum(log2FoldChange<0))
  data <- data[,c(column, "Significant_Changes", "neg")]
  data$neg <- data$neg * -1
  output <- subset(data, select = c(1,3))
  colnames(output)[2] <- "Significant_Changes"
  data <- rbind(data, output)
  data <- subset(data, select = -neg)

  data$colour <- ifelse(data$Significant_Changes < 0, "Decrease","Increase")

  data <- data[apply(data[2],1,function(z) !any(z==0)),]
  unique(data[])
}
