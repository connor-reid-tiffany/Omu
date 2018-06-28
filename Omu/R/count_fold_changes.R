#'count_fold_changes
#'
#' This function takes an input data frame that has been run in a statistical modeling function
#' and returns the number of compounds that significantly changed in each metabolite Class or Subclass.
#' @param count_data Output dataframe from the omu_summary function
#' @param ... Either a Class or Subclass column as a string, i.e. "Class
#' @param column The same value entered for the ... argument, i.e. column = "Class
#' @param sig_threshold Significance threshold for compounds that go towars the count,
#' sig_threshold = 0.05
#' @param keep_unknowns TRUE or FALSE for whether to drop compounds that weren't assigned
#' hierarchy metadata
#' @importFrom dplyr group_by_
#' @importFrom dplyr mutate
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#' t_test_df <- omu_summary(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE, p_adjust = "BH")
#' fold_change_counts <- count_fold_changes(count_data = t_test_df, "Class",
#' column = "Class", sig_threshold = 0.05, keep_unknowns = "FALSE")
#' @export

count_fold_changes <- function(count_data, ..., column, sig_threshold, keep_unknowns){
   log2FoldChange <- neg <- NULL

  count_data <- count_data[which(count_data[,"padj"] <= sig_threshold),]
  count_data <- count_data %>% group_by_(...) %>%
    mutate(Significant_Changes = sum(log2FoldChange>0),
           neg = sum(log2FoldChange<0))
  count_data <- count_data[,c(column, "Significant_Changes", "neg")]
  count_data$neg <- count_data$neg * -1
  output <- subset(count_data, select = c(1,3))
  colnames(output)[2] <- "Significant_Changes"
  count_data <- rbind(count_data, output)
  count_data <- subset(count_data, select = -neg)

  count_data$colour <- ifelse(count_data$Significant_Changes < 0, "Decrease","Increase")

  count_data <- count_data[apply(count_data[2],1,function(z) !any(z==0)),]


  if(keep_unknowns==TRUE){
    count_data = count_data
    }else if(keep_unknowns==FALSE){
     count_data <- drop_na(count_data)
    }
  unique(count_data[])

}
