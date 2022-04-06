#' Get counts for significant fold changes by metabolite class.
#'
#' @description Takes an input data frame from the output of omu_summary and creates a
#' data frame of counts for significantly changed metabolites by class hierarchy data.
#' @param count_data Output dataframe from the omu_summary function or omu_anova.
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
#' \dontshow{c57_nos2KO_mouse_countDF <- c57_nos2KO_mouse_countDF[1:30,]}
#' t_test_df <- omu_summary(count_data = c57_nos2KO_mouse_countDF,
#' metadata = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite",
#' Factor = "Treatment", log_transform = TRUE, p_adjust = "BH", test_type = "welch")
#'
#' fold_change_counts <- count_fold_changes(count_data = t_test_df,
#' column = "Class", sig_threshold = 0.05, keep_unknowns = "FALSE")
#' @export

count_fold_changes <- function(count_data, column, sig_threshold, keep_unknowns){

  if(is.null(count_data$padj)==TRUE){

    stop("count_data must be the output of omu_summary or omu_anova and have a padj column
    and a log2FoldChange column")

  }

  if(any(names(count_data) %in% column)==FALSE){

    stop("count_data is missing metabolite metadata. Did you forget to use assign_hierarchy?")

  }


class(count_data) <- "data.frame"

   log2FoldChange <- neg <- NULL

  count_data <- count_data[which(count_data[,"padj"] <= sig_threshold),]

  count_data <- count_data %>% group_by(count_data[,column]) %>%
   mutate(Significant_Changes = sum(log2FoldChange>0),
          neg = sum(log2FoldChange<0))

 count_data <- count_data[,c(column, "Significant_Changes", "neg")]
 count_data$neg <- count_data$neg * -1
 output <- count_data[,c(1,3)]
 count_data <- count_data[,c(1,2)]
 colnames(output)[2] <- "Significant_Changes"
 count_data <- rbind(count_data, output)

 count_data$color <- ifelse(count_data$Significant_Changes < 0, "Decrease","Increase")

 count_data <- count_data[apply(count_data[2],1,function(z) !any(z==0)),]


  if(keep_unknowns==TRUE){
    count_data = count_data
    }else if(keep_unknowns==FALSE){
     count_data <- drop_na(count_data)
    }
  unique(count_data[])

}
