#' Create a ratio table
#'
#' @title Creates a ratio table from the count_fold_changes function output.
#' @param fc_data data frame output from the count_fold_changes function
#' @param variable metadata from count_fold_changes, i.e. "Class"
#' @importFrom dplyr left_join
#' @importFrom plyr ddply
#' @importFrom plyr numcolwise
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#' \dontshow{c57_nos2KO_mouse_countDF <- c57_nos2KO_mouse_countDF[1:20,]}
#' t_test_df <- omu_summary(count_data = c57_nos2KO_mouse_countDF,
#' metadata = c57_nos2KO_mouse_metadata, numerator = "Strep", denominator = "Mock",
#' response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE, p_adjust = "BH", test_type = "welch")
#'
#' fold_change_counts <- count_fold_changes(count_data = t_test_df, "Class",
#' column = "Class", sig_threshold = 0.05, keep_unknowns = FALSE)
#'
#' ra_table(fc_data = fold_change_counts, variable = "Class")
#' @export

ra_table <- function(fc_data,variable){

    #Make decrease table
    data_dec = fc_data[fc_data$color %in% "Decrease",]
    data_dec$Decrease = abs(data_dec$Significant_Changes)
    data_dec$Decrease = prop.table(data_dec$Decrease)
    data_dec$Decrease = data_dec$Decrease * 100
    data_dec = data_dec[,c(1,4)]

    #Make increase table
    data_inc = fc_data[fc_data$color %in% "Increase",]
    data_inc$Increase = data_inc$Significant_Changes
    data_inc$Increase = prop.table(data_inc$Increase)
    data_inc$Increase = data_inc$Increase * 100
    data_inc = data_inc[,c(1,4)]

    #Make total table
    data_total = fc_data
    data_total$Significant_Changes = abs(data_total$Significant_Changes)
    data_total = ddply(data_total, variable, numcolwise(sum))
    data_total$Significant_Changes = prop.table(data_total$Significant_Changes)
    data_total$Significant_Changes = data_total$Significant_Changes * 100

    #Merge tables
    data_join <- left_join(data_total, data_dec, variable)
    data_join = left_join(data_join, data_inc, variable)

    return(data_join)
}
