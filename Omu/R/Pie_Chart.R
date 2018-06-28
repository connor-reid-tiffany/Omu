#' pie_chart
#' Makes pie chart as ggplot2 object from ra_table function output
#' @param ratio_data a dataframe object of percents. output from ra_table function
#' @param variable The metadata variable you are measuring, i.e. "Class"
#' @param column either "Increase", "Decrease", or "Significant_Changes"
#' @param color string denoting color for outline. use NA for no outline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 coord_polar
#' @importFrom stats reorder
#' @examples
#' c57_nos2KO_mouse_countDF <- assign_hierarchy(c57_nos2KO_mouse_countDF, TRUE, "KEGG")
#'
#' t_test_df <- omu_summary(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE, p_adjust = "BH")
#'
#' fold_change_counts <- count_fold_changes(count_data = t_test_df, "Class",
#' column = "Class", sig_threshold = 0.05, keep_unknowns = FALSE)
#'
#' ra_table <- ra_table(fc_data = fold_change_counts, variable = "Class")
#'
#' pie_chart(ratio_data = ra_table, variable = "Class", column = "Decrease", color = "black")
#' @export

pie_chart <- function(ratio_data, variable, column, color){

  variable <- reorder(variable, column)
  bar<- ggplot(ratio_data)+
    geom_bar(width = 1,aes(x="", y=ratio_data[,column], fill=ratio_data[,variable]),
             stat = "identity", color = color)
  pie <- bar + coord_polar("y", start=0) +
    theme_bw() + theme(panel.border = element_blank()) +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.ticks = element_blank())
}
