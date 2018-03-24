#'pie_chart
#'Makes pie chart as ggplot2 object from ra_table function output
#'@param data a dataframe object of percents. output from ra_table function
#'@param variable The meta data variable you are measuring, i.e. "Class"
#'@param column either "Increase", "Decrease", or "Significant_Changes"
#'@param color string denoting color for outline. use NA for no outline
#'@example pie_chart(data = ra_table, variable = "Increase",column = "Class", color = "black")
#'@export

pie_chart <- function(data,variable, column, color){
  variable <- reorder(variable, column)
  bar<- ggplot(data)+
    geom_bar(width = 1,aes(x="", y=data[,column], fill=data[,variable]),
             stat = "identity", color = color)
  pie <- bar + coord_polar("y", start=0) +
    theme_bw() + theme(panel.border = element_blank()) +
    theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) +
    theme(axis.title = element_blank()) +
    theme(axis.ticks = element_blank())
}
