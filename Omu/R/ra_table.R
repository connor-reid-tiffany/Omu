#'ra_table
#'
#'Make a relative abundance table from the count_fold_changes function output. For making pie charts
#'@param data data frame object from the count_fold_changes function
#'@param variable meta data from count_fold_changes, i.e. "Class"
#'@export
#'@example ra_table(data = count_fold_changes_output, variable = "Class")


ra_table <- function(data,variable){

    #Make decrease table
    data_dec = data[data$colour %in% "Decrease",]
    data_dec$Decrease = abs(data_dec$Significant_Changes)
    data_dec$Decrease = prop.table(data_dec$Decrease)
    data_dec$Decrease = data_dec$Decrease * 100
    data_dec = data_dec[,c(1,4)]

    #Make increase table
    data_inc = data[data$colour %in% "Increase",]
    data_inc$Increase = data_inc$Significant_Changes
    data_inc$Increase = prop.table(data_inc$Increase)
    data_inc$Increase = data_inc$Increase * 100
    data_inc = data_inc[,c(1,4)]

    #Make total table
    data_total = data
    data_total$Significant_Changes = abs(data_total$Significant_Changes)
    data_total = ddply(data_total, variable, numcolwise(sum))
    data_total$Significant_Changes = prop.table(data_total$Significant_Changes)
    data_total$Significant_Changes = data_total$Significant_Changes * 100

    #Merge tables
    data_join <- left_join(data_total, data_dec, variable)
    data_join = left_join(data_join, data_inc, variable)
}
