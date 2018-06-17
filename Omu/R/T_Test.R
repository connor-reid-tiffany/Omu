#' t_test
#' Performs t test, Standard Error, FDR correction, Fold Change, log2FoldChange. The order effects the fold change values
#' @param data should be a metabolomics count data frame
#' @param colData is meta data
#' @param numerator is the variable you wish to compare against the denominator, in quotes
#' @param denominator see above, in quotes
#' @param response_variable the name of the column with your response variables
#' @param Factor the column name for your independent variables
#' @param log_transform TRUE or FALSE value for whether or not log transformation of data is performed before the t test
#' @importFrom dplyr filter
#' @importFrom plyr ldply
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_all
#' @importFrom dplyr funs
#' @importFrom plotrix std.error
#' @importFrom magrittr %>%
#' @importFrom stats t.test
#' @importFrom stats p.adjust
#' @examples
#' t_test(data = c57_nos2KO_mouse_countDF, colData = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE)
#' @export

t_test <- function(data, colData, numerator, denominator, response_variable, Factor, log_transform){


  #Temporarily separate meta data from counts and store in other object
  rownames(data) <- data[,response_variable]
  data[,response_variable] <- NULL
  data_Int <- data[sapply(data, function(x) is.numeric(x))]


  #Transform for 'normalization' and T test
  data_Transpose <- as.data.frame(t(data_Int))
  data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), data_Transpose))
  Factor = colData[, Factor]
  data_Transpose$Factor = Factor[match(colData$Sample, data_Transpose$Sample)]
  data_Subset <- filter(data_Transpose, Factor==numerator|Factor==denominator)
  rownames(data_Subset) <- data_Subset[,"Sample"]
  data_Subset[,"Sample"] <- NULL

  #Treatment_vect <- as.data.frame(data_Subset$Treatment)
  data_Numeric <- data_Subset[sapply(data_Subset, function(x) is.numeric(x))]
  data_Numeric <- data.frame(lapply(data_Numeric, function(x) as.numeric(as.character(x))),check.names=F, row.names = rownames(data_Numeric))

  #Normalize the data using natural logarithm
  data_Log <- as.data.frame(log(data_Numeric))
  data_Log$Factor = data_Subset$Factor
  cols_to_test <- data_Log[sapply(data_Log, function(x) is.numeric(x))]
  Vect = colnames(cols_to_test)
  data_Numeric$Factor <- data_Subset$Factor

  #Create arguments
  if(log_transform==FALSE){
    data_mod = data_Numeric
  } else if (log_transform==TRUE){
    data_mod = data_Log
  }

  #Create arguments for T Test function
  model = data_mod[, "Factor"]

  #T Test function in function envir. Iterates T_Test across all response variables
  Run_T_Tests <- function(data_Log, Vect, model) {
    results <- ldply(
      Vect, function(Metabolite) {
        t_val = t.test(data_mod[[Metabolite]] ~ model)$statistic
        p_val = t.test(data_mod[[Metabolite]] ~ model)$p.value
        return(data.frame(Metabolite=Metabolite, t_value=t_val, pval = p_val))
      })
  }

  results <- Run_T_Tests(data_Log = data_Log, Vect = Vect, model = model)

  #Compute raw count means
  data_Log$Factor <- factor(data_Log$Factor, levels = c(numerator, denominator))
  data_Numeric$Factor = data_Log$Factor
  Means <- data_Numeric %>% group_by(Factor) %>% summarise_all(funs(mean))
  Means_T = as.data.frame(t(Means))
  colnames(Means_T) <- as.character(unlist(Means_T[1,]))
  Means_T = Means_T[-1,]
  colnames(Means_T)[1] <- "Numerator_Mean"
  colnames(Means_T)[2] <- "Denominator_Mean"

  #Compute log2FoldChange
  Means_T[,3] <- as.numeric(as.character(Means_T[,1])) / as.numeric(as.character(Means_T[,2]))
  colnames(Means_T)[3] <- "Fold_Change"
  Means_T[,4] = log2(as.numeric(as.character(Means_T[,3])))
  colnames(Means_T)[4] = "log2FoldChange"
  Means_T <- cbind(rownames(Means_T), data.frame(Means_T, row.names=NULL))
  colnames(Means_T)[1] <- response_variable

  #Compute standard error
  SE <- data_Log %>% group_by(Factor) %>% summarise_all(funs(std.error))
  SE <- as.data.frame(SE)
  rownames(SE) <- SE[,"Factor"]
  SE[,"Factor"] <- NULL
  SE_t = as.data.frame(t(SE))
  SE_t <- cbind(rownames(SE_t), data.frame(SE_t, row.names=NULL))
  colnames(SE_t)[1] <- response_variable

  #Merge metabo_path metadata, means & fold change, and t test results by metabolite name
  colnames(results)[1] <- response_variable
  results$padj = p.adjust(results$pval, method = "BH")
  data = cbind(rownames(data), data.frame(data, row.names=NULL))
  colnames(data)[1] <- "Metabolite"
  results = merge(results, data, by = response_variable, all = TRUE)
  results = merge(Means_T, results, by = response_variable, all.y = TRUE)
  results = merge(SE_t, results, by = response_variable, all.y = TRUE)

  class(results) = append(class(results), "cpd")

  return(results)
}
