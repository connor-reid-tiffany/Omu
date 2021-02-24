#' omu_summary
#' Performs comparison of means between two independent variables, standard deviation, standard error, FDR correction, fold change, log2FoldChange.
#' The order effects the fold change values
#' @param count_data should be a metabolomics count data frame
#' @param metadata is meta data
#' @param numerator is the variable you wish to compare against the denominator, in quotes
#' @param denominator see above, in quotes
#' @param response_variable the name of the column with your response variables
#' @param Factor the column name for your independent variables
#' @param log_transform TRUE or FALSE value for whether or not log transformation of data is performed
#' before the t test
#' @param p_adjust Method for adjusting the p value, i.e. "BH"
#' @param test_type One of "mwu", "students", or "welch" to determine which means comparison model to use
#' @importFrom dplyr filter
#' @importFrom plyr ldply
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_all
#' @importFrom magrittr %>%
#' @importFrom stats t.test
#' @importFrom stats wilcox.test
#' @importFrom stats p.adjust
#' @importFrom stats sd
#' @examples
#' \dontshow{c57_nos2KO_mouse_countDF <- c57_nos2KO_mouse_countDF[1:12,]}
#' omu_summary(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#' numerator = "Strep", denominator = "Mock", response_variable = "Metabolite", Factor = "Treatment",
#' log_transform = TRUE, p_adjust = "BH", test_type = "welch")
#' @export

omu_summary <- function(count_data, metadata, numerator, denominator, response_variable,
                        Factor, log_transform, p_adjust, test_type){

if(identical(as.character(colnames(test_h)[unlist(lapply(test_h, is.numeric))]), as.character(metadata$Sample))==FALSE){

  stop("Sample names in count_data and metadata do not match.")

}

if(any(colnames(metadata)=="Sample")==FALSE){

  stop("metadata is missing Sample column")

}

if (log_transform==TRUE){

find_zeros <- function(x){

x2 <- sapply(x, is.numeric)

x <- x[,x2]

xl <- sapply(x, function(x) any(x==0))

}

if (any(find_zeros(count_data)==TRUE)){

stop("Your data have zero values. If you trust these zeros are legitimate, set log_transform to FALSE and consider
using the square root to center your data.")

}

}

if (nrow(check_zeros(count_data = count_data, metadata = metadata, Factor = Factor,
                     numerator = numerator, denominator = denominator, response_variable = response_variable)) > 0) {

  warning("There are zero values in at least 25 percent of your samples within at least one of your Factor
          levels for these metabolites. Consider using check_zeros to subset your data.")

  print(unique(check_zeros(count_data = count_data, metadata = metadata, Factor = Factor,
                           numerator = numerator, denominator = denominator, response_variable = response_variable)[,1]))
}

  #Temporarily separate meta data from counts and store in other object
  rownames(count_data) <- count_data[,response_variable]
  count_data[,response_variable] <- NULL
  data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]


  #Transform for 'normalization' and T test
  data_Transpose <- as.data.frame(t(data_Int))
  data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), data_Transpose))
  Factor = metadata[, Factor]
  data_Transpose$Factor = Factor[match(metadata$Sample, data_Transpose$Sample)]
  data_Subset <- filter(data_Transpose, Factor==numerator|Factor==denominator)
  rownames(data_Subset) <- data_Subset[,"Sample"]
  data_Subset[,"Sample"] <- NULL

  #Treatment_vect <- as.data.frame(data_Subset$Treatment)
  data_Numeric <- data_Subset[sapply(data_Subset, function(x) is.numeric(x))]
  data_Numeric <- data.frame(lapply(data_Numeric,
                                    function(x) as.numeric(as.character(x))),check.names=F, row.names = rownames(data_Numeric))

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
  if(test_type == "students"){
  Run_Tests <- function(data_mod, Vect, model) {
    results <- ldply(
      Vect, function(Metabolite) {
        t_val = t.test(data_mod[[Metabolite]] ~ model)$statistic
        p_val = t.test(data_mod[[Metabolite]] ~ model)$p.value
        return(data.frame(Metabolite=Metabolite, t_value=t_val, pval = p_val))
      })
  }
  }else if(test_type == "mwu"){
    Run_Tests <- function(data_Log, Vect, model) {
      results <- ldply(
        Vect, function(Metabolite) {
          t_val = wilcox.test(data_mod[[Metabolite]] ~ model)$statistic
          p_val = wilcox.test(data_mod[[Metabolite]] ~ model)$p.value
          return(data.frame(Metabolite=Metabolite, t_value=t_val, pval = p_val))
        })
    }
  }else if(test_type == "welch"){
    Run_Tests <- function(data_Log, Vect, model) {
      results <- ldply(
        Vect, function(Metabolite) {
          t_val = t.test(data_mod[[Metabolite]] ~ model, var.equal = FALSE)$statistic
          p_val = t.test(data_mod[[Metabolite]] ~ model, var.equal = FALSE)$p.value
          return(data.frame(Metabolite=Metabolite, t_value=t_val, pval = p_val))
        })
  }
  }
  results <- Run_Tests(data_Log = data_Log, Vect = Vect, model = model)

  #Compute raw count means
  data_Log$Factor <- factor(data_Log$Factor, levels = c(numerator, denominator))
  data_Numeric$Factor = data_Log$Factor
  Means <- data_Numeric %>% group_by(Factor) %>% summarise_all(mean)
  Means_T = as.data.frame(t(Means))
  colnames(Means_T) <- as.character(unlist(Means_T[1,]))
  Means_T = Means_T[-1,]
  colnames(Means_T)[1] <- numerator
  colnames(Means_T)[1] <- paste(colnames(Means_T)[1], "mean", sep = ".")
  colnames(Means_T)[2] <- denominator
  colnames(Means_T)[2] <- paste(colnames(Means_T)[2], "mean", sep = ".")

  #Compute log2FoldChange
  Means_T[,3] <- as.numeric(as.character(Means_T[,1])) / as.numeric(as.character(Means_T[,2]))
  colnames(Means_T)[3] <- "Fold_Change"
  Means_T[,4] = log2(as.numeric(as.character(Means_T[,3])))
  colnames(Means_T)[4] = "log2FoldChange"
  Means_T <- cbind(rownames(Means_T), data.frame(Means_T, row.names=NULL))
  colnames(Means_T)[1] <- response_variable

  #Create standard error function
  st.err <- function(x) sd(x)/sqrt(length(x))

  #Compute standard deviation
  stdev <- data_Log %>% group_by(Factor) %>% summarise_all(sd)
  stdev <- as.data.frame(stdev)
  rownames(stdev) <- stdev[,"Factor"]
  stdev[,"Factor"] <- NULL
  stdev_t = as.data.frame(t(stdev))
  stdev_t <- cbind(rownames(stdev_t), data.frame(stdev_t, row.names=NULL))
  colnames(stdev_t)[1] <- response_variable
  colnames(stdev_t)[2] <- numerator
  colnames(stdev_t)[2] <- paste(colnames(stdev_t)[2], "stdev", sep = ".")
  colnames(stdev_t)[3] <- denominator
  colnames(stdev_t)[3] <- paste(colnames(stdev_t)[3], "stdev", sep = ".")

  #Compute standard error
  SE <- data_Log %>% group_by(Factor) %>% summarise_all(st.err)
  SE <- as.data.frame(SE)
  rownames(SE) <- SE[,"Factor"]
  SE[,"Factor"] <- NULL
  SE_t = as.data.frame(t(SE))
  SE_t <- cbind(rownames(SE_t), data.frame(SE_t, row.names=NULL))
  colnames(SE_t)[1] <- response_variable
  colnames(SE_t)[2] <- numerator
  colnames(SE_t)[2] <- paste(colnames(SE_t)[2], "std.err", sep = ".")
  colnames(SE_t)[3] <- denominator
  colnames(SE_t)[3] <- paste(colnames(SE_t)[3], "std.err", sep = ".")

  #Merge metabo_path metadata, means & fold change, and t test results by metabolite name
  colnames(results)[1] <- response_variable
  results$padj = p.adjust(results$pval, method = p_adjust)
  count_data = cbind(rownames(count_data), data.frame(count_data, row.names=NULL))
  colnames(count_data)[1] <- "Metabolite"
  results = merge(results, count_data, by = response_variable, all = TRUE)
  results = merge(Means_T, results, by = response_variable, all.y = TRUE)
  results = merge(SE_t, results, by = response_variable, all.y = TRUE)
  results = merge(stdev_t, results, by = response_variable, all.y = TRUE)

  class(results) = append(class(results), "cpd")

  return(results)
}
