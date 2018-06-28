#'omu_anova
#'Function to apply an anova across all response variables. Options for running the model
#'on 1 or 2 independent variables individually, or on the interaction between two
#'independent variables
#'@param count_data A metabolomics count data frame
#'@param metadata Metadata dataframe for the metabolomics count data frame
#'@param response_variable String of the column header for the response variables,
#'usually "Metabolite"
#'@param var1 String of the first independent variable you wish to test
#'@param var2 String of the second independent variable you wish to test. Optional parameter
#'@param interaction Boolean of TRUE or FALSE for whether or not you wish to model
#'an interaction between independent variables. Optional parameter
#'@param log_transform Boolean of TRUE or FALSE for whether or not you wish to log transform
#'your metabolite counts
#'@param p_adjust Method for p value adjustment, i.e. "BH"
#'@importFrom dplyr left_join
#'@importFrom plyr llply
#'@importFrom stats p.adjust
#'@importFrom stats anova
#'@importFrom stats lm
#'@examples
#'anova_df <- omu_anova(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#'response_variable = "Metabolite", var1 = "Treatment", var2 = "Background", log_transfor = TRUE,
#'p_adjust = "BH", interaction = TRUE)
#'@export


omu_anova <- function(count_data, metadata, response_variable, var1, var2, interaction, log_transform, p_adjust){

  #Make variables for setting column names
  variable1 = var1
  variable2 = var2

  #Convert column 1 to rownames
  rownames(count_data) <- count_data[,response_variable]
  count_data[,response_variable] <- NULL

  #Convert classes to numeric
  data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]

  #Transform for 'normalization'
  data_Transpose <- as.data.frame(t(data_Int))
  data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), data_Transpose))
  data_Transpose = as.data.frame(cbind(data_Transpose, metadata))
  nums <- sapply(data_Transpose, is.integer)
  factors <- sapply(data_Transpose, is.factor)

  #Convert class to numeric, natural log transform, and remerge num columns with factor columns
  #by sample name
  data_Num <- data.frame(lapply(data_Transpose[ , nums], function(x) as.numeric(as.integer(x))),
  check.names=F, row.names = rownames(data_Transpose))
  Vect = colnames(data_Num)
  data_Ln <- log(data_Num)
  data_Ln <- as.data.frame(cbind(Sample = rownames(data_Ln), data_Ln))
  data_Fact <- data_Transpose[, factors]
  data_Ln = merge(data_Ln, data_Fact, by = 'Sample')
  data_Ln <- data_Ln[, !names(data_Ln) %in% c('Sample', 'Sample.1')]
  data_Num <- as.data.frame(cbind(Sample = rownames(data_Num), data_Num))
  data_Num = merge(data_Num, data_Fact, by = 'Sample')
  data_Num <- data_Num[, !names(data_Num) %in% c('Sample', 'Sample.1')]

  #Create arguments
  if(log_transform==FALSE){
    data_mod = data_Num
  } else if (log_transform==TRUE){
    data_mod = data_Ln
  }

  Mod = data_Fact
  Mod = Mod[, !names(Mod) %in% c('Sample', "Sample.1")]

  if(missing(var2) & interaction==FALSE){
    var1 = metadata[,var1]
    results <- llply(
      Vect, function(x) {
        models <- lm(data_mod[[x]] ~ var1, data = Mod)
      })
    names(results) <- Vect
    results <- lapply(results, anova)

    results <- sapply(results, cbind)
    results <- t(results)
    results <- as.data.frame(results[,"Pr(>F)"])
    results <- as.data.frame(t(results))
    colnames(results)[1] <- variable1
    colnames(results)[1] <- paste(colnames(results)[1], "pval", sep = ".")

    #Compute adjusted p value
    results$padj = p.adjust(results[,1], method = p_adjust)
    colnames(results)[3] <- paste(colnames(results)[3], variable1, sep = ".")

    #Merge metadata, means & fold change, and t test results by metabolite name
    count_data <- cbind(rownames(count_data), data.frame(count_data, row.names=NULL))
    colnames(count_data)[1] <- response_variable
    results <- cbind(rownames(results), data.frame(results, row.names=NULL))
    colnames(results)[1] <- response_variable
    results = left_join(results, count_data, by = response_variable)
    results = results[, !(colnames(results) %in% c("V2"))]
    return(results)

  } else if(interaction==FALSE){
    var1 = metadata[, var1]
    var2 = metadata[, var2]
    results <- llply(
      Vect, function(x) {
        models <- lm(data_mod[[x]] ~ var1 + var2, data = Mod)
      })
    names(results) <- Vect
    results <- lapply(results, anova)

    results <- sapply(results, cbind)
    results <- t(results)
    results <- as.data.frame(results[,"Pr(>F)"])
    results <- as.data.frame(t(results))
    results <- results[,1:2]
    colnames(results)[1] <- variable1
    colnames(results)[1] <- paste(colnames(results)[1], "pval", sep = ".")
    colnames(results)[2] <- variable2
    colnames(results)[2] <- paste(colnames(results)[2], "pval", sep = ".")

    #Compute adjusted p value
    results$padj = p.adjust(results[,1], method = p_adjust)
    colnames(results)[3] <- paste(colnames(results)[3], variable1, sep = ".")
    results$padj = p.adjust(results[,2], method = p_adjust)
    colnames(results)[4] <- paste(colnames(results)[4], variable2, sep = ".")

    #Merge metadata, means & fold change, and t test results by metabolite name
    results <- cbind(rownames(results), data.frame(results, row.names=NULL))
    colnames(results)[1] <- response_variable
    count_data <- cbind(rownames(count_data), data.frame(count_data, row.names=NULL))
    colnames(count_data)[1] <- response_variable
    results = left_join(results, count_data, by = "Metabolite")
    return(results)

  } else if(interaction==TRUE & !missing(var2)){
    var1 = metadata[, var1]
    var2 = metadata[, var2]
    results <- llply(
      Vect, function(x) {
        models <- lm(data_mod[[x]] ~ var1 + var2 + var1 * var2, data = Mod)
      })
    names(results) <- Vect
    results <- lapply(results, anova)

    results <- sapply(results, cbind)
    results <- t(results)
    results <- as.data.frame(results[,"Pr(>F)"])
    results <- as.data.frame(t(results))
    results <- results[,1:3]
    colnames(results)[1] <- variable1
    colnames(results)[1] <- paste(colnames(results)[1], "pval", sep = ".")
    colnames(results)[2] <- variable2
    colnames(results)[2] <- paste(colnames(results)[2], "pval", sep = ".")
    colnames(results)[3] <- "Interaction.pval"

    #Compute adjusted p value
    results$padj = p.adjust(results[,1], method = p_adjust)
    colnames(results)[4] <- paste(colnames(results)[4], variable1, sep = ".")
    results$padj = p.adjust(results[,2], method = p_adjust)
    colnames(results)[5] <- paste(colnames(results)[5], variable2, sep = ".")
    results$padj = p.adjust(results[,3], method = p_adjust)
    colnames(results)[6] <- paste(colnames(results)[6], "Interaction", sep = ".")


    #Merge metadata, means & fold change, and t test results by metabolite name
    results <- cbind(rownames(results), data.frame(results, row.names=NULL))
    colnames(results)[1] <- response_variable

    count_data <- cbind(rownames(count_data), data.frame(count_data, row.names=NULL))
    colnames(count_data)[1] <- response_variable

    results$padj = p.adjust(results$`Interaction.pval`, method = "BH")
    results = left_join(results, count_data, by = "Metabolite")
    results <- results[ , !names(results) %in% c("padj")]
    return(results)
  }
}
