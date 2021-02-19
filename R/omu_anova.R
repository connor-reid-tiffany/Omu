#'Perform anova
#'Performs an anova across all response variables. The function can take a maximum of 2
#'independent variables and perform an interaction term between them.
#'@param count_data A metabolomics count data frame
#'@param metadata Metadata dataframe for the metabolomics count data frame
#'@param response_variable String of the column header for the response variables,
#'usually "Metabolite"
#'@param var1 String of the first independent variable you wish to test
#'@param var2 String of the second independent variable you wish to test. Default is NULL.
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
#'\dontshow{c57_nos2KO_mouse_countDF <- c57_nos2KO_mouse_countDF[1:12,]}
#'anova_df <- omu_anova(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#'response_variable = "Metabolite", var1 = "Treatment", var2 = "Background", log_transfor = TRUE,
#'p_adjust = "BH", interaction = TRUE)
#'@export


omu_anova <- function (count_data, metadata, response_variable, var1, var2 = NULL,
          interaction, log_transform, p_adjust)
{

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

if (is.null(var2)){

variable1 = var1

} else if (!is.null(var2)){

variable1 = var1
variable2 = var2}

rownames(count_data) <- count_data[, response_variable]
count_data[, response_variable] <- NULL
data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]
data_Transpose <- as.data.frame(t(data_Int))
data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose),
                                      data_Transpose))
data_Transpose <- merge(data_Transpose, metadata, "Sample")
rownames(data_Transpose) <- data_Transpose[,"Sample"]
nums <- sapply(data_Transpose, is.numeric)
factors <- sapply(data_Transpose, is.character)
data_Num <- data.frame(lapply(data_Transpose[, nums], function(x) as.numeric(as.integer(x))),
                       check.names = F, row.names = rownames(data_Transpose))
Vect = colnames(data_Num)
data_Ln <- log(data_Num)
data_Ln <- as.data.frame(cbind(Sample = rownames(data_Ln),
                               data_Ln))
data_Fact <- data.frame(Sample = data_Transpose[, factors])
data_Ln = merge(data_Ln, data_Fact, by = "Sample")
data_Ln <- data_Ln[, !names(data_Ln) %in% c("Sample", "Sample.1")]
data_Num <- as.data.frame(cbind(Sample = rownames(data_Num),
                                data_Num))
data_Num = merge(data_Num, data_Fact, by = "Sample")
data_Num <- data_Num[, !names(data_Num) %in% c("Sample",
                                               "Sample.1")]
if (log_transform == FALSE) {
  data_mod = data_Num
}
else if (log_transform == TRUE) {
  data_mod = data_Ln
}
Mod = data_Fact
Mod = Mod[, !names(Mod) %in% c("Sample", "Sample.1")]
if (is.null(var2) & interaction == FALSE) {
  var1 = metadata[, var1]
  results <- llply(Vect, function(x) {
    models <- lm(data_mod[[x]] ~ var1, data = Mod)
  })
  names(results) <- Vect
  results <- lapply(results, anova)
  results <- sapply(results, cbind)
  results <- t(results)
  results <- as.data.frame(results[, "Pr(>F)"])
  results <- as.data.frame(t(results))
  colnames(results)[1] <- variable1
  colnames(results)[1] <- paste(colnames(results)[1],
                                "pval", sep = ".")
  results$padj = p.adjust(results[, 1], method = p_adjust)
  colnames(results)[3] <- paste(colnames(results)[3],
                                variable1, sep = ".")
  count_data <- cbind(rownames(count_data), data.frame(count_data,
                                                       row.names = NULL))
  colnames(count_data)[1] <- response_variable
  results <- cbind(rownames(results), data.frame(results,
                                                 row.names = NULL))
  colnames(results)[1] <- response_variable
  results = left_join(results, count_data, by = response_variable)
  results = results[, !(colnames(results) %in% c("V2"))]
  return(results)
}
else if (interaction == FALSE) {
  var1 = metadata[, var1]
  var2 = metadata[, var2]
  results <- llply(Vect, function(x) {
    models <- lm(data_mod[[x]] ~ var1 + var2, data = Mod)
  })
  names(results) <- Vect
  results <- lapply(results, anova)
  results <- sapply(results, cbind)
  results <- t(results)
  results <- as.data.frame(results[, "Pr(>F)"])
  results <- as.data.frame(t(results))
  results <- results[, 1:2]
  colnames(results)[1] <- variable1
  colnames(results)[1] <- paste(colnames(results)[1],
                                "pval", sep = ".")
  colnames(results)[2] <- variable2
  colnames(results)[2] <- paste(colnames(results)[2],
                                "pval", sep = ".")
  results$padj = p.adjust(results[, 1], method = p_adjust)
  colnames(results)[3] <- paste(colnames(results)[3],
                                variable1, sep = ".")
  results$padj = p.adjust(results[, 2], method = p_adjust)
  colnames(results)[4] <- paste(colnames(results)[4],
                                variable2, sep = ".")
  results <- cbind(rownames(results), data.frame(results,
                                                 row.names = NULL))
  colnames(results)[1] <- response_variable
  count_data <- cbind(rownames(count_data), data.frame(count_data,
                                                       row.names = NULL))
  colnames(count_data)[1] <- response_variable
  results = left_join(results, count_data, by = "Metabolite")
  return(results)
}
else if (interaction == TRUE & !is.null(var2)) {
  var1 = metadata[, var1]
  var2 = metadata[, var2]
  results <- llply(Vect, function(x) {
    models <- lm(data_mod[[x]] ~ var1 + var2 + var1 *
                   var2, data = Mod)
  })
  names(results) <- Vect
  results <- lapply(results, anova)
  results <- sapply(results, cbind)
  results <- t(results)
  results <- as.data.frame(results[, "Pr(>F)"])
  colnames(results) <- Vect
  results <- as.data.frame(t(results))
  results <- results[, 1:3]
  colnames(results)[1] <- variable1
  colnames(results)[1] <- paste(colnames(results)[1],
                                "pval", sep = ".")
  colnames(results)[2] <- variable2
  colnames(results)[2] <- paste(colnames(results)[2],
                                "pval", sep = ".")
  colnames(results)[3] <- "Interaction.pval"
  results$padj = p.adjust(results[, 1], method = p_adjust)
  colnames(results)[4] <- paste(colnames(results)[4],
                                variable1, sep = ".")
  results$padj = p.adjust(results[, 2], method = p_adjust)
  colnames(results)[5] <- paste(colnames(results)[5],
                                variable2, sep = ".")
  results$padj = p.adjust(results[, 3], method = p_adjust)
  colnames(results)[6] <- paste(colnames(results)[6],
                                "Interaction", sep = ".")
  results <- cbind(rownames(results), data.frame(results,
                                                 row.names = NULL))
  colnames(results)[1] <- response_variable
  count_data <- cbind(rownames(count_data), data.frame(count_data,
                                                           row.names = NULL))
      colnames(count_data)[1] <- response_variable
  results <- left_join(results, count_data, by = "Metabolite")
  return(results)
}
}
