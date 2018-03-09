#'anova_function
#'Function to apply an anova across all response variables. Options for running the model on 1 or 2 independent
#'variables individually, or on the interaction between two independent variables
#'@export
#'@param data A metabolomics count data frame
#'@param colData Meta data dataframe for the metabolomics count data frame
#'@param response_variable String of the column header for the response variables, usually "Metabolite"
#'@param Var1 String of the first independent variable you wish to test
#'@param Var2 String of the second independent variable you wish to test. Optional parameter
#'@param interaction Boolean of TRUE or FALSE for whether or not you wish to model an interaction between 
#'independent variables. Optional parameter
#'@param log_transform Boolean of TRUE or FALSE for whether or not you wish to log transform your counts
#'@examples anova_function(data = metabolomics_counts, colData = metadata, response_variable = "Metabolite", Var1 = "Treatment", log_transform = TRUE)
#'anova_function()

anova_function <- function(data, colData, response_variable, Var1, Var2, interaction, log_transform){

data = column_to_rownames(df = data, var = response_variable)
data_Int <- data[sapply(data, function(x) is.numeric(x))]

#Transform for 'normalization' 
data_Transpose <- as.data.frame(t(data_Int))
data_Transpose <- as.data.frame(cbind(Sample = rownames(data_Transpose), data_Transpose))
data_Transpose = as.data.frame(cbind(data_Transpose, colData))
nums <- sapply(data_Transpose, is.integer)
factors <- sapply(data_Transpose, is.factor)

#Convert class to numeric, natural log transform, and remerge num columns with factor columns by sample name
data_Num <- data.frame(lapply(data_Transpose[ , nums], function(x) as.numeric(as.integer(x))),check.names=F, row.names = rownames(data_Transpose))
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

if(missing(Var2) & missing(interaction)){
Var1 = colData[,Var1]
results <- llply(
  Vect, function(x) {
    models <- lm(data_mod[[x]] ~ Var1, data = Mod)
    })
names(results) <- Vect
results <- lapply(results, anova)

results <- sapply(results, cbind)
results <- t(results)
results <- as.data.frame(results[,"Pr(>F)"])
results <- as.data.frame(t(results))
colnames(results)[1] <- "Var1 pval"

#Merge metadata, means & fold change, and t test results by metabolite name
results = rownames_to_column(results, "Metabolite")
results$padj = p.adjust(results$`Var1 pval`, method = "BH")
data = rownames_to_column(data, "Metabolite")
results = left_join(results, data, by = response_variable)
results = results[, !(colnames(results) %in% c("V2"))]
return(results)

} else if(missing(interaction)){
  Var1 = colData[, Var1]
  Var2 = colData[, Var2]
  results <- llply(
    Vect, function(x) {
      models <- lm(data_mod[[x]] ~ Var1 + Var2, data = Mod)
    })
  names(results) <- Vect
  results <- lapply(results, anova)
  
  results <- sapply(results, cbind)
  results <- t(results)
  results <- as.data.frame(results[,"Pr(>F)"])
  results <- as.data.frame(t(results))
  results <- results[,1:2]
  colnames(results)[1] <- "Var1 pval"
  colnames(results)[2] <- "Var2 pval"
  
  #Merge metadata, means & fold change, and t test results by metabolite name
  results = rownames_to_column(results, "Metabolite")
  data = rownames_to_column(data, "Metabolite")
  results$padj = p.adjust(results$`Var2 pval`, method = "BH")
  results = left_join(results, data, by = "Metabolite")
  return(results)
  
} else if(!missing(interaction) & !missing(Var2)){
  Var1 = colData[, Var1]
  Var2 = colData[, Var2]
  results <- llply(
    Vect, function(x) {
      models <- lm(data_mod[[x]] ~ Var1 + Var2 + Var1 * Var2, data = Mod)
    })
  names(results) <- Vect
  results <- lapply(results, anova)
  
  results <- sapply(results, cbind)
  results <- t(results)
  results <- as.data.frame(results[,"Pr(>F)"])
  results <- as.data.frame(t(results))
  results <- results[,1:3]
  colnames(results)[1] <- "Var1 pval"
  colnames(results)[2] <- "Var2 pval"
  colnames(results)[3] <- "Interaction pval"
  
  #Merge metadata, means & fold change, and t test results by metabolite name
  results = rownames_to_column(results, "Metabolite")
  data = rownames_to_column(data, "Metabolite")
  results$padj = p.adjust(results$`Interaction pval`, method = "BH")
  results = left_join(results, data, by = "Metabolite")
  return(results)
}
}



