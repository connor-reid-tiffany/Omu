#'Check data for zeros across samples within factor levels.
#'Will determine if there are more zeros than a user specified threshold within any given factor level(s). Returns
#'a vector of Metabolites that are 0 above the threshold in any given factor level.
#'@param count_data A metabolomics count data frame
#'@param metadata Metadata dataframe for the metabolomics count data frame
#'@param response_variable String of the column header for the response variables,
#'usually "Metabolite"
#'@param numerator String of the first independent variable you wish to test. Defualt is NULL
#'@param denominator String of the second independent variable you wish to test. Default is NULL.
#'@param threshold Integer. A percentage threshold for the number of zeros in a Metabolite. Default is 25.
#'@param Factor A factor with levels to test for zeros.
#'@importFrom dplyr filter
#'@examples
#'
#'check_zeros(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#'Factor = "Treatment")
#'\dontshow{c57_nos2KO_mouse_countDF[1,3:31] <- 0}
#'check_zeros(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#' Factor = "Treatment",numerator = "Strep", denominator = "Mock", threshold = 10)
#'@export

check_zeros <- function(count_data, metadata, numerator = NULL, denominator = NULL, threshold = 25,
                        response_variable = "Metabolite", Factor){

rownames(count_data) <- count_data[,response_variable]
count_data[,response_variable] <- NULL

data_Int <- count_data[sapply(count_data, function(x) is.numeric(x))]
data_t <- as.data.frame(t(data_Int))
data_t <- as.data.frame(cbind(Sample = rownames(data_t), data_t))
data_t <- data_t[order(data_t$Sample),]
metadata <- metadata[order(metadata$Sample),]
Factor = metadata[, Factor]
data_t$Factor = Factor[match(metadata$Sample, data_t$Sample)]

if(!is.null(numerator) && !is.null(denominator)){

data_t <- filter(data_t, Factor==numerator|Factor==denominator)
rownames(data_t) <- data_t[,"Sample"]

}

data_t[,"Sample"] <- NULL

data_t_list <- split(data_t,f = as.factor(data_t$Factor))

percent_zero_metabolite <- lapply(data_t_list, function(x) data.frame(percent_zero = colSums(x==0)/nrow(x) *100))

percent_zero_metabolite <- lapply(percent_zero_metabolite, function(x){ x$Metabolite <- rownames(x);return(x)})

percent_zero_metabolite <- lapply(percent_zero_metabolite, function(x) data.frame(Metabolite = x[x[1] > threshold,"Metabolite"]))

metabolites <- do.call(rbind, percent_zero_metabolite)

metabolites$level <- gsub(x = rownames(metabolites), pattern = "\\..*", replacement = "")

if (length(rownames(metabolites))==0) {

  print("Number of zero values is below threshold across all metabolites.")

} else if (length(rownames(metabolites)) > 0){

  return(metabolites)

}

}
