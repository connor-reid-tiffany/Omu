#'Perform anova
#'@description Performs an anova across all response variables, followed by a Tukeys test on every possible
#'contrast in your model and calculates group means and fold changes for each contrast. Returns a list of
#'data frames for each contrast, and includes a dataframe of model residuals
#'@param count_data A metabolomics count data frame
#'@param metadata Metadata dataframe for the metabolomics count data frame
#'@param response_variable String of the column header for the response variables,
#'usually "Metabolite"
#'@param model A formual class object, see ?formula for more info on formulas in R.
#'an interaction between independent variables. Optional parameter
#'@param log_transform Boolean of TRUE or FALSE for whether or not you wish to log transform
#'your metabolite counts
#'@param method A string of 'anova', 'kruskal', or 'welch'. anova performs an anova with a post hoc
#' tukeys test, kruskal performs a kruskal wallis with a post hoc dunn test, welch performs a
#' welch's anova with a post hoc games howell test
#'@importFrom plyr llply
#'@importFrom broom tidy
#'@importFrom stats p.adjust
#'@importFrom stats aov
#'@importFrom stats TukeyHSD
#'@importFrom stats complete.cases
#'@importFrom stats aggregate
#'@importFrom stats as.formula
#'@importFrom stats terms
#'@importFrom FSA dunnTest
#'@importFrom rstatix games_howell_test
#'@examples
#'\dontshow{c57_nos2KO_mouse_countDF <- c57_nos2KO_mouse_countDF[1:12,];
#'c57_nos2KO_mouse_metadata <- c57_nos2KO_mouse_metadata;}
#'anova_df <- omu_anova(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#'response_variable = "Metabolite", model = ~ Treatment, log_transform = TRUE)
#'
#'anova_df <- omu_anova(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#'response_variable = "Metabolite", model = ~ Treatment + Background, log_transform = TRUE)
#'
#'anova_df <- omu_anova(count_data = c57_nos2KO_mouse_countDF, metadata = c57_nos2KO_mouse_metadata,
#'response_variable = "Metabolite", model = ~ Treatment + Background + Treatment*Background,
#'log_transform = TRUE)
#'
#'@export


omu_anova <- function (count_data, metadata, response_variable = "Metabolite", model, log_transform = FALSE,method="anova")
{

  if(any(names(count_data) %in% response_variable)==FALSE){

    stop("metabolomics data are missing the response variable column. Did you make a typo?")

  }

 model_characters <- strsplit(gsub("[^[:alnum:] ]", "", as.character(model)[-1]), " +")[[1]]

 if(all(model_characters %in% names(metadata))==FALSE){

   stop("One or more model terms do not match column names in metadata. Did you make a typo?")

 }

 if(isTRUE(method %in% c("anova", "kruskal", "welch"))==FALSE){

   stop("method must be anova, kruskal, or welch. Did you make a typo?")

 }
  #extract variables from model object
  model_terms <- terms(model)

  variable_vector <- attr(attr(model_terms,which = "factors"),
                          which = "dimnames")[[1]]

  #conditionals to check for errors
  metadata <- as.data.frame(sapply(metadata, function(x) x <- as.character(x)))

  if(identical(sort(as.character(colnames(count_data)[unlist(lapply(count_data, is.numeric))])), sort(as.character(metadata$Sample))==FALSE)){

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

  #function to check for zeros in model terms
  check_zeros_anova <- function(count_data = count_data, metadata = metadata, variable_vector = variable_vector){

    check_zeros_list <- list()

    for (v in 1:length(variable_vector)){

      check_zeros_list[[v]] <- check_zeros(count_data = count_data, metadata = metadata, Factor = v)

    }

    return(check_zeros_list)

  }

  if (any(length(check_zeros_anova(count_data = count_data, metadata = metadata, variable_vector = variable_vector))) < 1) {

    warning("There are zero values in at least 25 percent of your samples within at least one of your Factor
            levels for these metabolites. Consider using check_zeros to subset your data.")

    print(unique(check_zeros_anova()[,1]))

  }
  #set aside original count_data
  count_data_o <- count_data[,!names(count_data) %in% c("KEGG", "Class", "Subclass_1", "Subclass_2", "Subclass_3", "Subclass_4")]
  #add backtick quotes to metabolites so that non-syntatic names can be converted into model formulae later
  count_data$Metabolite <- paste0("`", count_data$Metabolite)
  count_data$Metabolite <- paste0(count_data$Metabolite, "`")

  #format data for anova, make separate non-integer data to merge later
  rownames(count_data) <- count_data[, response_variable]
  count_data_character <- count_data[!sapply(count_data, is.numeric)]
  count_data[, response_variable] <- NULL
  data_Int <- count_data[sapply(count_data, is.numeric)]
  data_Transpose <- as.data.frame(t(data_Int))
  #make vector of metabolite names
  Vect = colnames(data_Transpose)

  if (log_transform==TRUE){

    data_mod <- log(data_Transpose)

  } else if (log_transform==FALSE){

    data_mod <- data_Transpose

  }

  #establish model terms as factors

  model_factors <- list()

  for (t in variable_vector){

    model_factors[[t]] <- metadata[,t]

  }

  if(method=="kruskal"& length(model_factors) > 1){

    stop("Method kruskal can only take a model with one term.")

  }

  if(method=="welch"& length(model_factors) > 1){

    stop("Method welch can only take a model with one term.")

  }

  if(method=="kruskal"&length(levels(as.factor(model_factors[[1]]))) < 3){

    stop("Method kruskal needs at least 3 levels in the model term.")

  }

  Vect_list <- as.list(Vect)

  #create list of models
  model_list <- lapply(Vect_list, function(x) {

    new_formula <- paste(paste(x, "~"), model)

    return(new_formula[2])

  })

  model_list <- lapply(model_list, as.formula)

  #remove backtick quotes from metabolites
  names(data_mod) <- gsub("`", "", names(data_mod), fixed = TRUE)
  Vect <- gsub(pattern = "`", replacement = "", x = Vect)
  count_data_character$Metabolite <- gsub(pattern = "`", replacement = "", x = count_data_character$Metabolite)

  #add metadata variables
  data_mod <- cbind(data_mod, metadata)
  #run anova
  if(method=="anova"){

    results_aov <- llply(model_list, function(x) { aov(x, data_mod)})

  }else if(method=="kruskal"){

    results_aov <- llply(model_list, function(x) { dunnTest(x,  data_mod)})

  }else if(method=="welch"){

    results_aov <- llply(model_list, function(x){games_howell_test(formula = x,data = data_mod)})

    }
  #combine information from anova model with info from tukeys post hoc test
  names(results_aov) <- Vect

  if(method=="anova"){

  results_tukey <- lapply(results_aov, TukeyHSD)
  results_aov <- lapply(results_aov, tidy)
  results_tukey <- lapply(results_tukey, tidy)
  results_tukey <- lapply(results_tukey, function(x){colnames(x)[length(x)] <- "padj"; return(x)})

  }else if(method=="kruskal"){

  #reduce to results dfs
  results <- lapply(results_aov, function(x) x$res)
  names(results) <- Vect
  results <- lapply(results, function(x){colnames(x)[1]<- "contrast"; return(x)})
  results <- lapply(results, function(x){colnames(x)[4]<- "padj";return(x)})
  results <- lapply(results, function(x){x$contrast <- gsub(pattern = " ", replacement = "", x = x$contrast); return(x)})
  results <- lapply(results, function(x){x$term <- names(model_factors)[1];return(x)})

  }else if(method=="welch"){

    results <- results_aov
    results <- lapply(results, function(x){colnames(x)[7] <- "padj"; return(x)})
    results <- lapply(results, function(x){x$contrast <- paste0(x$group1,"-",x$group2); return(x)})
    results <- lapply(results, function(x){x$term <- names(model_factors)[1];return(x)})

  }


  add_residuals_placeholder <- function(x){

    resid_empty_df <- data.frame(term = "Residuals", contrast = NA, null.value = NA,
                                 estimate = NA, conf.low = NA, conf.high = NA, padj = NA)
    x <- rbind(x, resid_empty_df)

    return(x)

  }

  if(method=="anova"){
  #change this part to add a residuals column by matching via metabolite
  results_tukey <- lapply(results_tukey, add_residuals_placeholder)

  #anova specific code
  results <- Map(merge, results_tukey, results_aov, by = "term")
  results <- lapply(results, as.data.frame)
  }
  #add metabolite column and then rbind all response variable dataframes into one tidy dataframe
  names_list <- as.list(names(results))
  names_list <- lapply(names_list, function(x) data.frame("Metabolite" = x))

  combined_list <- Map(c, results, names_list)
  combined_list <- lapply(combined_list, as.data.frame)

  combined_data_frame <- do.call("rbind", combined_list)
  combined_data_frame_merge <- merge(combined_data_frame, count_data_character, response_variable)
  #colnames(combined_data_frame_merge)[8] <- "padj"

  #calculate fold change for terms
  #for two way anova

  #create an index for each unique combination of grouped factor comparisons
  #observed in tukeys test, then use the index to split the transposed data frame
  #into a dataframe for each unique group comparison
  #then calculate means and fold changes for each unique group comparison
  #make dataframe with variable and contrast columns
  FC_groups <- combined_data_frame_merge[, c("term", "contrast")]
  #remove residuals
  FC_groups <- FC_groups[complete.cases(FC_groups),]
  #take unique contrasts
  FC_groups <- FC_groups[!duplicated(FC_groups$contrast),]
  FC_groups$vars <- strsplit(FC_groups$term, ":")
  #split into a list by contrast
  FC_groups_l <- split(FC_groups, FC_groups$contrast)

  #make appropriate metadata dataframes for each contrast
  FC_groups_metadata <- lapply(FC_groups_l, function(x){

    variables <- unlist(x$vars)

    x2 <- data.frame(matrix(nrow = nrow(metadata)))

    for (v in variables){

      x2[,v] <- metadata[, names(metadata) %in% v,drop = F]

    }

    rownames(x2) <- metadata$Sample
    x2 <- x2[,-1,drop = FALSE]

    return(x2)

  })

  #make grouped column by pasting all other columns in dataframe together
  FC_groups_metadata<- lapply(FC_groups_metadata, function(x)

  {x$Grouped <- do.call(paste, c(x[colnames(x)], sep=":")); return(x)})

  #make a list of corresponding constrasts
  contrasts_list <- data.frame(contrasts =names(FC_groups_metadata))
  contrasts_list <- split(contrasts_list, f = contrasts_list$contrasts)

  contrasts_list <- lapply(contrasts_list, function(x){

    x[rep(seq_len(nrow(x)), each = nrow(metadata)),,drop = FALSE ]


  })
  #add constrasts column to metadata
  contrasts_metadata <- Map("cbind", FC_groups_metadata, contrasts_list)
  contrasts_metadata <- lapply(contrasts_metadata, function(x)
  {x$contrasts <- strsplit(x$contrasts, "-"); return(x)})
  #subset metadata dataframes by contrast column
  contrasts_metadata <- lapply(contrasts_metadata, function(x)

    x <- x[which(x$Grouped %in% unlist(x$contrasts)),]
  )
  #merge metadata dataframes with dataTranspose via sample values
  contrasts_metadata <- lapply(contrasts_metadata, function(x){x$Sample <- rownames(x); return(x)})
  #remove ticks from dataTranspose column names and convert rownames into a column
  colnames(data_Transpose) <- gsub(pattern = "`", replacement = "", x = colnames(data_Transpose))
  data_Transpose$Sample <- rownames(data_Transpose)

  contrasts_metadata <- lapply(contrasts_metadata, function(x){

    x <- merge(x, data_Transpose, "Sample")

  })
  #reorder grouped factor levels based on contrasts
  contrasts_metadata <- lapply(contrasts_metadata, function(x){

    x$Grouped <- factor(x$Grouped, levels = unlist(x$contrasts[1]))
    x$Grouped <- as.factor(x$Grouped)

    return(x)
  })

  #calculate group means for each dataframe
  calc_means <- function(x) {

    cols_to_keep_num <- sapply(x, is.numeric)

    cols_to_keep_char <- names(x) %in% c("Grouped", "contrasts")

    x <- x[,cols_to_keep_num | cols_to_keep_char]

    x2 <- aggregate(x = x[,!names(x) %in% c("Grouped", "contrasts")],
                    list(x$Grouped), mean)

    return(x2)
  }

  metabo_list_means <- lapply(contrasts_metadata, calc_means)
  #calculate fold changes
  calc_FC <- function(x) {

    rownames(x) <- x$Group.1

    x <- x[,-1]

    x[paste0(rownames(x[1,]), "-",rownames(x[2,])),] <- x[1,]/x[2,]

    x["log2FoldChange",] <- log2(x[3,])

    return(x)

  }

  metabo_list_means <- lapply(metabo_list_means, calc_FC)
  metabo_list_means_t <- lapply(metabo_list_means, function(x) as.data.frame(t(x)))
  metabo_list_means_t <- lapply(metabo_list_means_t, function(x){x$Metabolite <- rownames(x) ; return(x)})
  #split combined_data_frame_merge into a list by contrasts, then merge those data frames
  #with metabo list means by metabolite
  if(method=="anova"){

  residuals_df <- combined_data_frame_merge[combined_data_frame_merge$term=="Residuals", c("Metabolite","df","sumsq","meansq","KEGG")]
  combined_data_frame_merge <- combined_data_frame_merge[!combined_data_frame_merge$term=="Residuals",]

  }

  combined_data_list <- split(combined_data_frame_merge, f = as.factor(combined_data_frame_merge$contrast))
  data_final <- Map(merge, combined_data_list, metabo_list_means_t, by='Metabolite', all=TRUE)

  #add residuals dataframe to the list, its done!
  if(method=="anova"){

  data_final$Residuals <- residuals_df

  }

  data_final <- lapply(data_final, function(x){x <- merge(x, count_data_o, by = "Metabolite"); return(x)})
  data_final <- lapply(data_final, function(x){class(x) <- append(class(x), "cpd"); return(x)})

  return(data_final)

}
