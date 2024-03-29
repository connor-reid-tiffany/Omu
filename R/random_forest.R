#' random_forest
#' Perform a classification or regression random forest model
#' @description a wrapper built around the randomForest function from package randomForest.
#' Returns a list with a randomForest object list, training data set, testing data set,
#' metabolite metadata, and confusion matrices for training and testing data
#' (if type was classification).
#' @param count_data Metabolomics data
#' @param metadata sample data
#' @param model a model of format variable ~.
#' @param training_proportion a numeric vector of length 2, first element is the percent of
#' samples to use for training the model, second element is the percent of samples used to
#' test the models accuracy
#' @param n_tree number of decision trees to create
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @importFrom caret confusionMatrix
#' @examples
#' rf_list <- random_forest(count_data = c57_nos2KO_mouse_countDF,metadata = c57_nos2KO_mouse_metadata,
#' model = Treatment ~.,training_proportion = c(60,40),n_tree = 500)
#' @export

random_forest <- function(count_data, metadata, model, training_proportion = c(80,20), n_tree = 500){

if(identical(sort(as.character(colnames(count_data)[unlist(lapply(count_data, is.numeric))])), sort(as.character(metadata$Sample))==FALSE)){

    stop("Sample names in count_data and metadata do not match.")

  }

  if(any(colnames(metadata)=="Sample")==FALSE){

    stop("metadata is missing Sample column")

  }

  model_characters <- strsplit(gsub("[^[:alnum:] ]", "", as.character(model)[-1]), " +")[[1]]

  if(all(model_characters %in% names(metadata))==FALSE){

    stop("One or more model terms do not match column names in metadata. Did you make a typo?")

  }
  #set RNG so model output is reproducible
  set.seed(123)
  #parse data to handle non syntatic metabolite names so they don't throw an error in the model
  count_data$Metabolite <- gsub(pattern = ",", replacement = "_", count_data$Metabolite)
  count_data$Metabolite <- gsub(pattern = " ", replacement = "_", count_data$Metabolite)
  count_data$Metabolite <- gsub(pattern = "'", replacement = "", count_data$Metabolite)
  count_data$Metabolite <- gsub(pattern = "\\(", replacement = "_", count_data$Metabolite)
  count_data$Metabolite <- gsub(pattern = "\\)", replacement = "_", count_data$Metabolite)
  count_data$Metabolite <- sub("^", "ID_", count_data$Metabolite)
  #generate metabolite metadata for use in downstream plotting functions
  metabo_meta <- assign_hierarchy(count_data = count_data, keep_unknowns = TRUE, identifier = "KEGG")
  metabo_meta <- metabo_meta[,!sapply(metabo_meta, is.numeric)]
  #generate the metabolite matrix and derive training and testing sets for the model
  rownames(count_data) <- count_data$Metabolite
  count_data_num <- count_data[,sapply(count_data, is.numeric)]
  count_data_num <- as.data.frame(t(count_data_num))
  count_data_num$Sample <- rownames(count_data_num)
  metadata_sub <- metadata[,c(as.character(model)[2], "Sample")]
  count_data_num <- merge(count_data_num, metadata_sub, by = "Sample")
  count_data_num <- count_data_num[,-1]
  index <- sample(2, nrow(count_data_num), replace = TRUE, prob = training_proportion)
  train <- count_data_num[index==1,]
  test <- count_data_num[index==2,]
  #fit the random forest model
  rfFit <- randomForest(model, data=train, proximity=TRUE, ntree=n_tree,importance = TRUE)
  #for classification, confusion matrices on the binary/multinomial variable can be produced
  if(rfFit$type=="classification"){

  pred_test <- predict(rfFit, newdata = test)

  test_conf_mat <- confusionMatrix(pred_test, test[,length(test)])

  pred_train <- predict(rfFit, train)

  train_conf_mat <- confusionMatrix(pred_train, train[,length(train)])

  rf_list <- list(rf = rfFit, train = train, test = test, cm_train = train_conf_mat, cm_test = test_conf_mat,
                  metabolite_meta = metabo_meta)
  #for regression confusion matrices can not be produced
  }else if(rfFit$type=="regression"){

    rf_list <- list(rf = rfFit, train = train, test = test,
                    metabolite_meta = metabo_meta)

  }
  return(rf_list)

}
