#' plot_variable_importance
#' @description Plot the variable importance from a random forest model. Mean Decrease Gini for Classification and
#' %IncMSE for regression.
#' @param rf_list The output from the random_forest function
#' @param color Metabolite metadata to color by
#' @param n_metabolites The number of metabolites to include. Metabolites are sorted by decreasing importance.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 labs
#' @examples
#' rf_list <- random_forest(c57_nos2KO_mouse_countDF,c57_nos2KO_mouse_metadata,
#' Treatment ~.,c(60,40),500)
#' plot_variable_importance(rf_list = rf_list, color = "Class", n_metabolites = 10)
#' @export

plot_variable_importance <- function(rf_list, color="Class", n_metabolites=10){

  if(is.null(rf_list$rf)==TRUE){

    stop("rf_list is mising randomForest output")

  }

  if(any(names(rf_list$metabolite_meta) %in% color)==FALSE){

    stop("color must be a metabolite metadata level, i.e. Class, Subclass_1, etc. did
    you make a typo?")

  }

  #address silly CRAN note
  Metabolite <- NULL
  MeanDecreaseGini <- NULL
  MeanDecreaseAccuracy <- NULL
  `%IncMSE` <- NULL

  #sort the importance data frame based off the appropriate variable
  importance <- as.data.frame(rf_list$rf$importance)
  if(rf_list$rf$type=="classification"){

  importance <- importance[order(rf_list$rf$importance[,4],decreasing=TRUE),]

  }else if(rf_list$rf$type=="regression"){

    importance <- importance[order(rf_list$rf$importance[,1],decreasing=TRUE),]

  }
  #subset data to selected number of top metabolites
  importance$Metabolite <- rownames(importance)
  importance <- importance[1:n_metabolites,]
  #merge with metabolite metadata
  importance <- merge(importance, rf_list$metabolite_meta, by = "Metabolite")
  #proudce plots
  if(rf_list$rf$type=="classification"){

  ggplot(data = importance, aes(x = reorder(Metabolite, MeanDecreaseGini), y = MeanDecreaseGini, color = importance[,color])) +
    geom_segment( aes(xend=Metabolite, y=0, yend=MeanDecreaseGini), size = 1) +
    coord_flip() +
    geom_point(aes(size = MeanDecreaseAccuracy), alpha=0.6) +
    theme(axis.title.y = element_blank(), axis.text = element_text(size = 14), axis.title.x = element_text(size = 16), legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))  +
    labs(y="Mean Decrease Gini",
         col=color, size="Mean Decrease Accuracy")

  }else if(rf_list$rf$type=="regression"){

    ggplot(data = importance, aes(x = reorder(Metabolite, `%IncMSE`), y = `%IncMSE`, color = importance[,color])) +
      geom_segment( aes(xend=Metabolite, y=0, yend=`%IncMSE`), size = 1) +
      coord_flip() +
      geom_point(aes(size = MeanDecreaseAccuracy), alpha=0.6) +
      theme(axis.title.y = element_blank(), axis.text = element_text(size = 14), axis.title.x = element_text(size = 16), legend.text = element_text(size = 14),
            legend.title = element_text(size = 16))  +
      labs(y="%IncMSE",
           col=color, size="IncNodePurity")

  }
}

#' plot_rf_PCA
#' @description PCA plot of the proximity matrix from a random forest classification model
#' @param rf_list The output from the random_forest function. This only works on classification models.
#' @param color A grouping factor. Use the one that was the LHS of your model parameter in the random_forest funciton
#' @param size The number for point size in the plot
#' @importFrom ggplot2 autoplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom stats prcomp
#' @examples
#' rf_list <- random_forest(c57_nos2KO_mouse_countDF,c57_nos2KO_mouse_metadata,
#' Treatment ~.,c(60,40),500)
#' plot_rf_PCA(rf_list = rf_list, color = "Treatment", size = 1.5)
#' @export

plot_rf_PCA <- function(rf_list, color, size){

  if(is.null(rf_list$rf)==TRUE){

    stop("rf_list is mising randomForest output")

  }

  if(any(names(rf_list$train) %in% color)==FALSE){

    stop("color variable not found in data. did you make a typo?")

  }

  #PCA is computed on the proximity matrix and regression models do not produce a proximity matrix
  if(rf_list$rf$type=="regression"){

    stop("PCA can't be performed on a random forest model of type regression")
  }

  PCA <- prcomp(rf_list$rf$proximity)

  autoplot(PCA, data = rf_list$train, colour = color, size = size) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
                                                                           legend.title = element_blank(), legend.text = element_text(size = 14))




}
