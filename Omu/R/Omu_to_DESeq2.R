#'Omu_to_DESeq2
#'@description Takes an Omu count matrix and generates a DESeq2 statistical model
#'@param countMatrix A metabolomics integer class matrix that has undergone feature scaling, as well as hierarchy annotation through the Omu python module.
#'@param colData A data frame of character vectors assigning metadata to the discrete variables of the countMatrix. Therefore, ncol of colData must equal nrow of countMatrix.
#'@param design The statistical model. I.e. ~ Treatment, ~Treatment + Condition, Treatment + Condition:Treatment (for interaction term), etc. See the DESeq2 manual for an in depth explanation of modeling.
#'@param comparison A character list for the comparisons you wish to view, listed as c("Factor of Interest", "Numerator level", "Denominator Level")
#'@param test Either "LRT" (likelihood ratio test) or "Wald". For LRT, a reduced design is required along with the full design
#'@param fitType "parametric", "local", or "mean"
#'@param reduced_design A reduced model missing a factor of interest. For LRT
#'@param group TRUE or FALSE value. Use group = TRUE if you want to combine Factors of interest into a single factor consisting of all level combinations. If group = TRUE, design must = ~ group. See Vignette for a more detailed explanation on modeling.
 #'@export
Omu_to_DESeq2 <- function(countMatrix, colData, design, comparison, test, fitType, reduced_design, group){
  if (!missing(reduced_design) & group == FALSE){
    countMatrix3 <- cbind(Metabolite = rownames(countMatrix), countMatrix)
    METABO = countMatrix3[, c("KEGG", "Metabolite", "Class", "Subclass_1", "Subclass_2", "Subclass_3", "Subclass_4")]
    countMatrix2 <- subset.matrix(countMatrix, select = -c(KEGG, Class, Subclass_1, Subclass_2, Subclass_3, Subclass_4))
    class(countMatrix2) <- "integer"  
      
      DESeq_object <- DESeqDataSetFromMatrix(countData = countMatrix2, colData = colData, design = design)
        ddsLRT <- DESeq(DESeq_object, test = test, fitType = fitType, reduced = reduced_design)
          res <- results(ddsLRT, cooksCutoff = FALSE, contrast = comparison)
          res = cbind(as(res, "data.frame"))
          res = cbind(Metabolite = rownames(res), res)
          res = merge(res, METABO, by="Metabolite")
          
    }else if (!missing(reduced_design) & group == TRUE){
        countMatrix3 <- cbind(Metabolite = rownames(countMatrix), countMatrix)
        METABO = countMatrix3[, c("KEGG","Metabolite", "Class", "Subclass_1", "Subclass_2", "Subclass_3", "Subclass_4")]
        countMatrix2 <- subset.matrix(countMatrix, select = -c(KEGG, Class, Subclass_1, Subclass_2, Subclass_3, Subclass_4))
        class(countMatrix2) <- "integer"  
        colData2 <- unite_(colData, from = colnames(colData)[-1] , col = 'group', sep = "")
        DESeq_object <- DESeqDataSetFromMatrix(countData = countMatrix2, colData = colData2, design = design)
          ddsLRT <- DESeq(DESeq_object, test = test, fitType = fitType, reduced = reduced_design)
            res <- results(ddsLRT, cooksCutoff = FALSE, contrast = comparison)
            res = cbind(as(res, "data.frame"))
            res = cbind(Metabolite = rownames(res), res)
            res = merge(res, METABO, by="Metabolite")
      
      }else{
        if (group == FALSE){
          countMatrix3 <- cbind(Metabolite = rownames(countMatrix), countMatrix)
          METABO = countMatrix3[, c("KEGG","Metabolite", "Class", "Subclass_1", "Subclass_2", "Subclass_3", "Subclass_4")]
          countMatrix2 <- subset.matrix(countMatrix, select = -c(KEGG, Class, Subclass_1, Subclass_2, Subclass_3, Subclass_4))
          class(countMatrix2) <- "integer" 
    
           DESeq_object <- DESeqDataSetFromMatrix(countData = countMatrix2, colData = colData, design = design)
            dds <- DESeq(DESeq_object, test = test, fitType = fitType)
            
              res <- results(dds, cooksCutoff = FALSE, contrast = comparison)
              res = cbind(as(res, "data.frame"))
              res = cbind(Metabolite = rownames(res), res)
              res = merge(res, METABO, by="Metabolite")
        
        }else if (group == TRUE){
          countMatrix3 <- cbind(Metabolite = rownames(countMatrix), countMatrix)
          METABO = countMatrix3[, c("KEGG","Metabolite", "Class", "Subclass_1", "Subclass_2", "Subclass_3", "Subclass_4")]
          countMatrix2 <- subset.matrix(countMatrix, select = -c(KEGG,Class, Subclass_1, Subclass_2, Subclass_3, Subclass_4))
          class(countMatrix2) <- "integer"  
          colData2 <- unite_(colData, from = colnames(colData)[-1] , col = 'group', sep = "")
            DESeq_object <- DESeqDataSetFromMatrix(countData = countMatrix2, colData = colData2, design = design)
              dds <- DESeq(DESeq_object, test = test, fitType = fitType)
                res <- results(dds, cooksCutoff = FALSE, contrast = comparison)
                res = cbind(as(res, "data.frame"))
                res = cbind(Metabolite = rownames(res), res)
                res = merge(res, METABO, by="Metabolite")}
  
  }
}
