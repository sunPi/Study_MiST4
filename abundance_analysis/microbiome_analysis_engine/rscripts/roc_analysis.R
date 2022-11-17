# ------------------ Functions -------------------------------------------------
requirements  <- function(pkgs){ # This function checks package requirements and install them if they are missing
  suppressMessages(if (!require("BiocManager", character.only = TRUE)) { # First check via R BioConductior
    install.packages("BiocManager")
    BiocManager::install()
  } else {
    ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
    if (any(!ipkgs)) {
      BiocManager::install(pkgs[!ipkgs])
      install.packages(pkgs[!ipkgs])
    } else {
      message("\n\nCool! your machine has everything is needed.\n\n")
    }
  })
  
  print("Loading required packages...")
  library(pacman)
  pacman::p_load(pkgs, install = TRUE, character.only = TRUE) # Check via RCran and other repositories
  return(pacman::p_loaded()) # Return loaded packages
} # Helper function for packages & requirements
checkIfBinary <- function(response){
  if(is.numeric(response)){
    for(n in response){
      if(n != 0 && n != 1){
        
        stop("Response variable is not in a binary form. Consider reformating the column into '0' for 'no response' 
           and '1' for 'response'!")
      } else{
        return(TRUE)
        # print(paste0(n, " is ok!"))
      }
    }
  } else{
    stop("Response variable is not numeric. Convert categorical to numeric and re-run the script!")
  }
} 
getYd         <- function(cm){
  tn <- cm$table[1]
  fn <- cm$table[2]
  fp <- cm$table[3]
  tp <- cm$table[4]
  
  return(cutpointr::youden(tp = tp, fp = fp, tn = tn, fn = fn))
} # Function that calculates the Youden's Index
autoThreshold <- function(roc.obj){
  th.sum <- roc.obj$sensitivities + roc.obj$specificities
  th.index <- which(th.sum %in% max(th.sum))
  
  thrs <- roc.obj$thresholds[th.index]
  
  return(thrs)
}
sampleDataset <- function(x, percent = 0.8){ # Creates a training and a testing partition based on the percent variable (def: 0.8)
  if(is.data.frame(x)){
    sqty      <- round(percent*nrow(x))
    
    finish <- 0
    while(finish == 0){
      indx      <- sample(nrow(x),sqty)
      sample    <- x[indx,]
      remainder <- x[-indx,]
      if(nlevels(as.factor(sample$response)) == 2 && nlevels(as.factor(remainder$response)) == 2){
        finish <- 1
      }
    }
    # indx      <- sample(nrow(x),sqty)
    # sample    <- x[indx,]
    # remainder <- x[-indx,]

    return(list(train = sample, test = remainder))
  }
  else if(is.vector(x)){
    sqty   <- round(percent*length(x))
    sample <- x[sample(nrow(x),sqty),]
    return(sample)
  }
}
buildModel <- function(data, train.data, test.data, xi = 3, yi = 2, family = NULL){
  
  if(!missing(train.data) && !missing(test.data)){
    
    f <- reformulate(termlabels = names(train.data)[xi], response = names(train.data)[yi])

    xname <- names(train.data)[xi]
    yname <- names(train.data)[yi]
    
    m <- glm(f, data = train.data, family = "binomial") # Use 'glm()' from native R to create a model of f(x) = an + b
    pp <- predict(m, newdata = test.data, type = "response") # Use the fitted model to predict and generate prediction probabilities
    
    roc.obj <- pROC::roc(test.data[,yi], pp) # Create a ROC object
    
    thrs <- autoThreshold(roc.obj)
    labels <- test.data$response
    predictions <- as.numeric(pp > thrs) 
      
  } else if(!missing(data)){
    
    formula <- reformulate(termlabels = names(data)[xi], response = names(data)[yi])
    
    xname <- names(data)[xi]
    yname <- names(data)[yi]
    
    m <- glm(formula, data = data, family = "binomial") # Use 'glm()' from native R to create a model of f(x) = an + b
    pp <- predict(m, type = "response") # Use the fitted model to predict and generate prediction probabilities
    
    roc.obj <- pROC::roc(data$response, pp) # Create a ROC object
    
    thrs <- autoThreshold(roc.obj)
    labels <- data$response
    predictions <- as.numeric(pp > thrs)
  }
  
  return(list(model        = m,
              preds.probas = pp,
              roc.obj      = roc.obj,
              threshold    = thrs,
              labels       = labels,
              predictions  = predictions,
              x            = xname,
              y            = yname
  ))
}
calcMetrics <- function(predictions, labels, pp, roc.obj){
  cm <- caret::confusionMatrix(data = as.factor(predictions), reference = as.factor(labels)) # Create a confusion matrix using the 'caret' packages
  
  metrics <- list(
    predictions      = predictions,
    preds.probas     = pp,
    true.values      = labels,
    accuracy         = cm$overall[1],  # Out of all the classes, how much we predicted correctly, which must be high as possible
    precision        = caret::precision(as.factor(predictions), as.factor(labels)),
    recall           = caret::recall(as.factor(predictions), as.factor(labels)), # Sensitivity
    F1               = caret::F_meas(as.factor(predictions), as.factor(labels)),  # Sensitivity / true positive rate : It measures the proportion of actual positives that are correctly identified.
    Youden.score     = getYd(cm = cm)
  ) # Create a list of metrics
  
  preds   <- as.data.frame(do.call(cbind, metrics[1:3]))
  metrics <- as.data.frame(do.call(cbind, metrics[4:7]))
  metrics <- cbind(metrics, "AUC" = roc.obj$auc)
  metrics <- cbind(metrics, "95% CI Lower" = pROC::ci.auc(roc.obj, conf.level = 0.95)[1])
  metrics <- cbind(metrics, "95% CI Upper" = pROC::ci.auc(roc.obj, conf.level = 0.95)[2])
  metrics <- cbind(metrics, "99% CI Lower" = pROC::ci.auc(roc.obj, conf.level = 0.99)[1])
  metrics <- cbind(metrics, "99% CI Upper" = pROC::ci.auc(roc.obj, conf.level = 0.99)[2])
  metrics <- cbind(metrics, "youden" = getYd(cm))
  metrics <- round(metrics, 2)
  
  roc.data <- data.frame("sensitivities" = roc.obj$sensitivities,
                         "specificities" = roc.obj$specificities,
                         "thresholds"    = roc.obj$thresholds)
  
  c.table  <- data.frame(cm$table)
  
  return(list("predictions" = preds,
              "metrics"     = metrics,
              "roc.data"    = roc.data,
              "c.table"     = c.table))
}
plotCurves <- function(model, metrics){
  roc.obj <- model$roc.obj
  labels  <- model$labels
  pp      <- model$preds.probas
  xname   <- model$x
  yname   <- model$y
  
  roc.curve <- pROC::plot.roc(roc.obj, print.auc = TRUE, main = paste0("ROC-AUC Curve of ",yname," vs. ",xname),
                              cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3, col = "red") # plot a ROC curve based on the ROC analysis
  pr.obj <- PRROC::pr.curve(labels, pp, curve = T)
  
  pr.curve <- plot(pr.obj, main = paste0("PR-ROC Curve of ",yname," vs. ",xname),
                   cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3) # Create a Precision-Recall Curve
  
}
crossValidate <- function(data, xi = 3, yi = 2, type = "k-fold", k = 20){
  if(type == "k-fold"){
    k.folds <- replicate(k, sampleDataset(data, percent = 0.8), simplify = FALSE)
    names(k.folds) <- to("fold_", k)
    
    folds.models  <- list()
    folds.metrics <- list()
    
    i <- 1
    for(fold in k.folds){
      # print(paste0("Using ", names(k.folds)[i], " to cross-validate..."))
      train.data <- fold$train
      test.data  <- fold$test
      
      m    <- buildModel(train.data = train.data, test.data = test.data, xi = xi, yi = yi)
      r    <- calcMetrics(m$predictions, m$labels, m$preds.probas, m$roc.obj)
      
      # folds.models  <- append(folds.models, list(m))
      folds.metrics <- append(folds.metrics, list(r$metrics))
      
      # writexl::write_xlsx(bm.metrics, paste0('./',mname,'/','metrics.xlsx'))
      
      i <- i + 1
    }
    
    folds.metrics <- do.call(rbind, folds.metrics)
    final.metrics <- round(apply(folds.metrics, 2, mean, na.rm = TRUE),2)
    return(as.data.frame(t(final.metrics)))
  }
}
# ------------------ Solve requirements (make script reusable) ----------------
pkgs <- c('caret', 'pROC', 'lessR') # The 'hilldiv' package requires the installation of 
                           # cmake package on linux. Install via "sudo apt install cmake"
requirements(pkgs) # Call the function

# ------------------ Introduce CLI Arguments -----------------
args <- commandArgs(trailingOnly=TRUE)

cv    <- args[1] # Input Number of column which is the dependent variable
folds <- as.numeric(args[2]) # Either 'TRUE' or 'FALSE'. If TRUE, model is k-fold CVed.
aname <- args[3] # Specify how many 'folds' for CV.
xi    <- as.numeric(args[4]) # Specify dependent variable column
yi    <- as.numeric(args[5]) # Specify independent variable column

# cv <- TRUE
# folds <- 1000
# aname <- "test_cv"
# xi <- 3
# yi <- 2


print(paste0("CV:",cv))
print(paste0("Folds:",folds))
print(paste0("Name:",aname))
print(paste0("DV:",xi))
print(paste0("IV:",yi))

# if(length(args) != 3){
#   stop(paste0("Script requires 3 arguments being set! Detected ", length(args), " set arguments."))
# } else {
#   cv    <- args[1] # Input Number of column which is the dependent variable
#   k     <- args[2] # Either 'TRUE' or 'FALSE'. If TRUE, model is k-fold CVed.
#   aname <- args[3] # Specify how many 'folds' for CV.
#   xi    <- args[4] # Specify dependent variable column
#   yi    <- args[5] # Specify independent variable column
#   
#   print(paste0("Looking for dependent variable in column ", dep.var))
#   if(do.split){
#     print(paste0("Iterations set to ", iter, ". Running a ", iter,"-fold CV."))
#   } else{
#     print("Training a logistic regression model on the whole dataset with no CV.")
#   }
# }


# ------------------ Data Process ---------------
# Test the predictive power of the ratio of log10(Alistipes_indistinctus/Parabacteroides_distasonis)
path <- "../datasets/biomarkers.xlsx" # This is your file with your ratios
bm.data <- as.data.frame(readxl::read_excel(path))

if(checkIfBinary(bm.data[,2])){
  print("Response column looks ok!")
}

# ------------------ Logistic Regression on Selected Biomarker -------------------------
# aname <- "results"
# dir.create(aname)
# 
# k <- 1000

# mnames <- to("model_", k)
# names(folds.models) <- to("model_", k)
# names(folds.metrics) <- to("model_", k)


responseModel <- function(data, xi = 3, yi = 2, qty = 0.8, folds = 20, aname = "results", cv = FALSE, type = "k-fold"){
  dir.create(aname)
  
  if(cv){
    print(paste0("Performing a ", folds, "-fold Cross-Validation..."))
    cv.res <- crossValidate(data = bm.data, xi = xi, yi = yi, k = folds)
    writexl::write_xlsx(cv.res, paste0('./',aname,'/','cv_metrics.xlsx'))
  } else{
    # -------------------------------- Model Construction --------------------
    m <- buildModel(data, xi = xi, yi = yi)
    
    # -------------------------------- Metrics -------------------------------
    bm.metrics <- calcMetrics(m$predictions, m$labels, m$pp, m$roc.obj)
    writexl::write_xlsx(bm.metrics, paste0('./',aname,'/','metrics.xlsx'))
    
    #--------------------------------- Plots ---------------------------------
    pdf(paste0('./',aname,'/','ROC-AUC_PR_curves.pdf'))
    print(plotCurves(m))
    dev.off()
  }
}


bm.data$response

responseModel(data = bm.data, xi = xi, yi = yi, folds = folds, aname = aname, cv = cv)



# if(do.split){
#   print("Splitting the data into a train-test sections...")
#   learning.set <- sampleDataset(mb.data)
# 
# 
#   m <- glm(response ~ ratio_alistipes_parabacteroides, data = learning.set$train, family = "binomial") # Use 'glm()' from native R to create a model of f(x) = an + b
#   pp <- predict(m, newdata = learning.set$test, type = "response") # Use the fitted model to predict and generate prediction probabilities
# 
#   roc.obj <- pROC::roc(learning.set$test$response, pp) # Create a ROC object
# 
#   thrs <- autoThreshold(roc.obj)
# 
#   labels <- learning.set$test$response
# 
#   predictions <- as.numeric(pp > thrs) # Define a threshold to transform probabilits into discrete predictions (default 0.5)
# 
#   #-------------------------------- Metrics -------------------------------
#   cm <- caret::confusionMatrix(data = as.factor(predictions), reference = as.factor(labels)) # Create a confusion matrix using the 'caret' packages
# 
#   metrics <- list(
#     predictions      = predictions,
#     preds.probas     = pp,
#     true.values      = labels,
#     accuracy         = cm$overall[1],  # Out of all the classes, how much we predicted correctly, which must be high as possible
#     precision        = caret::precision(as.factor(predictions), as.factor(labels)),
#     recall           = caret::recall(as.factor(predictions), as.factor(labels)), # Sensitivity
#     F1               = caret::F_meas(as.factor(predictions), as.factor(labels)),  # Sensitivity / true positive rate : It measures the proportion of actual positives that are correctly identified.
#     Youden.score     = getYd(cm = cm)
#   ) # Create a list of metrics
# 
#   preds   <- as.data.frame(do.call(cbind, metrics[1:3]))
#   metrics <- as.data.frame(do.call(cbind, metrics[4:7]))
#   metrics <- cbind(metrics, "AUC" = roc.obj$auc)
#   metrics <- cbind(metrics, "95% CI Lower" = pROC::ci.auc(roc.obj, conf.level = 0.95)[1])
#   metrics <- cbind(metrics, "95% CI Upper" = pROC::ci.auc(roc.obj, conf.level = 0.95)[2])
#   metrics <- cbind(metrics, "99% CI Lower" = pROC::ci.auc(roc.obj, conf.level = 0.99)[1])
#   metrics <- cbind(metrics, "99% CI Upper" = pROC::ci.auc(roc.obj, conf.level = 0.99)[2])
#   metrics <- cbind(metrics, "youden" = getYd(cm))
#   metrics <- round(metrics, 2)
# 
#   roc.data <- data.frame("sensitivities" = roc.obj$sensitivities,
#                          "specificities" = roc.obj$specificities,
#                          "thresholds"    = roc.obj$thresholds)
# 
#   c.table  <- data.frame(cm$table)
# 
#   writexl::write_xlsx(list("predictions" = preds,
#                            "metrics"     = metrics,
#                            "roc.data"    = roc.data,
#                            "c.table"     = c.table), paste0('./',mname,'/','metrics.xlsx')
#   )
# 
#   #-------------- Plots --------------
#   roc.curve <- pROC::plot.roc(roc.obj, print.auc = TRUE, main = "ROC-AUC Curve of a Logistic Regression model \n for the log(Alistipes/Parabacteroides) Ratio",
#                               cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3, col = "red") # plot a ROC curve based on the ROC analysis
#   pr.obj <- PRROC::pr.curve(labels, pp, curve = T)
# 
#   pr.curve <- plot(pr.obj, main = "PR-ROC Curve of a Logistic Regression model \n for the log(Parabacterioides/Alistipes) Ratio",
#                    cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3) # Create a Precision-Recall Curve
# 
# 
#   pdf(paste0('./',mname,'/','ROC-AUC_PR_curves.pdf'))
#   print(plot(roc.curve, print.auc = TRUE, main = "ROC-AUC Curve of a Logistic Regression model \n for the log(Alistipes/Parabacteroides) Ratio",
#              cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3, col = "red"))
#   print(plot(pr.obj, main = "PR-ROC Curve of a Logistic Regression model \n for the log(Parabacterioides/Alistipes) Ratio",
#              cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3))
#   dev.off()
# 
# } else {
#   print("Skipping train-test split...")
#   
#   # -------------------------------- Model Construction --------------------
#   m <- buildModel(bm.data)
#   
#   # -------------------------------- Metrics -------------------------------
#   bm.metrics <- calcMetrics(m$predictions, m$labels, m$pp, m$roc.obj)
#   writexl::write_xlsx(bm.metrics, paste0('./',mname,'/','metrics.xlsx'))
#   
#   #--------------------------------- Plots ---------------------------------
#   pdf(paste0('./',mname,'/','ROC-AUC_PR_curves.pdf'))
#   print(plotCurves(m))
#   dev.off()
#   
# }

test_roc_analysis <- function(data, x = 3, y = 2){
  m <- glm(response ~ ratio_alistipes_parabacteroides, data = mb.data, family = "binomial") # Use 'glm()' from native R to create a model of f(x) = an + b
  pp <- predict(m, type = "response") # Use the fitted model to predict and generate prediction probabilities
  roc.obj <- pROC::roc(mb.data$response, pp) # Create a ROC object
  
  thrs <- autoThreshold(roc.obj)
  
  labels <- mb.data$response
  
  predictions <- as.numeric(pp > thrs) # Define a threshold to transform probabilits into discrete predictions (default 0.5)
  
  # -------------------------------- Metrics -------------------------------
  cm <- caret::confusionMatrix(data = as.factor(predictions), reference = as.factor(labels)) # Create a confusion matrix using the 'caret' packages
  
  metrics <- list(
    predictions      = predictions,
    preds.probas     = pp,
    true.values      = labels,
    accuracy         = cm$overall[1],  # Out of all the classes, how much we predicted correctly, which must be high as possible
    precision        = caret::precision(as.factor(predictions), as.factor(labels)),
    recall           = caret::recall(as.factor(predictions), as.factor(labels)), # Sensitivity
    F1               = caret::F_meas(as.factor(predictions), as.factor(labels)),  # Sensitivity / true positive rate : It measures the proportion of actual positives that are correctly identified.
    Youden.score     = getYd(cm = cm)
  ) # Create a list of metrics
  
  pROC::ci.auc(roc.obj, conf.level = 0.95)[1]
  
  preds   <- as.data.frame(do.call(cbind, metrics[1:3]))
  metrics <- as.data.frame(do.call(cbind, metrics[4:7]))
  metrics <- cbind(metrics, "AUC" = roc.obj$auc)
  metrics <- cbind(metrics, "95% CI Lower" = pROC::ci.auc(roc.obj, conf.level = 0.95)[1])
  metrics <- cbind(metrics, "95% CI Upper" = pROC::ci.auc(roc.obj, conf.level = 0.95)[2])
  metrics <- cbind(metrics, "99% CI Lower" = pROC::ci.auc(roc.obj, conf.level = 0.99)[1])
  metrics <- cbind(metrics, "99% CI Upper" = pROC::ci.auc(roc.obj, conf.level = 0.99)[2])
  metrics <- cbind(metrics, "youden" = getYd(cm))
  metrics <- round(metrics, 2)
  
  roc.data <- data.frame("sensitivities" = roc.obj$sensitivities,
                         "specificities" = roc.obj$specificities,
                         "thresholds"    = roc.obj$thresholds)
  
  c.table  <- data.frame(cm$table)
  
  writexl::write_xlsx(list("predictions" = preds,
                           "metrics"     = metrics,
                           "roc.data"    = roc.data,
                           "c.table"     = c.table), paste0('./',mname,'/','metrics.xlsx')
  )
  
  #-------------- Plots --------------
  roc.curve <- pROC::plot.roc(roc.obj, print.auc = TRUE, main = "ROC-AUC Curve of a Logistic Regression model \n for the log(Alistipes/Parabacteroides) Ratio",
                              cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3, col = "red") # plot a ROC curve based on the ROC analysis
  pr.obj <- PRROC::pr.curve(labels, pp, curve = T)
  
  pr.curve <- plot(pr.obj, main = "PR-ROC Curve of a Logistic Regression model \n for the log(Parabacterioides/Alistipes) Ratio",
                   cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3) # Create a Precision-Recall Curve
  
  
  pdf(paste0('./',mname,'/','ROC-AUC_PR_curves.pdf'))
  print(plot(roc.curve, print.auc = TRUE, main = "ROC-AUC Curve of a Logistic Regression model \n for the log(Alistipes/Parabacteroides) Ratio",
             cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3, col = "red"))
  print(plot(pr.obj, main = "PR-ROC Curve of a Logistic Regression model \n for the log(Parabacterioides/Alistipes) Ratio",
             cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.3))
  dev.off()
}
# writexl::write_xlsx(list("predictions" = preds,
#                          "metrics"     = metrics,
#                          "roc.data"    = roc.data,
#                          "c.table"     = c.table), paste0('./',mname,'/','metrics.xlsx')
# )







#### REYCLCE CODE ####
# dplyr::pull(bm.data, 2)
# # Transform character data into numeric (binary 0,1)
# bm.data$response[mb.data$response %in% "G"] <- 0 # Growth is coded as 0
# bm.data$response[mb.data$response %in% "R"] <- 1 # Reduction is coded as 1
# bm.data$response <- as.numeric(mb.data$response) # Ensure numbers are actually numeric data types
# names(mb.data)[3] <- "ratio_alistipes_parabacteroides" # Name the column
# mb.data$ratio_alistipes_parabacteroides <- log(mb.data$ratio_alistipes_parabacteroides) # Transform data using log function




