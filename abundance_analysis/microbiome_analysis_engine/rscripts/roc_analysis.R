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
    indx      <- sample(nrow(x),sqty)
    sample    <- x[indx,]
    remainder <- x[-indx,]
    return(list(train = sample, test = remainder))
  }
  else if(is.vector(x)){
    sqty   <- round(percent*length(x))
    sample <- x[sample(nrow(x),sqty),]
    return(sample)
  }
}
# ------------------ Solve requirements (make script reusable) ----------------
pkgs <- c('caret', 'pROC') # The 'hilldiv' package requires the installation of 
                           # cmake package on linux. Install via "sudo apt install cmake"
requirements(pkgs) # Call the function

    
#------------ Data Process ---------------
# Test the predictive power of the ratio of log10(Alistipes_indistinctus/Parabacteroides_distasonis)
path <- "../datasets/ratio.xlsx" # This is your file with your ratios
mb.data <- readxl::read_excel(path)

# Transform character data into numeric (binary 0,1)
mb.data$response[mb.data$response %in% "G"] <- 0 # Growth is coded as 0
mb.data$response[mb.data$response %in% "R"] <- 1 # Reduction is coded as 1
mb.data$response <- as.numeric(mb.data$response) # Ensure numbers are actually numeric data types
names(mb.data)[3] <- "ratio_alistipes_parabacteroides" # Name the column

mb.data$ratio_alistipes_parabacteroides <- log(mb.data$ratio_alistipes_parabacteroides) # Transform data using log function

# ------------------ Logistic Regression on Log(A/P) -------------------------
args <- commandArgs(trailingOnly=TRUE)

iter     <- args[1]
do.split <- args[2]

mname <- paste0("model_",iter)
dir.create(mname)

if(do.split){
  print("Splitting the data into a train-test sections...")
  learning.set <- sampleDataset(mb.data)
  
  
  m <- glm(response ~ ratio_alistipes_parabacteroides, data = learning.set$train, family = "binomial") # Use 'glm()' from native R to create a model of f(x) = an + b
  pp <- predict(m, newdata = learning.set$test, type = "response") # Use the fitted model to predict and generate prediction probabilities
  
  roc.obj <- pROC::roc(learning.set$test$response, pp) # Create a ROC object
  
  thrs <- autoThreshold(roc.obj)
  
  labels <- learning.set$test$response
  
  predictions <- as.numeric(pp > thrs) # Define a threshold to transform probabilits into discrete predictions (default 0.5)
  
  #-------------------------------- Metrics -------------------------------
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
  
} else {
  print("Skipping train-test split...")
  
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










