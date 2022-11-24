# ------------------ Header ----------------------------------------------------
# sessionInfo()
'ROC Analysis script takes input data as .csv files, which has to have columns 
formated such as 1 sample id, 2 independent variable, 3...n dependent variable(s).

Usage:
  roc_analysis.R --dataset=<FILE_PATH> [options] 
  roc_analysis.R --version

Options:
  -c --cv                              [Default: FALSE]
  -x --xindex=<RESPONSE_COLUMN_INDEX>  [Default: 3]
  -y --yindex=<DEPENDENT_COLUMN_INDEX> [Default: 2]
  -n --aname=<EXPERIMENT_FOLDER_NME>
  -f --folds=<NO_OF_FOLDS>             [Default: 5]
  -h --help     Show this screen.
  --version     Show version.

' -> doc

library(docopt)
arguments <- docopt(doc)
# print(arguments)

path  <- arguments$dataset
xi    <- as.integer(arguments$xindex)
yi    <- as.integer(arguments$yindex)
aname <- arguments$aname
do.cv <- arguments$cv
folds <- as.integer(arguments$folds)

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

    m <- glm(f, data = train.data, family = "binomial")      # Use 'glm()' from native R to create a model of f(x) = an + b
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
responseModel <- function(data, xi = 3, yi = 2, qty = 0.8, folds = 5, aname, cv = FALSE, type = "k-fold"){
  if(cv){
    print(paste0("Performing a ", folds, "-fold Cross-Validation..."))
    cv.res <- crossValidate(data = bm.data, xi = xi, yi = yi, k = folds)
    if(!is.null(aname)){
      dir.create(aname)
      writexl::write_xlsx(cv.res, paste0('./',aname,'/','cv_metrics.xlsx'))
    } else{
      writexl::write_xlsx(cv.res, paste0('./cv_metrics.xlsx'))
    }
  } else{
    # -------------------------------- Model Construction 
    m <- buildModel(data, xi = xi, yi = yi)
    
    # -------------------------------- Metrics 
    bm.metrics <- calcMetrics(m$predictions, m$labels, m$pp, m$roc.obj)
    
    if(!is.null(aname)){
      dir.create(aname)
      writexl::write_xlsx(bm.metrics, paste0('./',aname,'/','metrics.xlsx'))
    } else{
      writexl::write_xlsx(bm.metrics, paste0('./metrics.xlsx'))
    }
    
    #--------------------------------- Plots 
    if(!missing(aname)){
      pdf(paste0('./',aname,'/','ROC-AUC_PR_curves.pdf'))
      print(plotCurves(m))
      dev.off()
    } else{
      pdf(paste0('./ROC-AUC_PR_curves.pdf'))
      print(plotCurves(m))
      dev.off()
    }
  }
}
# ------------------ Solve requirements (make script reusable) ----------------
pkgs <- c('caret', 'pROC', 'lessR')
requirements(pkgs)

# ------------------ Data Process ---------------
bm.data <- as.data.frame(readxl::read_excel(path))

if(checkIfBinary(bm.data[,2])){
  print("Response column looks ok!")
} else{
  stop("Data must in a binary form!")
}

# ------------------ Logistic Regression on Selected Biomarker -------------------------
# k <- 1000
# mnames <- to("model_", k)
# names(folds.models) <- to("model_", k)

responseModel(data = bm.data, xi = xi, yi = yi, folds = folds, aname = aname, cv = do.cv)
