# Descriptive Statistics in R
# https://statsandr.com/blog/descriptive-statistics-in-r/#introduction
# ------------------ Define Functions
requirements <- function(pkgs){ # This function checks package requirements and install them if they are missing
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
}
load_data <- function(path){ # Function that takes path/to/the/data.xlsx in excel format and returns a loaded data frame
  df <- do.call(cbind.data.frame, readxl::read_excel(path))
  return(df)
}
replace_nan <- function(x){
  x[is.nan(x)] <- 0
  return(x)
}
compareMeans <- function(x,y,paired = FALSE){
  norm.t <- list(x,y)
  
  for(n in norm.t){
    if(length(n) < 3 || length(n) > 500){
      print(paste0("Can not do a normality test since sample size is ", length(n), ". It must be between 3 and 5000"))
      norm <- FALSE
    } else{
      if(shapiro.test(n)$p.value < 0.05){
        norm <- FALSE
        print("Detected a dataset that follows a non-normal distribution. Stopping the test...")
        break
      } else{
        norm <- TRUE
      }
    }
  }
  
  print("Succesfully tested all data for normality. All data follows normal distribution!")
  
  if(norm){
    p <- var.test(x, y)$p.value # Variance statistical test
    
    if(p < 0.05){
      print("Variances are not equal, performing a Welch Two Sample t-test.")
      welch <- t.test(x, y, var.equal = FALSE)
      return(welch)
    } else{
      print("Variances are equal, performing a standard t-test.")
      stud <- t.test(x, y, var.equal = TRUE, paired = paired)
      return(stud)
    }
  } else{
    print("Performing a Wilcox Rank-Sum test...")
    wcox <- wilcox.test(x, y, paired = paired)
    return(wcox)
  }
}
runPCA2 <- function(x, sid = NULL, groups = NULL, ellipse = TRUE){
  x <- x[, colSums(x != 0) > 0]
  pca <- prcomp(x, center = TRUE, scale. = TRUE)
  p   <- ggbiplot(pca, ellipse=ellipse, labels = NULL, groups = groups)
  return(list(pca    = pca,
              biplot = p))
}
# runPCA <- function(x, groups = NULL, ellipse = TRUE){
  # http://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining
  x <- x[, colSums(x != 0) > 0]
  pca <- prcomp(x, center = TRUE, scale. = TRUE)
  p   <- fviz_pca_biplot(pca, col.var="black", habillage=groups,
                         addEllipses=ellipse, ellipse.level=0.95, repel=TRUE, 
                         title = "Responders vs. Non-Responders PCA Biplot")
  return(list(pca    = pca,
              biplot = p))
}
run_abun_analysis <- function(data, name = "../results.xlsx"){
  results  <- list()
  dir.create(name) # Creates a new directory for results
  #-------------------- Data Process --------------------
  # Split the data into responders and non-responders
  responders    <- filter(data, data$response %in% 1)
  nonresponders <- filter(data, data$response %in% 0)
  
  sid        <- data[,1:4]
  sid.res    <- responders[,1:4]
  sid.nonres <- nonresponders[,1:4]
  
  #-------------------- Shannon Entropy Calculation --------------------
  # Calculate Entropy (H(x)) 
  # Measures the uncertainty of a probability distribution
  cat(green("Calculating Shanon Entropy...\n"))
  
  entp.res    <- apply(responders[5:ncol(responders)], 1, entropy)
  entp.nonres <- apply(nonresponders[5:ncol(nonresponders)], 1, entropy)
  
  entp.res    <- replace_nan(entp.res)
  entp.nonres <- replace_nan(entp.nonres)
  
  repets <- length(c(entp.res, entp.nonres))
  Entropy <- as.numeric(c(entp.res, entp.nonres))
  sls <- data.frame(cbind(Response=c(rep("Responders", repets/2), rep("Non-Responders", repets/2)), Entropy = as.numeric(c(entp.res, entp.nonres))))
  sls$Entropy <- as.numeric(sls$Entropy)
  
  Responders    <- entp.res
  Nonresponders <- entp.nonres
  
  pdf(paste0(name,"/","entropy_plots.pdf"), height = 11.69, width = 8.27)
  # par(mfrow=c(2,2))
  par(cex.lab=2.0) # y-axis
  par(cex.axis=2.0) # is for x-axis
  boxplot(cbind("Responders" = entp.res, "Non-Responders" = entp.nonres))
  dev.off()
  
  ggplot(data=sls, aes(x=Entropy,fill = Response)) + geom_histogram(bins = 10)
  ggsave(paste0(name,"/","entropy_stacked.pdf"))
  
  cor.test(entp.res, entp.nonres)
  
  entropy.t <- compareMeans(entp.res, entp.nonres) # Use Welch Two Sample t-test to compare means
  summary(entp.res) # Resulting means are not significantly different from each other
  summary(entp.nonres) # ...Responders have a higher diversity meaen, which is not statistically significant
  entp.res.sd    <- sd(entp.res)
  entp.nonres.sd <- sd(entp.nonres) 
  
  #-------------------- Calculate Simpsons D --------------------
  cat(green("Calculating Simpsons distribution index...\n"))
  
  simpsons.res    <- apply(responders[5:ncol(responders)], 1, simpson)
  simpsons.nonres <- apply(nonresponders[5:ncol(nonresponders)], 1, simpson)
  
  simpsons.res    <- replace_nan(simpsons.res)
  simpsons.nonres <- replace_nan(simpsons.nonres)
  
  # shapiro.test(simpsons.res) # Data appears to be normally distributed
  # shapiro.test(simpsons.nonres)
  # var.test(simpsons.res, simpsons.nonres, alternative = "two.sided") # Variances are equal
  
  # simps.t         <- t.test(simpsons.res, simpsons.nonres) # Use Welch Two Sample t-test to compare means
  simps.t         <- compareMeans(simpsons.res, simpsons.nonres) 
  simps.res.sd    <- sd(simpsons.res)
  simps.nonres.sd <- sd(simpsons.nonres)
  #-------------------- Calculate Hill numbers --------------------
  #https://www.uvm.edu/~ngotelli/manuscriptpdfs/ChaoHill.pdf
  # q = 0 indicates species richess index
  # q = 1 exponential of Shannon's entropy index
  # q = 2 the inverse of Simpson's concentration index
  cat(green("Calculating Hill's Index...\n"))
  
  resp.matrix    <- as.matrix(responders[,5:ncol(responders)])
  nonresp.matrix <- as.matrix(nonresponders[,5:ncol(nonresponders)])
  
  hill.res <- as.data.frame(Hill(resp.matrix, window = 3, alpha = 1, rasterOut=TRUE, 
                                 np = 1, na.tolerance=1.0, cluster.type = "SOCK", 
                                 debugging = FALSE))
  
  hill.nonres <- as.data.frame(Hill(nonresp.matrix, window = 3, alpha = 1, rasterOut=TRUE, 
                                    np = 1, na.tolerance=1.0, cluster.type = "SOCK", 
                                    debugging = FALSE))
  
  names(hill.res)    <- names(as.data.frame(resp.matrix))
  names(hill.nonres) <- names(as.data.frame(nonresp.matrix))
  
  hill.index.res    <- apply(hill.res, 2, mean)
  hill.index.nonres <- apply(hill.nonres, 2, mean)
  
  hill.index.res    <- replace_nan(hill.index.res)
  hill.index.nonres <- replace_nan(hill.index.nonres)
  
  repets <- length(c(hill.index.res, hill.index.nonres))
  sls    <- data.frame(cbind(Response=c(rep("Responders", repets/2), rep("Non-Responders", repets/2)), Hill_Index = as.numeric(c(hill.index.res, hill.index.nonres))))
  sls$Hill_Index <- as.numeric(sls$Hill_Index)
  
  # sls$Response[sls$Response %in% 0] <- "Nonresponders"
  # sls$Response[sls$Response %in% 1] <- "Responders"
  # print(data$SampleID)
  # print(data$response)
  
  Responders    <- hill.index.res
  Nonresponders <- hill.index.nonres
  
  pdf(paste0(name,"/","hill_ind_plots.pdf"), height = 11.69, width = 8.27)
  # par(mfrow = c(2,2))
  par(cex.lab=2.0) # y-axis
  par(cex.axis=2.0) # is for x-axis
  boxplot(cbind("Responders" = entp.res, "Non-Responders" = entp.nonres))
  dev.off()
  
  # bw <- 2 * IQR(sls$Hill_Index) / length(sls$Hill_Index)^(1/3)
  
  ggplot(data=sls, aes(x=Hill_Index,fill = Response)) + geom_histogram(bins = 5)
  ggsave(paste0(name,"/","hill_stacked.pdf"))
  
  if(length(boxplots$out) > 0){
    hill_outliers <- boxplots$out
    table <- cbind("Species Boxplot Outliers" = names(hill_outliers),"Hill Index" = hill_outliers)
    rownames(table) <- NULL
    
    pdf(paste0(name,"/","hill_outliers.pdf"), height = 11.69, width = 8.27)
    grid.table(table)
    dev.off()
  }
  
  hill.res.sd       <- sd(hill.index.res)
  hill.nonres.sd    <- sd(hill.index.nonres)
  # hill.index.wilcox <- wilcox.test(hill.index.res, hill.index.nonres, paired = FALSE) # Rank-Sum Test or Mann-Whitney U Test
  hill.index.wilcox <- compareMeans(hill.index.res, hill.index.nonres)
  #-------------------- Gamma diversity --------------------
  # Q parameter corresponds to how the Hill function calculates mean per group
  # 0 = weighted harmonic mean
  # 1 = geometric mean
  # 2 = weighted arithmetic mean
  
  # Prepare inputs for gamma
  
  input.res    <- as.matrix(responders[,5:ncol(responders)])
  input.nonres <- as.matrix(nonresponders[,5:ncol(nonresponders)])
  
  gamma.res    <- gamma_div(responders[,5:ncol(responders)], q = 0)
  gamma.nonres <- gamma_div(input.nonres, q = 0)
  
  #-------------------- PCA Analysis ----
  cat(green("Running PCA analysis...\n"))
  
  groups <- c(rep("Responders", 11), rep("Non-Responders", 11))
  # microbiome.pca <- runPCA(data[,5:ncol(data)], sid = sid, groups = groups) 
  microbiome.pca <- runPCA(data[,5:ncol(data)], groups = groups) 
  ggsave(paste0(name,"/","pca_analysis_pc1_pc2.pdf"))
  
  # # Responders
  # res.pca <- runPCA(responders[,5:ncol(responders)])
  # ggsave(paste0(name,"/","pca_analysis_responders.pdf"))
  # 
  # # Nonresponders
  # nonres.pca <- runPCA(nonresponders[,5:ncol(nonresponders)])
  # ggsave(paste0(name,"/","pca_analysis_nonresponders.pdf"))
  
  #--------------------Build the results data frame --------------------
  cat(green("Building the results dataframe...\n"))
  groups    <- rbind("responders", "nonresponders")
  entropy   <- rbind("responders" = mean(entp.res), "nonresponders" = mean(entp.nonres))
  entropy.p <- rbind("p-value"    = entropy.t$p.value, NA)
  entp.sd   <- rbind(entp.res.sd, entp.nonres.sd)
  
  simps     <- rbind("responders" = mean(simpsons.res), "nonresponders" = mean(simpsons.nonres))
  simps.p   <- rbind("p-value"    = simps.t$p.value, NA)
  simps.sd  <- rbind(simps.res.sd, simps.nonres.sd)
  
  hillind   <- rbind("responders" = mean(hill.index.res), "nonresponders" = mean(hill.index.nonres))
  hill.p    <- rbind("p-value"    = hill.index.wilcox$p.value, NA)
  hill.sd   <- rbind(hill.res.sd, hill.nonres.sd)
  
  gamma     <- rbind(gamma.res, gamma.nonres)
  
  c.names   <- c("Groups", 
                 "Shannon Entropy", "p-value", "Standard Deviation", 
                 "Simpsons Index", "p-value", "Standard Deviation",
                 "Hill Index", "p-value", "Standard Deviation",
                 "Gamma Index")
  results <- data.frame(groups, 
                        round(entropy,4), round(entropy.p,4), round(entp.sd,4), 
                        round(simps,4), round(simps.p,4), round(simps.sd,4),
                        round(hillind,4), round(hill.p,4), round(hill.sd,4),
                        gamma)
  names(results) <- c.names
  writexl::write_xlsx(results, paste0(name,"/","stats.xlsx"))
}

# ------------------ Solve requirements (make script reusable)
pkgs <- c('rasterdiv', 'agrmt', 'dplyr', 'hilldiv', 'abdiv', 'gridExtra', 'grid', 'ggbiplot', 'crayon', 'factoextra') # The 'hilldiv' package requires the installation of 
                                                                                              # cmake package on linux. Install via "sudo apt install cmake"
requirements(pkgs) # Call the function

# ------------------ Read in the Microbiome data & Prepare data ------------------

SA.data <- load_data("../asv/abun_species.absolute.xlsx") # HERE IS WHERE YOU LOAD YOUR DATASET

# Check for duplicates ----
# SA.data.duplicates <- as.data.frame(read.csv("../datasets/raw/asv_table.species.absolute.xls", sep='\t'))
# length(SA.data.duplicates$Taxonomy) # Length of original taxonomy names 
# length(unique(SA.data.duplicates$Taxonomy)) # Lenght of taxonomy names with removed duplicates
# ordered <- SA.data.duplicates[order(SA.data.duplicates$Taxonomy),]
# rownames(ordered) <- NULL
# ordered # Ordered dataframe makes it easier to spot duplicates

#-------------------- Basic Analysis -------------------------
cat(magenta("First microbiome abundance analysis run...\n"))

run_abun_analysis(SA.data, name = "species_absolute_baseline_results")

#-------------------- Boruta/XGBoost Analysis -------------------------
cat(magenta("Analysing microbiome using Boruta/XGBoost overlaping features...\n"))

boruta      <- read.csv("../aleks_boruta/bor_sel_absolute.csv")
boruta.feat <- boruta$X

xgb <- readxl::read_excel("../../MesoXGBoost/outfolder/0.05/species_absolute_0.08e_600i_0.5g/species_absolute_0.08e_600i_0.5g_imatrix.xlsx")
xgb.feat <- xgb$Feature

sel.feat <- intersect(xgb.feat,boruta.feat)
SA.data.sel <- SA.data %>% select(SampleID,best_percentage_change,order,response,sel.feat)

run_abun_analysis(SA.data.sel, name = "species_absolute_xgb&boruta_results")

#-------------------- Only Boruta -------------------------
cat(magenta("Analysing microbiome using Boruta features...\n"))

boruta.feat <- boruta.feat[-c(1,8)]

SA.data.bor <- SA.data %>% select(SampleID,best_percentage_change,order,response,boruta.feat, grep("Clostridiales_bacterium", names(SA.data)))
run_abun_analysis(SA.data.bor, name = "species_absolute_boruta_selected_results")

#-------------------- Only XGBoost -------------------------
cat(magenta("Analysing microbiome using XGBoost features...\n"))

SA.data.xgb <- SA.data %>% select(SampleID,best_percentage_change,order,response,xgb.feat)

run_abun_analysis(SA.data.xgb, name = "species_absolute_xgb_selected_results")
