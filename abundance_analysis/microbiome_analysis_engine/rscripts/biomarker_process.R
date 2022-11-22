path <- "../datasets/biomarkers_raw.xlsx" # This is your file with your ratios
bm.data <- as.data.frame(readxl::read_excel(path))

#### UNORGANIZED DATA PROCESSING ####
str(bm.data)

class(bm.data)
typeof(bm.data)

cor.mat <- cor(as.matrix(bm.data[,-1]))
heatmap(cor(as.matrix(bm.data[,-1])), Colv = NA, Rowv = NA)

MASS::boxcox()


# bm.data$logratio_alistipes_parabacteroides <- log(bm.data$logratio_alistipes_parabacteroides + 1)
# bm.data$cytotoxic_myeloid_ratio <- log(bm.data$cytotoxic_myeloid_ratio + 1)
# bm.data$UPD_ratio <- log(bm.data$UPD_ratio + 1)
MIST_ID    <- bm.data[,1]
response   <- bm.data[,2]
d     <- bm.data[,-c(1:2)]
log.d <- mapply(log1p, d)

log.bm.data <- data.frame(MIST_ID,response,log.d)

inv.logap <- log.bm.data$logratio_alistipes_parabacteroides * (-1)

log.bm.data$combined_marker <- inv.logap + log.bm.data$hrd_sum_marker + log.bm.data$cytotoxic_myeloid_ratio +log.bm.data$UPD_ratio


# bm.data$logratio_alistipes_parabacteroides <- log1p(bm.data$logratio_alistipes_parabacteroides)
# bm.data$cytotoxic_myeloid_ratio <- log1p(bm.data$cytotoxic_myeloid_ratio)
# bm.data$UPD_ratio <- log1p(bm.data$UPD_ratio )
# log1p(bm.data$hrd_sum_marker)

writexl::write_xlsx(log.bm.data, "../datasets/biomarkers.xlsx")

