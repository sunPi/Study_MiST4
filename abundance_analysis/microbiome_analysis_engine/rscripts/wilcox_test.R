microbiome <- readRDS("/home/jan1/bioinf-tools/pipelines/Shared/parp-pipe/src/MesoXGBoost/outfolder/relative_0.1_150_development/robj/pipe.data.RDS")
microbiome$abun_species.relative <- microbiome$abun_species.relative[order(microbiome$abun_species.relative$response, decreasing = TRUE),]
rownames(microbiome$abun_species.relative) <- NULL

microbiome$abun_species.relative$response[1:11]
microbiome$abun_species.relative$response[12:length(microbiome$abun_species.relative$response)]

top10boosted <- select(microbiome$abun_species.relative, top10)
top10boosted[1]

i <- 1
res <- list()
n   <- list()

for(col in top10boosted){
  R <- col[1:11]
  G <- col[12:length(col)]
  
  pv <- wilcox.test(R, G)$p.value
  pv
  if(pv <= 0.05){
    res <- append(res, list(pv))
    n   <- append(n, top10[i])
  }
  i <- i + 1
  }

names(res) <- n
res
