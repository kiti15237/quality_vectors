## ---- setup ----
setwd("C:/Users/ctata/Documents/Lab/quality_vectors/microbiome_plvm-master/src/antibiotics-study")
library("rstan")
library("reshape2")
library("tidyverse")
library("stringr")
library("phyloseq")
library("ggscaffold")
library("feather")
theme_set(min_theme(list(text_size = 7, subtitle_size = 9)))
source("./posterior_check_funs.R")
dir.create("../../data/fits/", recursive = TRUE)
dir.create("../../data/figure-input/", recursive = TRUE)
dir.create("../../doc/figure/", recursive = TRUE)
set.seed(11242016)

## ---- get-data ----
setwd("C:/Users/ctata/Documents/Lab/quality_vectors/microbiome_plvm-master/src/AG")
abt <- get(load("../../data/antibiotics-study/abt.rda"))
ag <- readRDS("ag_data.rds")
abt <- abt %>%
  subset_samples(ind == argv$subject) #The 54 samples here refer to 54 timepoints of the same individual

releveled_sample_data <- abt %>%
  sample_data %>%
  mutate(
    condition = recode(
      condition,
      "Pre Cp" = "Pre",
      "1st Cp" = "1st Course",
      "1st WPC" = "1st Course",
      "2nd Cp" = "2nd Course",
      "2nd WPC" = "2nd Course",
      "Post Cp" = "Post"
    )
  )
rownames(releveled_sample_data) <- abt %>%
  sample_data %>%
  rownames
sample_data(abt) <- releveled_sample_data


## ---- lda ----
ag_small <- subset_samples(ag, BMI > 10 & BMI < 70)

#Severly restrict samples by sequencing depth
ag_small <- prune_samples(colSums(otu_table(ag_small)) > 13000, ag_small)
ag_small <- prune_samples(colSums(otu_table(ag_small)) < 15000, ag_small)
#Remove any taxa not appearing in those samples
ag_small <- filter_taxa(ag_small, function(abund) return(sum(abund) > 10), prune = T)

x <- t(get_taxa(ag_small))
dimnames(x) <- NULL



fitModel <- function(K){
  print(paste("Number of topics: ", K))
  stan_data <- list(
    K = K,
    V = ncol(x),
    D = nrow(x),
    n = x,
    alpha = rep(1, K),
    gamma = rep(0.5, ncol(x))
  )
  
  start_fit <- Sys.time()
  f <- stan_model(file = "../stan/lda_counts.stan")
  
  possibleError <- tryCatch(
    
    stan_fit <- vb(
      f,
      data = stan_data,
      output_samples = 50,
      eta = 1,
      adapt_engaged = FALSE
    ),
    error=function(e) e
    
  )
  
  if(inherits(possibleError, "error")){
    print(paste(K, " topics did not converge"))
    next
  }
  else{
    cat(sprintf(
      "Finished in %f minutes\n",
      Sys.time() - start_fit, 4)
    )
  }
  return(stan_fit)
}
Ks <- c(2, 3, 4, 5, 6, 7, 8, 9)
for(K in Ks){
  fitModel(K)
}
stan_fit <- fitModel(K = 7)



#In this scenario, each time point gets its own topics - it's as if every time point is its own sample

save(
  stan_fit,
  file = sprintf("../../data/fits/lda-%s-%s.rda", "AG", gsub("[:|| ||-]", "", Sys.time()))
)
samples <- rstan::extract(stan_fit)
rm(stan_fit)


## ---- extract_beta ----
# underlying RSV distributions
beta_logit <- samples$beta

for (i in seq_len(nrow(beta_logit))) {
  for (k in seq_len(stan_data$K)) {
    beta_logit[i, k, ] <- log(beta_logit[i, k, ])
    beta_logit[i, k, ] <- beta_logit[i, k, ] - mean(beta_logit[i, k, ])
  }
}

beta_hat <- beta_logit %>%
  melt(
    varnames = c("iterations", "topic", "rsv_ix"),
    value.name = "beta_logit"
  ) %>%
  as_data_frame()

beta_hat$rsv <- rownames(tax_table(abt))[beta_hat$rsv_ix]
taxa <- as_data_frame(tax_table(abt)@.Data)
taxa$rsv <- rownames(tax_table(abt))
taxa$Taxon_5[which(taxa$Taxon_5 == "")] <- taxa$Taxon_4[which(taxa$Taxon_5 == "")]

beta_hat <- beta_hat %>%
  left_join(taxa) %>%
  mutate(
    topic = paste("Topic", topic),
    Taxon_5 = str_extract(Taxon_5, "[^_]+")
  )

sorted_taxa <- names(sort(table(beta_hat$Taxon_5), decreasing = TRUE))
beta_hat$Taxon_5 <- factor(beta_hat$Taxon_5, levels = sorted_taxa)
beta_hat$rsv <- factor(beta_hat$rsv, levels = taxa$rsv)
































#Read in AG data and make phyloseq object
library(phyloseq)
otu_filtered_AG_07perc_feces <- read.delim("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/feces/otu_filtered_AG_07perc_feces.csv")
otu <- otu_filtered_AG_07perc_feces
colnames(otu) <- gsub("X", "", colnames(otu))

AG_mapping <- read.delim("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/AG_mapping.txt")
map <- AG_mapping

map <- map[map$BODY_SITE == "UBERON:feces", ]
map$BMI <- as.numeric(as.character(map$BMI))
map <- map[!is.na(map$BMI), ]
map <- map[map$BMI < 100, ]
map <- map[map$BMI > 1, ]
rownames(map) <- map$X.SampleID
otu <- otu[ , colnames(otu) %in% rownames(map)]
biom <- import_biom("C:/Users/ctata/Documents/Lab/quality_vectors/data/AG_new/07-taxa/07-taxa/notrim/otu_table_fecal_L3.biom")


ps <- phyloseq(otu_table(otu, taxa_are_rows = T), sample_data(map))
setwd("C:/Users/ctata/Documents/Lab/quality_vectors/microbiome_plvm-master/src/AG")
saveRDS(ps, "ag_data.rds")
