# Genomic prediction pipeline for Camelina (yield)
# Models: Bayesian regression (BGLR: BRR, BayesB, BayesC), Ridge regression (glmnet), rrBLUP, Random Forest (ranger)
# Inputs:
#  - phenotype CSV: columns: ID, Yield (or change names below)
#  - genotype VCF file: samples must match phenotype IDs (or provide mapping)
# Usage: edit file paths below and run in R (R >= 4.0)

# ---------------------------
# Install required packages (run once)
# ---------------------------
# install.packages(c('vcfR','glmnet','ranger','caret','data.table','tidyverse'))
# if (!requireNamespace('BGLR', quietly = TRUE)) install.packages('BGLR', repos='https://cloud.r-project.org')
# install.packages('rrBLUP')

# ---------------------------
# Load libraries
# ---------------------------
library(vcfR)
library(data.table)
library(glmnet)
library(ranger)
library(BGLR)
library(caret)
library(tidyverse)
library(rrBLUP)

# ---------------------------
# Parameters - edit these
# ---------------------------
phen_file <- 'phenotypes.csv'
vcf_file  <- 'genotypes.vcf.gz'
pheno_id_col <- 'ID'
pheno_trait_col <- 'Yield'
n_folds <- 5
n_repeats <- 3
seed <- 42

# ---------------------------
# Read phenotype
# ---------------------------
pheno <- fread(phen_file)
colnames(pheno)[colnames(pheno)==pheno_id_col] <- 'ID'
colnames(pheno)[colnames(pheno)==pheno_trait_col] <- 'trait'
pheno <- pheno %>% dplyr::select(ID, trait)

# ---------------------------
# Convert VCF to genotype matrix (0/1/2)
# ---------------------------
vcf <- read.vcfR(vcf_file, verbose=FALSE)
gt <- extract.gt(vcf, element = 'GT', as.numeric = FALSE)
gt_t <- t(gt)

convert_gt_row <- function(gt_vec) {
  out <- sapply(gt_vec, function(g) {
    if (is.na(g) || g==".") return(NA_real_)
    g0 <- strsplit(g, ':')[[1]][1]
    g0 <- gsub('|','/', g0, fixed=TRUE)
    if (g0 %in% c('./.','.|.','.', '.')) return(NA_real_)
    alleles <- strsplit(g0, '/')[[1]]
    a_num <- suppressWarnings(as.numeric(alleles))
    if (any(is.na(a_num))) return(NA_real_)
    return(sum(a_num, na.rm=TRUE))
  })
  return(as.numeric(out))
}

message('Converting genotype strings to allele counts...')
G_list <- apply(gt_t, 1, convert_gt_row)
G <- do.call(rbind, G_list)
colnames(G) <- colnames(gt)
rownames(G) <- rownames(gt_t)

common_ids <- intersect(pheno$ID, rownames(G))
pheno2 <- pheno %>% filter(ID %in% common_ids) %>% arrange(ID)
G2 <- G[pheno2$ID, , drop=FALSE]

marker_means <- colMeans(G2, na.rm=TRUE)
inds_na <- which(is.na(G2), arr.ind = TRUE)
if (nrow(inds_na) > 0) {
  for (i in seq_len(nrow(inds_na))) {
    r <- inds_na[i,1]; c <- inds_na[i,2]
    G2[r,c] <- marker_means[c]
  }
}

G_scaled <- scale(G2, center = TRUE, scale = TRUE)
G_scaled <- as.matrix(G_scaled)

# ---------------------------
# Cross-validation setup
# ---------------------------
set.seed(seed)
cv_results <- list(ridge=c(), bayesBRR=c(), bayesB=c(), bayesC=c(), rrblup=c(), rf=c())

for (rep in seq_len(n_repeats)) {
  folds <- createFolds(pheno2$trait, k = n_folds, list = TRUE, returnTrain = FALSE)
  for (f in seq_along(folds)) {
    test_idx <- folds[[f]]
    train_idx <- setdiff(seq_len(nrow(pheno2)), test_idx)

    y_train <- pheno2$trait[train_idx]
    y_test  <- pheno2$trait[test_idx]
    X_train <- G_scaled[train_idx, , drop=FALSE]
    X_test  <- G_scaled[test_idx, , drop=FALSE]

    # ---- Ridge regression (glmnet) ----
    cv_glmnet <- cv.glmnet(x=X_train, y=y_train, alpha=0, nfolds=5)
    pred_ridge <- predict(cv_glmnet, newx=X_test, s='lambda.min')
    cv_results$ridge <- rbind(cv_results$ridge, c(cor=cor(as.numeric(pred_ridge), y_test),
                                                  rmse=sqrt(mean((as.numeric(pred_ridge)-y_test)^2))))

    # ---- Bayesian Ridge Regression (BRR) ----
    ETA <- list(list(X=G2[train_idx, , drop=FALSE], model='BRR'))
    fm <- BGLR(y=y_train, ETA=ETA, nIter=15000, burnIn=5000, verbose=FALSE)
    pred_brr <- as.numeric(G2[test_idx, , drop=FALSE] %*% fm$ETA[[1]]$b + fm$mu)
    cv_results$bayesBRR <- rbind(cv_results$bayesBRR, c(cor=cor(pred_brr,y_test), rmse=sqrt(mean((pred_brr-y_test)^2))))

    # ---- BayesB ----
    ETA_B <- list(list(X=G2[train_idx, , drop=FALSE], model='BayesB'))
    fm_B <- BGLR(y=y_train, ETA=ETA_B, nIter=15000, burnIn=5000, verbose=FALSE)
    pred_B <- as.numeric(G2[test_idx, , drop=FALSE] %*% fm_B$ETA[[1]]$b + fm_B$mu)
    cv_results$bayesB <- rbind(cv_results$bayesB, c(cor=cor(pred_B,y_test), rmse=sqrt(mean((pred_B-y_test)^2))))

    # ---- BayesC ----
    ETA_C <- list(list(X=G2[train_idx, , drop=FALSE], model='BayesC'))
    fm_C <- BGLR(y=y_train, ETA=ETA_C, nIter=15000, burnIn=5000, verbose=FALSE)
    pred_C <- as.numeric(G2[test_idx, , drop=FALSE] %*% fm_C$ETA[[1]]$b + fm_C$mu)
    cv_results$bayesC <- rbind(cv_results$bayesC, c(cor=cor(pred_C,y_test), rmse=sqrt(mean((pred_C-y_test)^2))))

    # ---- rrBLUP ----
    kinship <- A.mat(G2[train_idx, , drop=FALSE])
    fit_rr <- kin.blup(data=data.frame(ID=pheno2$ID[train_idx], trait=y_train),
                       geno='ID', pheno='trait', K=kinship)
    pred_rr <- A.mat(G2[test_idx, , drop=FALSE], impute.method='mean') %*% fit_rr$u + mean(y_train)
    cv_results$rrblup <- rbind(cv_results$rrblup, c(cor=cor(pred_rr,y_test), rmse=sqrt(mean((pred_rr-y_test)^2))))

    # ---- Random Forest ----
    df_train <- as.data.frame(X_train); df_train$trait <- y_train
    rf_fit <- ranger(trait~., data=df_train, num.trees=500)
    pred_rf <- predict(rf_fit, data=as.data.frame(X_test))$predictions
    cv_results$rf <- rbind(cv_results$rf, c(cor=cor(pred_rf,y_test), rmse=sqrt(mean((pred_rf-y_test)^2))))

    message(sprintf('Rep %d Fold %d done.', rep, f))
  }
}

# ---------------------------
# Summarize CV results
# ---------------------------
summarize <- function(mat) {
  mat <- as.matrix(mat)
  data.frame(metric=colnames(mat), mean=apply(mat,2,mean), sd=apply(mat,2,sd))
}

print('Ridge:'); print(summarize(cv_results$ridge))
print('Bayes BRR:'); print(summarize(cv_results$bayesBRR))
print('BayesB:'); print(summarize(cv_results$bayesB))
print('BayesC:'); print(summarize(cv_results$bayesC))
print('rrBLUP:'); print(summarize(cv_results$rrblup))
print('Random Forest:'); print(summarize(cv_results$rf))

# End of script
