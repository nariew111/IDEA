custom_r_lib <- Sys.getenv("IDEA_GP_R_LIBS", unset = "")
if (nzchar(custom_r_lib)) {
  .libPaths(c(strsplit(custom_r_lib, .Platform$path.sep)[[1]], .libPaths()))
}

suppressPackageStartupMessages(library(rrBLUP))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(genio))
# suppressPackageStartupMessages(library(KAML))

read_plink_bed_internal <- function(bfile) {
  # Use genio for fast reading if available
  if (requireNamespace("genio", quietly = TRUE)) {
    # genio expects path without extension for read_plink
    # but let's be safe and construct full paths
    bim_file <- paste0(bfile, ".bim")
    fam_file <- paste0(bfile, ".fam")
    bed_file <- paste0(bfile, ".bed")
    
    # Read plink files
    # read_plink returns a list with X (genotypes), bim, fam
    # X is an integer matrix (0, 1, 2, NA)
    data <- genio::read_plink(bfile)
    
    # genio returns 0/1/2/NA. 
    # Standard coding: 0=hom1, 1=het, 2=hom2.
    # We need to check if this matches the custom reader's expectation.
    # The custom reader returned: 0, 1, 2.
    # genio returns m x n matrix.
    
    X <- data$X
    fam <- data$fam
    bim <- data$bim
    
    # Convert to standard matrix if needed, but genio returns a matrix.
    # Ensure dimensions are correct (n_snp x n_ind)
    
    fam_df <- data.frame(fam=fam$fam, id=fam$id, stringsAsFactors=FALSE)
    # Adjust column names if genio uses different ones
    if("id" %in% names(fam)) {
        fam_df$id <- fam$id
    } else {
        fam_df$id <- fam$V2
    }
    
    bim_df <- data.frame(id=bim$id, stringsAsFactors=FALSE)
     if("id" %in% names(bim)) {
        bim_df$id <- bim$id
    } else {
        bim_df$id <- bim$V2
    }

    return(list(X=X, fam=fam_df, bim=bim_df))
  }

  fam_path <- paste0(bfile, ".fam")
  bim_path <- paste0(bfile, ".bim")
  bed_path <- paste0(bfile, ".bed")
  
  fam <- read.table(fam_path, header=FALSE, stringsAsFactors=FALSE)
  bim <- read.table(bim_path, header=FALSE, stringsAsFactors=FALSE)
  
  n_ind <- nrow(fam)
  n_snp <- nrow(bim)
  
  con <- file(bed_path, "rb")
  on.exit(close(con))
  
  header <- readBin(con, what="raw", n=3)
  if (length(header) != 3 ||
      as.integer(header[1]) != 108L ||
      as.integer(header[2]) != 27L ||
      as.integer(header[3]) != 1L) {
    stop("BED file is not in SNP-major mode or has invalid header")
  }
  
  bytes_per_snp <- ceiling(n_ind / 4)
  X <- matrix(NA_integer_, nrow=n_snp, ncol=n_ind)
  
  for (s in 1:n_snp) {
    buf <- readBin(con, what="raw", n=bytes_per_snp)
    if (length(buf) != bytes_per_snp) {
      stop("Unexpected end of BED file")
    }
    ibuf <- as.integer(buf)
    for (b in 1:bytes_per_snp) {
      val <- ibuf[b]
      for (k in 0:3) {
        ind <- (b - 1) * 4 + k + 1
        if (ind > n_ind) break
        code <- bitwAnd(bitwShiftR(val, 2 * k), 3L)
        geno <- switch(as.character(code),
                       "0" = 2L,
                       "1" = NA_integer_,
                       "2" = 1L,
                       "3" = 0L,
                       NA_integer_)
        X[s, ind] <- geno
      }
    }
  }
  
  fam_df <- data.frame(fam=fam$V1, id=fam$V2, stringsAsFactors=FALSE)
  bim_df <- data.frame(id=bim$V2, stringsAsFactors=FALSE)
  
  list(X=X, fam=fam_df, bim=bim_df)
}

parse_args <- function(args) {
  options <- list()
  for (i in seq(1, length(args), by=2)) {
    key <- sub("^--", "", args[i])
    value <- args[i+1]
    options[[key]] <- value
  }
  return(options)
}

args <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(args)

if (!is.null(opts$libpath)) {
  .libPaths(c(opts$libpath, .libPaths()))
}

# Helper function to ensure absolute paths
make_abs <- function(path) {
  if (is.null(path)) return(NULL)
  if (substr(path, 1, 1) == "/") return(path)
  return(normalizePath(path, mustWork=FALSE))
}

if (is.null(opts$bfile) && (is.null(opts$train_bfile) || is.null(opts$test_bfile))) {
  stop("Missing required arguments: either --bfile or both --train_bfile and --test_bfile")
}

if (is.null(opts$pheno) || is.null(opts$out)) {
  stop("Missing required arguments: --pheno and --out")
}

method <- if (!is.null(opts$method)) opts$method else "rrblup"
cat(paste("Using method:", method, "\n"), file=stderr())

bfile <- make_abs(opts$bfile)
train_bfile <- make_abs(opts$train_bfile)
test_bfile <- make_abs(opts$test_bfile)
pheno_file <- make_abs(opts$pheno)
train_file <- make_abs(opts$train)
test_file <- make_abs(opts$test)
extract_file <- make_abs(opts$extract)
weights_file <- make_abs(if (!is.null(opts$weights_file)) opts$weights_file else opts$weights)
out_file <- make_abs(opts$out)

pheno_df <- read.table(pheno_file, header=FALSE, stringsAsFactors=FALSE)
if (class(pheno_df[1,3]) == "character") {
  pheno_df <- read.table(pheno_file, header=TRUE, stringsAsFactors=FALSE)
}
colnames(pheno_df)[1:3] <- c("FID", "IID", "Pheno")

train_ids <- if(!is.null(train_file)) read.table(train_file, header=FALSE, stringsAsFactors=FALSE) else NULL
test_ids <- if(!is.null(test_file)) read.table(test_file, header=FALSE, stringsAsFactors=FALSE) else NULL

# Function to read genotype matrix
read_geno <- function(bfile, ids=NULL, extract_file=NULL) {
  if (requireNamespace("genio", quietly=TRUE)) {
    cat(paste("Reading", bfile, "using genio::read_plink\n"), file=stderr())
    X_list <- genio::read_plink(bfile)
    fam_data <- X_list$fam
    bim_data <- X_list$bim
    G <- X_list$X
  } else {
    cat(paste("Reading", bfile, "using internal reader\n"), file=stderr())
    X_list <- read_plink_bed_internal(bfile)
    fam_data <- X_list$fam
    bim_data <- X_list$bim
    G <- X_list$X
  }
  
  if (!is.null(ids)) {
    keep_indices <- which(paste(fam_data$fam, fam_data$id) %in% paste(ids[,1], ids[,2]))
  } else {
    keep_indices <- seq_len(nrow(fam_data))
  }
  
  if (!is.null(extract_file) && file.exists(extract_file)) {
    snps <- read.table(extract_file, header=FALSE, stringsAsFactors=FALSE)[,1]
    extract_indices <- which(bim_data$id %in% snps)
  } else {
    extract_indices <- seq_len(nrow(bim_data))
  }
  
  M <- t(G)
  M <- M[keep_indices, extract_indices, drop=FALSE]
  fam_data <- fam_data[keep_indices, , drop=FALSE]
  
  # Impute NA with mean
  for (j in 1:ncol(M)) {
    if (any(is.na(M[,j]))) {
      M[is.na(M[,j]), j] <- mean(M[,j], na.rm=TRUE)
    }
  }
  
  return(list(M=M, fam=fam_data, bim=bim_data[extract_indices, , drop=FALSE]))
}

if (!is.null(train_bfile) && !is.null(test_bfile)) {
  # Separate files
  cat("Loading separate train/test files...\n", file=stderr())
  res_train <- read_geno(train_bfile, train_ids, extract_file)
  res_test <- read_geno(test_bfile, test_ids, extract_file)
  
  M_train <- res_train$M
  M_test <- res_test$M
  
  train_ids_actual <- res_train$fam
  test_ids_actual <- res_test$fam
  
  # Ensure columns (SNPs) align
  # We assume they are aligned if extracted with same SNP list or if --extract is provided
  # If not, we intersect by ID
  snps_train <- res_train$bim$id
  snps_test <- res_test$bim$id
  
  common_snps <- intersect(snps_train, snps_test)
  if (length(common_snps) < length(snps_train) || length(common_snps) < length(snps_test)) {
      cat(paste("Aligning SNPs: Train", length(snps_train), "Test", length(snps_test), "Common", length(common_snps), "\n"), file=stderr())
      idx_tr <- match(common_snps, snps_train)
      idx_te <- match(common_snps, snps_test)
      M_train <- M_train[, idx_tr, drop=FALSE]
      M_test <- M_test[, idx_te, drop=FALSE]
  }

  current_ids <- rbind(train_ids_actual, test_ids_actual)
  # Phenotype matching
  y_train <- pheno_df$Pheno[match(paste(train_ids_actual$fam, train_ids_actual$id), paste(pheno_df$FID, pheno_df$IID))]
  y_test <- pheno_df$Pheno[match(paste(test_ids_actual$fam, test_ids_actual$id), paste(pheno_df$FID, pheno_df$IID))]
  
  if (!is.null(weights_file) && file.exists(weights_file)) {
      cat(paste("Applying weights from", weights_file, "\n"), file=stderr())
      weights_df <- read.table(weights_file, header=FALSE, stringsAsFactors=FALSE)
      # Assume col 1 is SNP, col 2 is weight
      # Normalize weights? Not strictly necessary for relative importance, but good for scale.
      # Align with M_train columns
      
      # M_train columns are SNPs
      # If aligned previously, fine.
      
      # We need to make sure weights are aligned with common_snps
      # common_snps is the column names for M_train/M_test if we aligned them.
      # Wait, read_geno returns M with no column names set explicitly in the code above (it's just a matrix).
      # But we have snps_train and snps_test. And we subsetted using idx_tr/idx_te.
      
      # Let's set colnames to be safe
      colnames(M_train) <- common_snps
      colnames(M_test) <- common_snps
      
      w_idx <- match(common_snps, weights_df$V1)
      
      # Fill missing weights with 1? or 0? 1 is safer default (no weighting).
      w_vec <- rep(1.0, length(common_snps))
      valid_w <- !is.na(w_idx)
      w_vec[valid_w] <- weights_df$V2[w_idx[valid_w]]
      
      # Apply sqrt(weight) scaling
      # u* ~ N(0, I), u = D^0.5 u*
      # y = Z D^0.5 u* + e
      # So Z_new = Z * sqrt(W)
      
      sqrt_w <- sqrt(w_vec)
      # Use matrix multiplication for scaling columns efficiently or sweep
      # M is n x p. We want to scale columns.
      # M %*% diag(sqrt_w) is slow.
      # sweep(M, 2, sqrt_w, "*") is better.
      
      M_train <- sweep(M_train, 2, sqrt_w, "*")
      M_test <- sweep(M_test, 2, sqrt_w, "*")
  }

} else {
  # Single file logic
  cat("Loading single bfile...\n", file=stderr())
  keep_ids <- rbind(train_ids, test_ids)
  res <- read_geno(bfile, keep_ids, extract_file)
  
  M <- res$M
  fam <- res$fam
  
  # Phenotype
  y <- pheno_df$Pheno[match(paste(fam$fam, fam$id), paste(pheno_df$FID, pheno_df$IID))]
  
  is_train <- paste(fam$fam, fam$id) %in% paste(train_ids[,1], train_ids[,2])
  is_test <- paste(fam$fam, fam$id) %in% paste(test_ids[,1], test_ids[,2])
  
  M_train <- M[is_train, , drop=FALSE]
  y_train <- y[is_train]
  
  M_test <- M[is_test, , drop=FALSE]
  y_test <- y[is_test]
  
  train_ids_actual <- fam[is_train, ]
  test_ids_actual <- fam[is_test, ]

  if (!is.null(weights_file) && file.exists(weights_file)) {
      cat(paste("Applying weights from", weights_file, "\n"), file=stderr())
      weights_df <- read.table(weights_file, header=FALSE, stringsAsFactors=FALSE)
      snps <- res$bim$id
      w_idx <- match(snps, weights_df$V1)
      
      w_vec <- rep(1.0, length(snps))
      valid_w <- !is.na(w_idx)
      w_vec[valid_w] <- weights_df$V2[w_idx[valid_w]]
      
      sqrt_w <- sqrt(w_vec)
      
      # Apply to M_train and M_test
      # Note: M_train and M_test are subsets of M, so columns are same
      M_train <- sweep(M_train, 2, sqrt_w, "*")
      M_test <- sweep(M_test, 2, sqrt_w, "*")
  }
}


if (any(is.na(y_train))) {
  keep_train <- !is.na(y_train)
  M_train <- M_train[keep_train, , drop=FALSE]
  y_train <- y_train[keep_train]
  train_ids_actual <- train_ids_actual[keep_train, ]
  cat(paste("Removed", sum(!keep_train), "training individuals due to missing phenotype.\n"), file=stderr())
}

pred_test <- NULL

if (method == "rrblup") {
  # rrBLUP (Ridge Regression)
  # Using mixed.solve or simple ridge regression on SNPs
  # Since X is n_snp x n_ind, we transpose it to n_ind x n_snp
  # mixed.solve solves y = Xb + u + e
  
  # We use kinship BLUP or RR-BLUP. Here RR-BLUP is faster for p >> n if implemented right,
  # but rrBLUP package mixed.solve is generally robust.
  # To use mixed.solve for RR-BLUP: K = NULL, Z = M_train
  
  # Standardize M
  # M_train <- scale(M_train) # optional, mixed.solve centers it
  
  ans <- mixed.solve(y=y_train, Z=M_train, K=NULL, SE=FALSE, return.Hinv=FALSE)
  marker_effects <- ans$u
  beta <- ans$beta
  
  # Predict
  pred_test <- as.vector(M_test %*% marker_effects) + as.vector(beta)
  
} else if (method == "lasso") {
  # Lasso using glmnet
  fit <- cv.glmnet(x=M_train, y=y_train, family="gaussian", alpha=1, type.measure="mse")
  pred_test <- predict(fit, newx=M_test, s="lambda.min")
  pred_test <- as.vector(pred_test)

} else if (method == "kaml") {
  suppressPackageStartupMessages(library(KAML))
      # KAML requires bigmemory format, so we can't use X directly if it's not file-backed.
      # But KAML.Data() can take PLINK files.
  kaml_prefix <- paste0(out_file, "_kaml_tmp")
  kaml_dir <- dirname(kaml_prefix)
  kaml_base <- basename(kaml_prefix)
  
  if (!dir.exists(kaml_dir)) dir.create(kaml_dir, recursive = TRUE)
  
  old_wd <- getwd()
  setwd(kaml_dir)
  
  KAML.Data(bfile=bfile, out=kaml_base)
  
  setwd(old_wd)
  
  desc_path <- file.path(kaml_dir, paste0(kaml_base, ".desc"))
  bin_path  <- file.path(kaml_dir, paste0(kaml_base, ".bin"))
  if (!file.exists(desc_path) || !file.exists(bin_path)) {
    desc_candidates <- list.files(kaml_dir, pattern="\\.desc$", full.names=FALSE)
    if (length(desc_candidates) > 0) {
      kaml_base <- sub("\\.desc$", "", basename(desc_candidates[1]))
      desc_path <- file.path(kaml_dir, paste0(kaml_base, ".desc"))
      bin_path  <- file.path(kaml_dir, paste0(kaml_base, ".bin"))
    }
  }
  
  kaml_pheno_file <- paste0(out_file, "_kaml_pheno.txt")
  
  df_train <- data.frame(Taxa=train_ids_actual$id, Trait=y_train)
  df_test <- data.frame(Taxa=test_ids_actual$id, Trait=NA)
    
  plink_fam_file <- paste0(bfile, ".fam")
  kaml_fam <- read.table(plink_fam_file, header=FALSE, stringsAsFactors=FALSE)
  kaml_ids <- kaml_fam$V2
    
  full_pheno <- rep(NA, length(kaml_ids))
  names(full_pheno) <- kaml_ids
    
  matches <- match(df_train$Taxa, kaml_ids)
  valid_matches <- !is.na(matches)
  full_pheno[matches[valid_matches]] <- df_train$Trait[valid_matches]
    
  kaml_pheno_df <- data.frame(kaml_ids, full_pheno)
  write.table(kaml_pheno_df, kaml_pheno_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    
  kaml_out_prefix <- paste0(out_file, "_kaml_res")
    
  old_wd <- getwd()
  setwd(kaml_dir)
    
  tryCatch({
      if (!file.exists(paste0(kaml_base, ".desc")) || !file.exists(paste0(kaml_base, ".bin"))) {
        desc_candidates <- list.files(".", pattern="\\.desc$", full.names=FALSE)
        bin_candidates  <- list.files(".", pattern="\\.bin$",  full.names=FALSE)
        base_candidates <- intersect(sub("\\.desc$", "", desc_candidates), sub("\\.bin$", "", bin_candidates))
        if (length(base_candidates) > 0) {
          kaml_base <- base_candidates[1]
        } else {
          stop("KAML.Data bigmemory (.desc/.bin) not found")
        }
      }
      
      gbase <- sub("\\.geno$", "", kaml_base)
      res <- KAML(pfile=kaml_pheno_file, gfile=gbase, file.output=FALSE, cpu=1)
      
      pred_test_local <- NULL
      if (!is.null(res$predicted)) {
        preds <- res$predicted
        test_rows <- preds[preds$Taxa %in% test_ids_actual$id, , drop=FALSE]
        m <- match(test_ids_actual$id, test_rows$Taxa)
        test_rows <- test_rows[m, , drop=FALSE]
        col_name <- "Predicted"
        if (!col_name %in% colnames(test_rows)) {
          nums <- sapply(test_rows, is.numeric)
          cnames <- colnames(test_rows)
          possibles <- cnames[nums & cnames != "Trait"]
          if (length(possibles) > 0) {
            col_name <- possibles[length(possibles)]
          }
        }
        pred_test_local <- as.vector(test_rows[[col_name]])
      }
      
      if (is.null(pred_test_local) || all(is.na(pred_test_local))) {
        stop("KAML returned empty predictions")
      }
      pred_test <<- pred_test_local
      
    }, error = function(e) {
      cat(paste("KAML Error:", e$message, "\n"), file=stderr())
      pred_test <<- rep(NA, length(y_test))
    })

    setwd(old_wd)
    
} else if (method %in% c("bayesb", "bayesc")) {
  if (!requireNamespace("BGLR", quietly=TRUE)) {
      stop("BGLR package not installed")
  }
  suppressPackageStartupMessages(library(BGLR))
  
  nIter <- 3000
  burnIn <- 500
  
  model_type <- if (method == "bayesb") "BayesB" else "BayesC"
  cat(paste("Running", model_type, "with BGLR...\n"), file=stderr())
  
  ETA <- list(list(X=M_train, model=model_type))
  
  fit <- BGLR(y=y_train, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  
  b_hat <- fit$ETA[[1]]$b
  mu_hat <- fit$mu
  
  pred_test <- as.vector(M_test %*% b_hat) + mu_hat

} else if (method == "rf") {
  if (!requireNamespace("ranger", quietly=TRUE)) {
      stop("ranger package not installed")
  }
  suppressPackageStartupMessages(library(ranger))
  
  cat("Running Random Forest (ranger)...\n", file=stderr())
  
  df_train <- as.data.frame(M_train)
  df_train$y_pheno_target <- y_train
  
  fit <- ranger(y_pheno_target ~ ., data=df_train, num.trees=300, write.forest=TRUE)
  
  df_test <- as.data.frame(M_test)
  pred <- predict(fit, data=df_test)
  pred_test <- pred$predictions

} else if (method == "svm") {
  if (!requireNamespace("e1071", quietly=TRUE)) {
      stop("e1071 package not installed")
  }
  suppressPackageStartupMessages(library(e1071))
  
  cat("Running SVM (e1071)...\n", file=stderr())
  
  fit <- svm(x=M_train, y=y_train, type="eps-regression", kernel="radial")
  pred_test <- predict(fit, newdata=M_test)
  pred_test <- as.vector(pred_test)

} else {
  stop(paste("Unknown method:", method))
}

acc <- cor(pred_test, y_test, use="complete.obs")
write(acc, file=out_file)

if (!is.null(opts$save_preds) && opts$save_preds != "0") {
  pred_file <- if (!is.null(opts$pred_out)) opts$pred_out else paste0(out_file, ".pred")
  out_df <- data.frame(FID=test_ids_actual$fam, IID=test_ids_actual$id, Pred=pred_test, True=y_test)
  write.table(out_df, pred_file, row.names=FALSE, quote=FALSE, sep="\t")
}
