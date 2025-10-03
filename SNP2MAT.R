library(vcfR)
library(data.table)

# Input VCF
vcf_file <- "genotypes.vcf.gz"

# Read VCF
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Extract genotype calls
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
gt_t <- t(gt)  # samples x markers

# Function to convert GT to numeric (0/1/2/NA)
convert_gt_row <- function(gt_vec) {
  out <- sapply(gt_vec, function(g) {
    if (is.na(g) || g == ".") return(NA_real_)
    g0 <- strsplit(g, ":")[[1]][1]
    g0 <- gsub("|", "/", g0, fixed = TRUE)
    if (g0 %in% c("./.", ".|.", ".", ".")) return(NA_real_)
    alleles <- strsplit(g0, "/")[[1]]
    a_num <- suppressWarnings(as.numeric(alleles))
    if (any(is.na(a_num))) return(NA_real_)
    return(sum(a_num, na.rm = TRUE))
  })
  return(as.numeric(out))
}

# Convert all samples
G_list <- apply(gt_t, 1, convert_gt_row)
G <- do.call(rbind, G_list)
colnames(G) <- colnames(gt)   # marker IDs
rownames(G) <- rownames(gt_t) # sample IDs

# Save to CSV
write.csv(G, file = "genotypes_matrix.csv", row.names = TRUE)
