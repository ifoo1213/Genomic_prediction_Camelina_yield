library(vcfR)
library(data.table)

# Input VCF
vcf_file <- "Camelina.imputated.vcf"

# Read VCF
vcf <- read.vcfR(vcf_file, verbose = FALSE)


# Extract genotypes (GT field only)
geno <- extract.gt(vcf, element = "GT")

# Convert to numeric (0,1,2)
geno_num <- apply(geno, 2, function(x) {
  x <- gsub("0/0", "0", x)
  x <- gsub("0/1", "1", x)
  x <- gsub("1/0", "1", x)   # just in case
  x <- gsub("1/1", "2", x)
  x <- gsub("\\./\\.", NA, x)  # missing
  as.numeric(x)
})


# Now geno_num is SNPs × individuals
geno_num <- t(geno_num)  # transpose to individuals × SNPs

# Add IDs
rownames(geno_num) <- colnames(geno)   # sample IDs
colnames(geno_num) <- vcf@fix[,3]      # SNP IDs (the ID column in VCF)

# Save to CSV
write.csv(geno_num, "camelina_genotypes.csv", row.names = TRUE)
