##### This script is used to convert an unpolarized VCF file to polarized SNAPP format. The VCF file is not directly loaded, but rather the key information is printed into separate files before running this script;

# Author use only: GitHub version derived from custom KL19 script vcf_to_SNAPP_KL19_AKBCWA_maf05_10kbp_polarized.R, 2025/12/02;

#####

### Prerequisites:

# User has already printed key information from the VCF file to create the .geno and .samples input files;

# The .geno file contains the genotypes and can be generated as follows using bcftools:

# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' KL19_AKBCWA_forPsi_maf05_10kbp.vcf > KL19_AKBCWA_forPsi_maf05_10kbp.geno

# The .samples file contains the sample  names and can be generated as follows using bcftools:

# bcftools query -l KL19_AKBCWA_forPsi_maf05_10kbp.vcf > KL19_AKBCWA_forPsi_maf05_10kbp.samples

# User has already generated the ancestral and derived allele defintions (ancestralDerived_GERP_KL19.txt) using the get_ancestralDerived_GERP.R script;

#####

# load the genotypes;

df <- read.table("KL19_AKBCWA_forPsi_maf05_10kbp.geno", sep = "\t", header = F, stringsAsFactors = F);

# load the sample names and append to df;

samples <- read.table("KL19_AKBCWA_forPsi_maf05_10kbp.samples", sep = "\t", header = F, stringsAsFactors = F);

colnames(df) <- c("chr", "pos", "ref", "alt", samples$V1);

# load the ancestral-derived definitions;

ancDer <- read.table("ancestralDerived_GERP_KL19.txt", sep = "\t", header = T, stringsAsFactors = F);

##########

# add the ancestral-derived definitions;

df$order <- 1:nrow(df);

df <- merge(df, ancDer[ , c("chr", "pos", "ref", "alt", "ancestral", "derived")], all.x = T);

df <- df[order(df$order), ]; # reorder;

# visually confirm the ref/alt and anc/der always match one another;

unique(df[ , c("ref", "alt", "ancestral", "derived")]);

# polarize the SNPs relative to the ancestral-derived definitions;

df$refIsAncestral <- df$ref == df$ancestral;

df[df == "0/1"] <- "1";
df[df == "./."] <- "?";

df[df == "0/0" & df$refIsAncestral == T] <- "0";
df[df == "0/0" & df$refIsAncestral == F] <- "2";

df[df == "1/1" & df$refIsAncestral == T] <- "2";
df[df == "1/1" & df$refIsAncestral == F] <- "0";

# check derived allele frequency by individual;

daf <- apply(df, 2, function(x) {(sum(x == "2")*2 + sum(x == "1")) / (2*sum(x != "?"))});

daf[order(daf)];

# reformat and retain only columns necessary for SNAPP format;

snapp <- df[ , !(colnames(df) %in% c("chr", "pos", "ref", "alt", "ancestral", "derived", "refIsAncestral", "order"))];

snapp <- t(snapp);

snapp[1:10, 1:10];

#

write.table(snapp, "KL19_AKBCWA_forPsi_maf05_10kbp_polarized.snapp", sep = ",", row.names = T, col.names = F, quote = F);
