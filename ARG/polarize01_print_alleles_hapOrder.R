##### This script is used prior to polarizing the sites in a .hap file for running Relate. This is an initial step that is used to print information about the ancestral and derived alleles corresponding to each site in the same order that it appears in the initial unpolarized .hap file;

# Author use only: GitHub version derived from custom KL19 script 01_print_alleles_hapOrder_KL19_v5.R, 2025/12/03;

#####

### Prerequisites:

# User has already identified ancestral and derived alleles in the genome, which can be done using the Github script get_ancestralDerived_GERP.R in the directionality_index folder. There are more columns printed than we need, for here we only need to retain five columns (chr, pos, ref, alt, ancestral_allele). The file formatted with these five columns only is as follows:

# KL19_nDNA_allFilters_rangewide_8x_noRelatives_ancDer.txt;

# User has printed the first five columns of the unpolarized .hap file and saved to a separate file. This can be done as follows:

# awk '{print $1, $2, $3, $4, $5}' KL19_forARG_v5_nonMiss90_variantCheck.hap > KL19_forARG_v5_nonMiss90_variantCheck_firstFiveColumns.txt

#####

# load the ancestral and derived alleles file;
# format is five tab-delimited columns: chr = scaffold name, pos = basepair position on the scaffold, ref = reference allele, alt = alternate allele, ancestral_allele = ancestral allele;

alleles <- read.table("KL19_nDNA_allFilters_rangewide_8x_noRelatives_ancDer.txt", sep = "\t", header = T, stringsAsFactors = F); # takes a while;

# load the first five columns of the unpolarized .hap file;

hap <- read.table("KL19_forARG_v5_nonMiss90_variantCheck_firstFiveColumns.txt", sep = " ", header = F, stringsAsFactors = F); # takes a while;
colnames(hap) <- c("chr", "uniq_id", "pos", "ref", "alt");

#####

# ensure that the alleles file has all the columns we need, including derived;

colnames(alleles)[colnames(alleles) == "ancestral_allele"] <- "ancestral";

alleles$derived <- NA;
alleles$derived[alleles$ancestral == alleles$ref] <- alleles$alt[alleles$ancestral == alleles$ref];
alleles$derived[alleles$ancestral == alleles$alt] <- alleles$ref[alleles$ancestral == alleles$alt];

# add the ancestral/derived definitions to the haplotypes by merging, and DO NOT sort the hap dataframe;

hap$order <- 1:nrow(hap);

hap2 <- merge(hap, alleles[ , c("chr", "pos", "ref", "alt", "ancestral", "derived")], sort = F); # takes a while;

nrow(hap);
nrow(hap2); # should be equal;

sum(hap$order == hap2$order);
sum(hap$order != hap2$order); # should be zero;

# just for interest's sake:
sum(hap2$ref == hap2$ancestral, na.rm = T);
sum(hap2$ref == hap2$derived, na.rm = T);

# confirm no unexpected combinations;
unique(hap2[ , c("ref", "alt", "ancestral", "derived")]);

# since there may be sites with no info (but which we still need to retain because they are present in the haplotype file and we want to print the sites in exactly the same order), we need to create three categories of a new column refIsAnc;
  # 0 ref is not ancestral (i.e., ref == derived);
  # 1 ref is ancestral (i.e., ref == ancestral);
  # -9 ref is not defined (i.e., ancestral == NA);

hap2$refIsAnc <- NA;
hap2$refIsAnc[!(is.na(hap2$ancestral)) & hap2$ref == hap2$derived] <- 0;
hap2$refIsAnc[!(is.na(hap2$ancestral)) & hap2$ref == hap2$ancestral] <- 1;
hap2$refIsAnc[is.na(hap2$ancestral)] <- -9;

unique(hap2$refIsAnc);
sum(hap2$refIsAnc == 0);
sum(hap2$refIsAnc == 1);
sum(hap2$refIsAnc == -9);

# save the output;

write.table(hap2, "alleles_withAncestralDefinitions_inHapOrder_KL19_v5.txt", sep = "\t", row.names = F, col.names = T, quote = F);
