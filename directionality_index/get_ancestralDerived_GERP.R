##### This script is used to process output from the modified gerpcol script of Taylor et al. (2024) and explicitly identify the ancestral vs. derived allele;

# Author use only: GitHub version derived from custom KL19 script get_ancestralDerived_GERP_KL19.R, 2025/12/02;

#####

### Prerequisites:

# User has already run the modified gerpcol script from Taylor et al. (2024) https://doi.org/10.1016/j.cub.2024.02.002;

# The modified gerpcol script is available from https://github.com/BeckySTaylor/Phylogenomic_Analyses/ within the gerp.tar.gz file, and the version used here was dated 2023/11/20;

#####

# get a list of all scaffolds;

autoScf <- read.table("Macpyr2_autosomalScaffolds.txt", header = F, sep = "\t", stringsAsFactors = F); # list of autosomal scaffolds with three tab-delimited columns: scaffold name, initial base-pair position of scaffold (always 1), final base-pair position of scaffold (i.e., scaffold length);

# read in all the variable sites;

df <- read.table("refAlt_alleles_KL19_nDNA_allFilters_rangewide_1x_noRelatives.txt", header = F, sep = "\t", stringsAsFactors = F); # list of variant sites in the genome, one per line, each line with four tab-delimited columns: scaffold name, bp position within scaffold, reference allele, alternative allele;
colnames(df) <- c("chr", "pos", "ref", "alt");

#

for (i in 1:nrow(autoScf)) {
  
  # counter;
  print(paste0("Reading in gerpcol for ", autoScf$V1[i]));
  
  # read in the gerpcol output;
  # the gerpcol output is from the modified gerpcol script mentioned at the beginning of this file;
  # the gerpcol output has a single file per scaffold;
  # the format is six-column tab-delimited: 0-indexed base-pair position, N parameter from gerpcol, S parameter from gerpcol, allele for outgroup 1, allele for outgroup 2, allele for outgroup 3;
  gerpcol <- read.table(paste0("gerpcol_output_fasta500/", autoScf$V1[i], "_GERP_formatted.mfa.rates"), header = T, sep = "\t", stringsAsFactors = F);
  colnames(gerpcol) <- c("pos", "N", "S", "Nerelu", "Saccja", "Lamidi"); # the focal species in this script is Macrocystis, so we have "Nerelu" or Nereocystis luetkeana as outgroup 1;
  
  # remove lines where N = -1, which indicates that there was no outgroup alignment;
  gerpcol <- gerpcol[!(gerpcol$N == -1), ];
  
  # remove the columns N and S, which we don't need anymore and were used in GERP analyses;
  gerpcol <- gerpcol[ , c("pos", "Nerelu", "Saccja", "Lamidi")];
  
  # change the gerpcol pos to be 1-indexed instead of 0-indexed, to match the format of the original VCF indexing and df;
  gerpcol$pos <- gerpcol$pos + 1;
  
  # add chromosome to gerpcol;
  gerpcol$chr <- autoScf$V1[i];
  
  # subset df to only the scaffold of interest;
  cur_df <- df[df$chr == autoScf$V1[i], ];
  
  # merge with gerpcol;
  # IMPORTANT, use "all = F, sort = F" to ensure that the order of alleles is not changed;
  cur_df <- merge(cur_df, gerpcol, all = F, sort = F);
  
  # combine into output;
  if (i == 1) {
    out <- cur_df;
  } else {
    out <- rbind(out, cur_df);
  }
  
}

#

#####

# define the ancestral and derived allele - in order for an allele to be derived, it must not be present in any of the three outgroups (but missing data is allowed, i.e., we do not need data for all three outgroups);
# if neither the ref nor the alt is non-present in any of the three outgroups, then we cannot define ancestral and derived alleles;

out$n_outgroup_alleles <- apply(out[ , c("Nerelu", "Saccja", "Lamidi")], 1, function(x) {length(unique(x[!is.na(x) & x != "N"]))});

out$unique_outgroup_allele <- NA;
out$unique_outgroup_allele[out$n_outgroup_alleles == 1] <- apply(out[out$n_outgroup_alleles == 1, c("Nerelu", "Saccja", "Lamidi")], 1, function(x) {unique(x[!is.na(x) & x != "N"])});

out$ancestral <- NA;
out$derived <- NA;

out$ancestral[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)] <- out$ref[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)];
out$ancestral[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)] <- out$alt[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)];

out$derived[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)] <- out$alt[out$unique_outgroup_allele == out$ref & !is.na(out$unique_outgroup_allele)];
out$derived[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)] <- out$ref[out$unique_outgroup_allele == out$alt & !is.na(out$unique_outgroup_allele)];

# remove sites where the ancestral and derived could not be defined (it is still NA) because the unique outgroup allele present was neither the ref nor the alt allele;
out <- out[!(is.na(out$ancestral)), ];

#####

write.table(out, "ancestralDerived_GERP_KL19.txt", sep = "\t", row.names = F, col.names = T, quote = F);

