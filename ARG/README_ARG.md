# Ancestral Recombination Graphs

### polarize01_print_alleles_hapOrder.R

This _R_ script is an initial preparatory script used to print the ancestral and derived allele definitions for each site an unpolarized .hap file that you wish to polarize. The ancestral/derived allele definitions are printed in exactly the same order as the sites that appear in the unpolarized .hap file.

### polarize02_polarize_haplotypes.py

This _python_ script is used to perform the actual polarization of the .hap file. It takes an unpolarized .hap file plus the output of polarize01_print_alleles_hapOrder.R as its input files, and outputs a polarized .hap file that can be used as input in Relate when constructing an Ancestral Recombination Graph (ARG).

### recombinationMap_reformatFastEPRR_forRelate.R

This _R_ script is used to reformat a recombination map generated with FastEPRR into the format required for running Relate.
