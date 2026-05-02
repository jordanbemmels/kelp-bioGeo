This repository contains custom scripts and files for the manuscript:

Bemmels et al. 2026. Spatial inference of ancestor locations suggests northern refugia for canopy-forming kelps in the Pacific Northwest. <i>Manuscript in review</i>.

There was no substantial new program or software developed for this manuscript; rather, the focus of this repository is providing a few helper scripts showing how to process and format data, as well illustrating how some of the more invovled pipelines were run. These scripts are intended to supplement to the verbal explanations in the manuscript.

See the readme files in the subdirectories (contents also pasted below) for further info about the subdirectory contents.

The pages and documents in this repository were developed by Jordan Bemmels (jbemmels@umich.edu) except for the Spacetrees files, which were developed by Matt Osmond (University of Toronto).

# ARG

These are helper scripts to assist with preparing necessary requirements for constructing Ancestral Recombination Graphs (ARGs) in Relate.

### polarize01_print_alleles_hapOrder.R

This _R_ script is an initial preparatory script used to print the ancestral and derived allele definitions for each site an unpolarized .hap file that you wish to polarize. The ancestral/derived allele definitions are printed in exactly the same order as the sites that appear in the unpolarized .hap file.

### polarize02_polarize_haplotypes.py

This _python_ script is used to perform the actual polarization of the .hap file. It takes an unpolarized .hap file plus the output of polarize01_print_alleles_hapOrder.R as its input files, and outputs a polarized .hap file that can be used as input in Relate when constructing an Ancestral Recombination Graph (ARG).

### recombinationMap_reformatFastEPRR_forRelate.R

This _R_ script is used to reformat a recombination map generated with FastEPRR into the format required for running Relate.

# ENM

These are scripts that illustrate how occurrence records were filtered and ecological niche models (ENMs) were constructed usingn MaxNet.

### filter_GBIF_occurrences.R

This _R_ script illustrates how the raw occurrence records from GBIF (gbif.org) were filtered prior to constructing ecological niche models (ENMs). The code subsets all occurrences to an appropriate extent and background area, removes occurrences flagged with severe issues, removes occurrences with low precision or extreme rounding of latitude/longitude, removes individual occurrences that are suspicious or outside the known range of either species, and removes occurrences that fall on land. It also converts the points to an Albers equal area projection focused on Vancouver Island.

### niche_models.R

This _R_ script illustrates how to run the ecological niche models in MaxNet, and produce the evaluation metrics as described in the manuscript. There is no conceptually new code or software presented here, but rather the focus is on illustrating how existing programs were run.

# Spacetrees

This directory contains Snakemake pipelines for running Spacetrees (Osmond and Coop 2004) to infer the locations of genetic ancestors. The files here do not provide conceptually novel code, but rather, illustrate how Spacetrees was run for the manuscript.

All scripts in this directory are authored by Matt Osmond (University of Toronto).

Spacetrees is available from https://github.com/osmond-lab/spacetrees

Osmond M, Coop G. 2024. Estimating dispersal rates and locating genetic ancestors with genome-wide genealogies. _eLife_ **13**: e72177.

### Snakefile-KL19

Snakemake pipeline to run Spacetrees for Macrocystis.

### Snakefile-KL20

Snakemake pipeline to run Spacetrees for Nereocystis.

### requirements.txt, startup.sh, and utils.py

Associated required files.

### spacetrees.py

Spacetrees code and functions.

# Directionality index

These are helper scripts to process ancestral/derived allele definitions and polarize SNPs prior to calculating directionality index (psi).

### get_ancestralDerived_GERP.R

This _R_ script is used to process the output from the modified gerpcol script of Taylor et al. (2024) and explicitly define the ancestral and derived allele for all SNP sites.

The modified gerpcol script that must be run prior to running this script is available from

https://github.com/BeckySTaylor/Phylogenomic_Analyses/blob/main/gerp.tar.gz

Taylor, R.S., Manseau, M., Keobouasone, S., Liu, P., Mastromonaco, G., Solmundson, K., Kelly, A., Larter, N.C., Gamberg, M., Schwantje, H., et al. (2024). High genetic load without purging in caribou, a diverse species at risk. _Curr. Biol._ **34**, 1234-1246.e7. https://doi.org/10.1016/j.cub.2024.02.002.

### vcf_to_SNAPP_polarized.R

This _R_ script is used to convert an unpolarized VCF file into a polarized SNAPP file, prior to calculating directionality index (psi).

The alleles are polarized using the ancestral-derived definitions previously generated in the get_ancestralDerived_GERP.R script.
