# Original script by Matt Osmond, edited by Jordan Bemmels

##### This script is used to polarize the alleles in a .hap file, which will be used downstream as input for Relate (to construct an ARG)

# Author use only: GitHub version derived from custom KL19 script 02_polarize_haplotypes_KL19_v5.py, 2025/12/03

#####

### Prerequisites:

# User has already run the script polarize01_print_alleles_hapOrder.R. In this example, the output is as follows:

# alleles_withAncestralDefinitions_inHapOrder_KL19_v5.txt

# User already has an unpolarized .hap file that the user wishes to polarize. In this example, it is called:

# KL19_forARG_v5_nonMiss90_variantCheck.hap

#####

nflips = 0 #counting how many times the reference allele is derived
ntot = 0 #total number of alleles processed
nwritten = 0 #total number of alleles written (excludes sites where ancestral allele is not defined)
with open('alleles_withAncestralDefinitions_inHapOrder_KL19_v5.txt','r') as f: #alleles file with polarization 
	next(f) # skip the first header line
	with open('KL19_forARG_v5_nonMiss90_variantCheck.hap','r') as h: #unpolarized hap file
		with open('KL19_forARG_v5_nonMiss90_variantCheck_polarized.hap','w') as out: #polarized hap output file
			for i, (allele, hap) in enumerate(zip(f, h)): #one line at a time
				# get the refIsAnc column, which tells us if a polarization flip is needed
				refIsAnc = allele.split()[8] # do not specify a split, as it is the last column, if specify '\t' it doesn't work as the last character '\n' is appended then
				#
				# get ref and alt alleles
				haps = hap.split(' ')
				reference = haps[3] #reference allele
				alternate = haps[4] #alternate allele
				# only need to flip allele if reference is not ancestral
				if refIsAnc == '0':
					nflips += 1
					haps[3] = alternate #make the alternate allele ancestral
					haps[4] = reference #make the reference allele derived
					haps[5:] = [1 - int(i) for i in haps[5:]] #flip the haps (make 0 -> 1 and 1 -> 0)
					hap = (' ').join(map(str,haps)) + '\n'
				# only print if the ancestral allele is defined, otherwise do not print and the SNP site is totally removed from output
				if refIsAnc != '-9':
					# write to outfile
					out.write(hap)
					nwritten += 1
				ntot += 1
print('retained',nwritten,'out of',ntot,'alleles, as a fraction:',nwritten/ntot)
print('flipped',nflips,'out of',nwritten,'retained alleles, as a fraction:',nflips/nwritten)
