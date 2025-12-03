##### Code is to reformat the genetic map from FastEPRR output into the input format required by Relate. We will also fill in any gaps or windows with zero recombination rate by merging gaps/zeroes with the two non-zero recombination windows on either side, and recalculating a single average recombination rate applicable to the entire merged window (as we do not want any areas of zero recombination);

# Author use only: GitHub version derived from custom KL15 script reformatFastEPRR_forRelate_KL15_v4.R, 2025/12/03;

#####

### Prerequisites:

# User has already run the program FastEPRR. In this example, the FastEPRR output files are contained within the directory step_05_Ne10K, and there is one file per chromosome named chr1, chr2, chr3, etc. Each FastEPRR file is three-column space-delimited with the following columns: Start = starting position in base-pairs along the chromosome of the current window, End = ending position in base-pairs along the chromosome, r(cM/Mb) = recombination rate in cM/Mb;

# User has a standard .fai file resulting from indexing the fasta file for the reference genome. This can be generated with samtools as follows:

# samtools faidx Macpyr2_AssemblyScaffolds.fasta;

# User has created the "output" directory to contain the results, one file per chromosome;

#####

# get the names of all the chromosomes we have available from FastEPRR output;

allChr <- list.files("step_05_Ne10K");

# we need also the lengths of the chromsomes, from the fasta index;
# and need to convert the name to match the allChr names, as the allChr names begin with "chr" whereas the fasta file names begin with "scaffold_" or "contig_";

fai <- read.table("Macpyr2_AssemblyScaffolds.fasta.fai", sep = "\t", header = F);
fai$chrName <- gsub("scaffold_", "", fai$V1);
fai$chrName <- gsub("contig_", "", fai$chrName);
fai$chrName <- paste0("chr", fai$chrName);

#

# begin processing each chromosome;

for (i in 1:length(allChr)) {
  
  cur_chr <- allChr[i];
  cur_length <- fai$V2[fai$chrName == cur_chr];
  
  fE <- read.table(paste0("step_05_Ne10K/", cur_chr), sep = " ", header = T, stringsAsFactors = F);
  colnames(fE) <- c("Start", "End", "r");
  
  ##### fill in any gaps temporarily with a zero recombination rate - note that the first and last windows of the chromosome do not align with the exact chromosome ends, so those are always gaps (unless it literally starts at 1 or ends at exactly the length of the chromosome) and for these special positions, rather than giving a new window with zero, I am extending the first and last windows to cover the chromosome ends;
  
  if (fE$Start[1] == 1) {
    fE_mod <- fE[1, ];
  } else {
    # here we are re-calculating r for the first window to include the initial part of the chromosome;
    fE_mod <- data.frame(Start = 1, End = fE$End[1], r = fE$r[1] * (fE$End[1] - fE$Start[1])/fE$End[1]);
  }
  
  for (j in 2:(nrow(fE))) {
    
    # if the next window begins where the previous left off, simply add the next line - if not, then add an intermediate window (here temporarily with a zero recombination rate, which we'll fix later);
    if (fE$Start[j] == fE_mod$End[nrow(fE_mod)] + 1) {
      fE_mod <- rbind(fE_mod, fE[j, ]);
    } else {
      fE_mod <- rbind(fE_mod, c(fE_mod$End[nrow(fE_mod)] + 1, fE$Start[j] - 1, 0));
    }
    
  }
  
  # adjust the last line by re-calculating r for the last window to include the last part of the chromosome, if necessary;
  if (fE_mod$End[nrow(fE_mod)] != cur_length){
    fE_mod$r[nrow(fE_mod)] <- fE_mod$r[nrow(fE_mod)] * (fE_mod$End[nrow(fE_mod)] - fE_mod$Start[nrow(fE_mod)])/(cur_length - fE_mod$Start[nrow(fE_mod)]);
    fE_mod$End[nrow(fE_mod)] <- cur_length;
  }
  
  ##### if there are any windows with zero recombination rate, make the rate of that region plus its two surrounding windows equal to the average rate across all three windows (accounting for the total size of the three windows);
  ##### note that if the gap includes *multiple* windows, we have to smooth over *more than 3* windows and figure out how many to smooth it over;
  for (j in 1:nrow(fE_mod)) {
    if (fE_mod$r[j] == 0) {
      start_line <- j - 1;
      possible_end_lines <- which(fE_mod$r != 0);
      end_line <- min(possible_end_lines[possible_end_lines > j], nrow(fE_mod)); # the possibility of the end_line being nrow(fE_mod) is to account for the situation when the final line ends in zero - in that case we simply end on the final line;
      #
      avg_rate <- mean(fE_mod$r[start_line:end_line]);
      fE_mod$r[start_line:end_line] <- avg_rate;
    }
  }
  
  ##### the windows now cover the entire chromosome, do not have any gaps, and have been adjusted (if necessary) so that there are no zeroes for recombination rate;
  ##### ready to convert to relate format;
  
  # the first position begins at zero, has a recombination rate equal to the first window
  genmap <- data.frame("pos" = 0, "COMBINED_rate" = fE_mod$r[1], "Genetic_Map" = 0);
  
  # continue with the rest of the positions;
  for (j in 2:nrow(fE_mod)) {
    new_p <- fE_mod$Start[j];
    new_r <- fE_mod$r[j];
    new_rdist <- genmap$COMBINED_rate[j-1]/1000000*(new_p-genmap$pos[j-1])+genmap$Genetic_Map[j-1];
    genmap <- rbind(genmap, c(new_p, new_r, new_rdist));
  }
  
  # write;
  write.table(genmap, file = paste0("output/", cur_chr, ".txt"), sep = " ", row.names = F, col.names = T, quote = F);

}

##############################
##############################

##### Also summarize the genetic map size of each chromosome, and overall;

gmsum <- data.frame(chr = allChr, size_bp = NA, size_cM = NA);

for (i in 1:nrow(gmsum)) {
  
  cur_chr <- gmsum$chr[i];
  
  gmsum$size_bp[i] <- fai$V2[fai$chrName == cur_chr];
  
  cur_genmap <- read.table(paste0("output/", cur_chr, ".txt"), sep = " ", header = T, stringsAsFactors = F);
  
  gmsum$size_cM[i] <- cur_genmap$Genetic_Map[nrow(cur_genmap)];
  
}

gmsum <- gmsum[order(gmsum$size_bp, decreasing = T), ];

gmsum <- rbind(gmsum, c("totalAuto", sum(gmsum$size_bp), sum(gmsum$size_cM)));

gmsum$r_cMperMb <- as.numeric(gmsum$size_cM) / as.numeric(gmsum$size_bp) * 1000000;

write.table(gmsum, "genetic_map_Ne10K_summary.txt", sep = "\t", row.names = F, col.names = T, quote = F);
