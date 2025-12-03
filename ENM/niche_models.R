##### Code illustrates how ecological niche models were constructed using MaxNet;

# Author use only: GitHub version derived from custom ENMs script kelp_niche_models_241113.R, 2025/12/03;

#####

### Prerequisites:

# Filtered environmental rasters from MARSPEC (https://marspec.weebly.com/) for current and LGM (ensemble, and CCSM3) are available that have been cropped to the NE Pacific and converted to an Albers Equal Area Projection focused on Vancouver Island (suffix "_vi_aea"), i.e., CRS = "+proj=aea +lat_0=40 +lon_0=-125 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";

# Filtered occurrence records from GBIF (gbif.org) area available, previously filtered with the script filter_GBIF_occurrences.R;

######

#require(maxnet); # we will be using MaxNet - but see note for SDMtune;
require(terra);
require(maps);
require(ENMeval); # this package does not work at all for evaluating MaxNet models;
require(SDMtune); # this package works with MaxNet models, though note that the maxnet model must be run *within* SDMtools;
require(zeallot); # required for some parts of SDMtune and not loaded automatically;

##### load current environmental data as a raster stack;
current_files <- paste0("current_5m_vi_aea/", list.files("current_5m_vi_aea"));
current <- rast(current_files);

##### load LGM environmental data as a raster stack;
# note that there are both generic geophysical environmental data plus those specific to a particular global circulation model (ensemble and CCSM3, repsectively);
LGM_ensemble_files <- c(paste0("LGM_5m_geophysical_vi_aea/", list.files("LGM_5m_geophysical_vi_aea")), paste0("LGM_5m_ensemble_noCCSM_vi_aea/", list.files("LGM_5m_ensemble_noCCSM_vi_aea")));
LGM_ensemble <- rast(LGM_ensemble_files);

LGM_ccsm_files <- c(paste0("LGM_5m_geophysical_vi_aea/", list.files("LGM_5m_geophysical_vi_aea")), paste0("LGM_5m_CCSM_vi_aea/", list.files("LGM_5m_CCSM_vi_aea")));
LGM_ccsm <- rast(LGM_ccsm_files);

##### load occurrences and background points;
macro <- read.table("GBIF_filtered/Macrocystis_global_241010_filtered.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "", comment = "");
nereo <- read.table("GBIF_filtered/Nereocystis_global_241010_filtered.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "", comment = "");
phaeo <- read.table("GBIF_filtered/Phaeophyceae_global_241010_filtered.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "", comment = "");

########## EXTRACT ENVIRONMENTAL DATA FOR RELEVANT CELLS ##########

##### add cell number to the occurrences and background points;
macro$cell <- cellFromXY(current[[1]], macro[ , c("Long_vi_aea", "Lat_vi_aea")]);
nereo$cell <- cellFromXY(current[[1]], nereo[ , c("Long_vi_aea", "Lat_vi_aea")]);
phaeo$cell <- cellFromXY(current[[1]], phaeo[ , c("Long_vi_aea", "Lat_vi_aea")]);

n_occs_all_macro <- length(unique(macro$cell));
n_occs_all_nereo <- length(unique(nereo$cell));
n_occs_all_phaeo <- length(unique(phaeo$cell));

n_occs_all_macro; # 227;
n_occs_all_nereo; # 198;
n_occs_all_phaeo; # 1648;

##### add distance from shore;
macro$biogeo05_5m <- unlist(extract(current[["biogeo05_5m"]], macro$cell));
nereo$biogeo05_5m <- unlist(extract(current[["biogeo05_5m"]], nereo$cell));
phaeo$biogeo05_5m <- unlist(extract(current[["biogeo05_5m"]], phaeo$cell));

########## CREATE A VERSION OF RASTERS THAT ONLY COVERS COASTLINE ##########

##### for some purposes, such as assessing correlations among variables, we want to consider only the pixels adjacent to the coastline - pixels out in the open ocean are irrelevant to us and should not be used;
# note that biogeo05_5m is distance to shore, in km;

# check how far away from the coastline our species are actually found;

hist(macro$biogeo05_5m); # there is one outlier at remote rocks in the middle of the ocean, we can ignore this one;
hist(nereo$biogeo05_5m);
hist(phaeo$biogeo05_5m);

max(macro$biogeo05_5m[macro$biogeo05_5m < 100]); # 11.6 km maximum;
max(nereo$biogeo05_5m); # 10.8 km maximum;

# how would a 12-km buffer look;
plot(current$biogeo05_5m <= 12);

# what percent of unique cell records are retained;

length(unique(macro$cell[macro$biogeo05_5m <= 12])) / length(unique(macro$cell)); #99.6%;
length(unique(nereo$cell[nereo$biogeo05_5m <= 12])) / length(unique(nereo$cell)); #100%;
length(unique(phaeo$cell[phaeo$biogeo05_5m <= 12])) / length(unique(phaeo$cell)); #93.4%;

# create a filtered version of data with a 12-km buffer;

current_buffer12km <- current;
LGM_ensemble_buffer12km <- LGM_ensemble;
LGM_ccsm_buffer12km <- LGM_ccsm;

current_buffer12km[current_buffer12km$biogeo05_5m > 12] <- NA;
LGM_ensemble_buffer12km[LGM_ensemble_buffer12km$biogeo05_5m > 12] <- NA;
LGM_ccsm_buffer12km[LGM_ccsm_buffer12km$biogeo05_5m > 12] <- NA;

########## FILTER OCCURRENCES TO ONLY THOSE ALONG COASTLINE ##########

macro <- macro[macro$biogeo05_5m <= 12, ];
nereo <- nereo[nereo$biogeo05_5m <= 12, ];
phaeo <- phaeo[phaeo$biogeo05_5m <= 12, ];

##### retain only required variables and extract variables;

macro_occs <- macro[ , c("Long_vi_aea", "Lat_vi_aea", "cell")];
for (i in 1:nlyr(current)) {
	macro_occs[ , names(current[[i]])] <- extract(current[[i]], macro_occs$cell);
}

nereo_occs <- nereo[ , c("Long_vi_aea", "Lat_vi_aea", "cell")];
for (i in 1:nlyr(current)) {
	nereo_occs[ , names(current[[i]])] <- extract(current[[i]], nereo_occs$cell);
}

phaeo_occs <- phaeo[ , c("Long_vi_aea", "Lat_vi_aea", "cell")];
for (i in 1:nlyr(current)) {
	phaeo_occs[ , names(current[[i]])] <- extract(current[[i]], phaeo_occs$cell);
}

##### final number of unique occurrences;

n_occs_all_macro <- length(unique(macro_occs$cell));
n_occs_all_nereo <- length(unique(nereo_occs$cell));
n_occs_all_phaeo <- length(unique(phaeo_occs$cell));

n_occs_all_macro; # 226;
n_occs_all_nereo; # 198;
n_occs_all_phaeo; # 1540;

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

########## MACROCYSTIS ##########

##### Pre-process data #####

# prepare 20 different datasets where 80% of unique occurrences and 80% of unique backgrounds are retained;
# these can be used for having 20 different replicates of the model, and averaging the results, to generate some uncertainty and to ensure that the results do not depend on the exact set of presence / background records;
# note that we can NOT simply use the randomFolds() function from SDMtune, because this function would divide the entire set of records *including duplicates*, which we do not want to include (more than one observation per raster cell);
# however, we cannot simply remove duplicates first and then use randomFolds() to address that issue, because then we lose information about sampling effot in certain pixels that are extremely overrepresented, distoring the sampling effort;
# thus, I need to customize the datasets myself;

# logic is to randomly resort the dataframe then remove duplicates for cell number, then take the first 100 samples;

macro_pool <- macro_occs;
macro_pool$ID <- 1:nrow(macro_occs);

macro_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(macro_pool), ncol = 20));

for (i in 1:20) {
	
	macro_randomized <- macro_pool[sample(nrow(macro_pool)), ];
	macro_randomized_noDuplicates <- macro_randomized[!(duplicated(macro_randomized$cell)), ];
	macro_randomized_IDs <- macro_randomized_noDuplicates$ID[1:round(nrow(macro_randomized_noDuplicates) * 0.8)];
	
	rownames(macro_randomized_noDuplicates[1:round(nrow(macro_randomized_noDuplicates) * 0.8), ]);
	
	macro_sampling_20reps[macro_pool$ID %in% macro_randomized_IDs, i] <- T;

}
 
#

phaeo_pool <- phaeo_occs;
phaeo_pool$ID <- 1:nrow(phaeo_occs);

phaeo_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(phaeo_pool), ncol = 20));

for (i in 1:20) {
	
	phaeo_randomized <- phaeo_pool[sample(nrow(phaeo_pool)), ];
	phaeo_randomized_noDuplicates <- phaeo_randomized[!(duplicated(phaeo_randomized$cell)), ];
	phaeo_randomized_IDs <- phaeo_randomized_noDuplicates$ID[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8)];
	
	rownames(phaeo_randomized_noDuplicates[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8), ]);
	
	phaeo_sampling_20reps[phaeo_pool$ID %in% phaeo_randomized_IDs, i] <- T;

}

#####
 
##### Select top correlated variables to retain #####

# FOR THIS VERSION HERE, USE CORRELATED VARIABLE RETENTION BASED ON BIOLOGICAL PRINCIPLES;

#

# start with all variables;

retained_var <- names(current_buffer12km);

# inspect correlations - only within our 12km buffer where likely to be relevant;

bg_sample <- spatSample(current_buffer12km, size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg <- prepareSWD(species = "Background points", a = bg_sample, env = current_buffer12km);

plotCor(bg, method = "spearman", cor_th = 0.8);

# retain biogeo08_5m (mean annual SSS) - remove all variables that are correlated;

env_cor <- corVar(bg, method = "spearman", cor_th = 0.8);

env_cor$Var1 <- as.character(env_cor$Var1); # don't want factor;
env_cor$Var2 <- as.character(env_cor$Var2); # don't want factor;

to_remove <- c(env_cor$Var2[env_cor$Var1 == "biogeo08_5m"], env_cor$Var1[env_cor$Var2 == "biogeo08_5m"]);

retained_var <- retained_var[!(retained_var %in% to_remove)];

# re-examine correlations;

bg_sample_retained <- spatSample(current_buffer12km[[retained_var]], size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg_retained <- prepareSWD(species = "Background points", a = bg_sample_retained, env = current_buffer12km[[retained_var]]);

plotCor(bg_retained, method = "spearman", cor_th = 0.8);

# next, retain biogeo15_5m (SST of the warmest ice-free month) - remove all variables that are correlated;

to_remove <- c(env_cor$Var2[env_cor$Var1 == "biogeo15_5m"], env_cor$Var1[env_cor$Var2 == "biogeo15_5m"]);

retained_var <- retained_var[!(retained_var %in% to_remove)];

# re-examine correlations;

bg_sample_retained <- spatSample(current_buffer12km[[retained_var]], size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg_retained <- prepareSWD(species = "Background points", a = bg_sample_retained, env = current_buffer12km[[retained_var]]);

plotCor(bg_retained, method = "spearman", cor_th = 0.8);

# we have biogeo16_5m and biogeo17_5m at 0.99 correlation - it likely does not matter which we remove with this level of correlation - arbitrarily retain biogeo16_5m (range in SST) rather than biogeo17_5m (variance in SST) as range is a bit more intuitive;

retained_var <- retained_var[!(retained_var == "biogeo17_5m")];

# we also have bathy_5m and biogeo06_5m correlated at -0.85 - bathy_5m is bathymetry whereas biogeo06_5m is bathymetric slope - the raw bathymetry is a priori much more biologically relevant to kelp distribution than the bathymetric slope (especially at large pixel size where slope over an enormous area is likely irrelevant) - remove bathymetric slope and retain raw bathymetry;

retained_var <- retained_var[!(retained_var == "biogeo06_5m")];

# confirm that no retained variables are correlated above the 0.8 threshold;

bg_sample_retained <- spatSample(current_buffer12km[[retained_var]], size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg_retained <- prepareSWD(species = "Background points", a = bg_sample_retained, env = current_buffer12km[[retained_var]]);

plotCor(bg_retained, method = "spearman", cor_th = 0.8);

# final list of retained uncorrelated variables from this step;

uncor_var_macro <- retained_var;

# save for later to not have to repeat;
#write.table(as.data.frame(uncor_var_macro), file = "Model_output_241113/Macrocystis/env_var_uncorrelated.txt", sep = "\t", row.names = F, col.names = F, quote = F);

uncor_var_macro <- unlist(as.vector(read.table("Model_output_241113/Macrocystis/env_var_uncorrelated.txt", header = F, sep = "\t", stringsAsFactors = F)));

#####

##### Loop through the 20 sampling replicates to select top correlated variables to retain #####

# note that we will use the same logic as reduceVar(), except that we must use a custom function because we want to average across 20 sampling replicates, instead of using only the full dataset;

minImp = 0;
retained_var_macro <- uncor_var_macro;

while(minImp < 5) {
	
	for (i in 1:20) {
				
		dat <- prepareSWD(species = "Macrocystis", env = current_buffer12km[[retained_var_macro]], p = macro_occs[macro_sampling_20reps[ , i], 1:2], a = phaeo_occs[phaeo_sampling_20reps[ , i], 1:2]);
		
		# again, we need to find which fc doesn't result in errors;
		#model <- train(method = "Maxnet", data = dat, fc = "lqph");
		model <- train(method = "Maxnet", data = dat, fc = "lqp");
		
		res_imp <- varImp(model, permut = 20);
		
		if (i == 1) {
			imp_20reps <- res_imp[ , c("Variable", "Permutation_importance")];
		} else {
			colnames(res_imp)[2] <- paste0("Permutation_importance_", i);
			imp_20reps <- merge(imp_20reps, res_imp[ , 1:2]);
		}
		
	}
	
	imp_20reps$mean <- apply(imp_20reps[ , 2:21], 1, mean);
	
	# sort so that the first row is the least important variable;
	imp_20reps <- imp_20reps[order(imp_20reps$mean), ];

	minImp <- imp_20reps$mean[1];
			
	if (minImp < 5) {
		retained_var_macro <- retained_var_macro[!(retained_var_macro == imp_20reps$Variable[1])];
		cat(paste0("Retained variables: ", retained_var_macro, collapse = ""));
	}	
	
}

# final list of (previously filtered to be uncorrelated and) important variables;

final_var_macro <- retained_var_macro;

# save for later to not have to repeat;
#write.table(as.data.frame(final_var_macro), file = "Model_output_241113/Macrocystis/env_var_final.txt", sep = "\t", row.names = F, col.names = F, quote = F);

final_var_macro <- unlist(as.vector(read.table("Model_output_241113/Macrocystis/env_var_final.txt", header = F, sep = "\t", stringsAsFactors = F)));

#####

##### Tune model hyperparameters on the full dataset #####

# We are using the full dataset instead of 20 replicates because it is conceptually difficult to combine the results of 20 different hyperparameter tunings - instead, just use the full data;

# Use cross-validation rather than testing and training data;

# reload variables if necessary - continuing from where left off previously;

final_var_macro <- unlist(as.vector(read.table("Model_output_241113/Macrocystis/env_var_final.txt", header = F, sep = "\t", stringsAsFactors = F)));

#

dat <- prepareSWD(species = "Macrocystis", env = current_buffer12km[[final_var_macro]], p = macro_occs[!(duplicated(macro_occs$cell)), 1:2], a = phaeo_occs[!(duplicated(phaeo_occs$cell)), 1:2]);

folds <- randomFolds(dat, k = 5, only_presence = F); # use only_presence = F to also permute the background observations;

cv_model <- train("Maxnet", data = dat, folds = folds, fc = "lqpht");

h_grid <- list(reg = seq(0.2, 3, 0.2), fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"));

exp_grid <- gridSearch(cv_model, hypers = h_grid, metric = "auc"); # do not specify test data, as it is contained within the folds of cv_model_test;

plot(exp_grid, title = "Grid search results");

exp_grid@results;

# get the best model - the one that has the highest AUC in testing data;

head(exp_grid@results[order(-exp_grid@results$test_AUC), ]);

index <- which.max(exp_grid@results$test_AUC);

best_hyperparam <- exp_grid@results[index, ];

# save for later to not have to repeat;
#write.table(best_hyperparam, file = "Model_output_241113/Macrocystis/hyperparam_final.txt", sep = "\t", row.names = F, col.names = T, quote = F);

best_hyperparam <- read.table("Model_output_241113/Macrocystis/hyperparam_final.txt", sep = "\t", header = T, stringsAsFactors = F);

#########################
#########################
#########################

##### Ready to run the final models, evaluate, and make predictions #####

# We want to get predictions for the 20 resampled datasets, so that the final results will incorporate uncertainty and not rely too strongly on the exact sets of presences and absences;

# reload variables and re-create 20 sampling reps, as needed;

final_var_macro <- unlist(as.vector(read.table("Model_output_241113/Macrocystis/env_var_final.txt", header = F, sep = "\t", stringsAsFactors = F)));

best_hyperparam <- read.table("Model_output_241113/Macrocystis/hyperparam_final.txt", sep = "\t", header = T, stringsAsFactors = F);

macro_pool <- macro_occs;
macro_pool$ID <- 1:nrow(macro_occs);

macro_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(macro_pool), ncol = 20));

for (i in 1:20) {
	
	macro_randomized <- macro_pool[sample(nrow(macro_pool)), ];
	macro_randomized_noDuplicates <- macro_randomized[!(duplicated(macro_randomized$cell)), ];
	macro_randomized_IDs <- macro_randomized_noDuplicates$ID[1:round(nrow(macro_randomized_noDuplicates) * 0.8)];
	
	rownames(macro_randomized_noDuplicates[1:round(nrow(macro_randomized_noDuplicates) * 0.8), ]);
	
	macro_sampling_20reps[macro_pool$ID %in% macro_randomized_IDs, i] <- T;

}
 
#

phaeo_pool <- phaeo_occs;
phaeo_pool$ID <- 1:nrow(phaeo_occs);

phaeo_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(phaeo_pool), ncol = 20));

for (i in 1:20) {
	
	phaeo_randomized <- phaeo_pool[sample(nrow(phaeo_pool)), ];
	phaeo_randomized_noDuplicates <- phaeo_randomized[!(duplicated(phaeo_randomized$cell)), ];
	phaeo_randomized_IDs <- phaeo_randomized_noDuplicates$ID[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8)];
	
	rownames(phaeo_randomized_noDuplicates[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8), ]);
	
	phaeo_sampling_20reps[phaeo_pool$ID %in% phaeo_randomized_IDs, i] <- T;

}

# get the environmental datasets;

current_buffer12km_macro <- current_buffer12km[[final_var_macro]];
LGM_ensemble_buffer12km_macro <- LGM_ensemble_buffer12km[[final_var_macro]];
LGM_ccsm_buffer12km_macro <- LGM_ccsm_buffer12km[[final_var_macro]];

# set up a dataframe for the validation results;

val_macro <- data.frame(rep = 1:20, AUC_training = NA, AUC_testing = NA, TSS_training = NA, TSS_testing = NA, AICc_training = NA, AICc_testing = NA);

# ready to run;

for (i in 1:20) {
		
	##### Prepare the data for each replicate;
	
	dat <- prepareSWD(species = "Macrocystis", env = current_buffer12km_macro, p = macro_occs[macro_sampling_20reps[ , i], 1:2], a = phaeo_occs[phaeo_sampling_20reps[ , i], 1:2]);

	##### Split the data into folds, to use with cross-validation;
	
	folds <- randomFolds(dat, k = 5, only_presence = F); # use only_presence = F to also permute the background observations;

	##### Train the model (i.e., run the model) with maxnet;

	model <- train(method = "Maxnet", data = dat, folds = folds, reg = best_hyperparam$reg, fc = best_hyperparam$fc);

	##### get the projections on a map;

	map_current_buffer12km <- predict(model, data = current_buffer12km_macro, type = "cloglog");
	map_LGM_ensemble_buffer12km <- predict(model, data = LGM_ensemble_buffer12km_macro, type = "cloglog");
	map_LGM_ccsm_buffer12km <- predict(model, data = LGM_ccsm_buffer12km_macro, type = "cloglog");

	##### save the projections;
	
	writeRaster(map_current_buffer12km, paste0("Model_output_241113/Macrocystis/projections/current_", i, ".tif"), overwrite = T);
	writeRaster(map_LGM_ensemble_buffer12km, paste0("Model_output_241113/Macrocystis/projections/LGM_ensemble_", i, ".tif"), overwrite = T);
	writeRaster(map_LGM_ccsm_buffer12km, paste0("Model_output_241113/Macrocystis/projections/LGM_ccsm_", i, ".tif"), overwrite = T);

	##### get the model evaluation metrics on the training and test data - note that AICc doesn't run on cross-validation models (and doesn't seem relevant in this situation);
	
	val_macro$AUC_training[i] <- auc(model);
	val_macro$AUC_testing[i] <- auc(model, test = T);
	
	val_macro$TSS_training[i] <- tss(model);
	val_macro$TSS_testing[i] <- tss(model, test = T);
	
	##### evaluate variable importance and variable response curves;
	
	if (i == 1) {
		pi_macro <- varImp(model, permut = 20)[ , 1:2];
		colnames(pi_macro) <- c("Variable", "PI_1");
		#
		jk_macro <- doJk(model, metric = "auc"); # no test data are specified with a cross-validation model;
		colnames(jk_macro) <- c("Variable", "Train_AUC_without_1", "Train_AUC_withonly_1");
	} else {
		new_pi_macro <- varImp(model, permut = 20)[ , 1:2];
		colnames(new_pi_macro) <- c("Variable", paste0("PI_", i));
		pi_macro <- merge(pi_macro, new_pi_macro);
		#
		new_jk_macro <- doJk(model, metric = "auc"); # no test data are specified with a cross-validation model;
		colnames(new_jk_macro) <- c("Variable", paste0("Train_AUC_without_", i), paste0("Train_AUC_withonly_", i));
		jk_macro <- merge(jk_macro, new_jk_macro);
	}
	
	for (j in 1:length(final_var_macro)) {
		outname <- paste0("Model_output_241113/Macrocystis/response_curves/current_", final_var_macro[j], "_", i, ".pdf")
		pdf(outname);
		print(plotResponse(model, var = final_var_macro[j], type = "cloglog", only_presence = T, marginal = F, rug = T));
		dev.off();
	}

	##### get common thresholds for presence-absence;
	
	# note that the thresholds() function requires a single model, which does not work with a SDMmodelCV object - instead, we can get the five individual models of SDMmodel objects and calculate thresholds individually on each of them, then get the mean (the mean will correspond to and be compatible with the mean predicted habitat suitability from the SDMmodelCV, which is simply the mean of the five model predictions);
	
	ths_1 <- thresholds(model@models[[1]], type = "cloglog");
	ths_2 <- thresholds(model@models[[2]], type = "cloglog");
	ths_3 <- thresholds(model@models[[3]], type = "cloglog");
	ths_4 <- thresholds(model@models[[4]], type = "cloglog");
	ths_5 <- thresholds(model@models[[5]], type = "cloglog");
	
	ths_combined <- cbind(ths_1[1:2], ths_2[2], ths_3[2], ths_4[2], ths_5[2]);
				
	if (i == 1) {
		thresholds_macro <- data.frame(rep = 1, MTP = mean(as.numeric(ths_combined[1, 2:6])), ETSS = mean(as.numeric(ths_combined[2, 2:6])), MTSS = mean(as.numeric(ths_combined[3, 2:6])));
	} else {
		thresholds_macro <- rbind(thresholds_macro, c(i, mean(as.numeric(ths_combined[1, 2:6])), mean(as.numeric(ths_combined[2, 2:6])), mean(as.numeric(ths_combined[3, 2:6]))));
	}

}

write.table(thresholds_macro, "Model_output_241113/Macrocystis/results_presenceAbsence_thresholds.txt", sep = "\t", row.names = F, col.names = T, quote = F);

write.table(val_macro, "Model_output_241113/Macrocystis/results_modelMetrics.txt", sep = "\t", row.names = F, col.names = T, quote = F);

#

pi_macro$PI_mean <- apply(pi_macro[ , 2:21], 1, mean);
pi_macro$PI_sd <- apply(pi_macro[ , 2:21], 1, sd);

jk_macro$Train_AUC_without_mean <- apply(jk_macro[ , paste0("Train_AUC_without_", 1:20)], 1, mean);
jk_macro$Train_AUC_without_sd <- apply(jk_macro[ , paste0("Train_AUC_without_", 1:20)], 1, sd);

jk_macro$Train_AUC_withonly_mean <- apply(jk_macro[ , paste0("Train_AUC_withonly_", 1:20)], 1, mean);
jk_macro$Train_AUC_withonly_sd <- apply(jk_macro[ , paste0("Train_AUC_withonly_", 1:20)], 1, sd);

write.table(pi_macro, "Model_output_241113/Macrocystis/results_permutation_importance.txt", sep = "\t", row.names = F, col.names = T, quote = F);
write.table(jk_macro, "Model_output_241113/Macrocystis/results_jackknife.txt", sep = "\t", row.names = F, col.names = T, quote = F);

# get a single final raster, plus uncertainty, and make a map;

macro_stack_current_files <- paste0("Model_output_241113/Macrocystis/projections/", list.files("Model_output_241113/Macrocystis/projections/", pattern = "current"));
macro_stack_current <- rast(macro_stack_current_files);

macro_current_mean_projection <- mean(macro_stack_current);

plot(macro_current_mean_projection);
plot(macro_current_mean_projection, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_current_mean_projection, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

# 

macro_stack_LGM_ensemble_files <- paste0("Model_output_241113/Macrocystis/projections/", list.files("Model_output_241113/Macrocystis/projections/", pattern = "LGM_ensemble"));
macro_stack_LGM_ensemble <- rast(macro_stack_LGM_ensemble_files);

macro_LGM_ensemble_mean_projection <- mean(macro_stack_LGM_ensemble);

plot(macro_LGM_ensemble_mean_projection);
plot(macro_LGM_ensemble_mean_projection, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_LGM_ensemble_mean_projection, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

macro_stack_LGM_ccsm_files <- paste0("Model_output_241113/Macrocystis/projections/", list.files("Model_output_241113/Macrocystis/projections/", pattern = "LGM_ccsm"));
macro_stack_LGM_ccsm <- rast(macro_stack_LGM_ccsm_files);

macro_LGM_ccsm_mean_projection <- mean(macro_stack_LGM_ccsm);

plot(macro_LGM_ccsm_mean_projection);
plot(macro_LGM_ccsm_mean_projection, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_LGM_ccsm_mean_projection, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

# plot presence/absence using:
	# MTP - minimum training presence;
	# ETSS - equal training sensitivity and specificity;
	# MTSS - maximum training sensitivity plus specificity;
	# Q05 - 5 percent quantile (my own method);

thresholds_macro <- read.table("Model_output_241113/Macrocystis/results_presenceAbsence_thresholds.txt", sep = "\t", header = T, stringsAsFactors = F);

#

macro_stack_current_mtp <- macro_stack_current;
for (i in 1:20) {
	macro_stack_current_mtp[[i]] <- macro_stack_current_mtp[[i]] >= thresholds_macro$MTP[i];
}
macro_current_mtp_sum <- sum(macro_stack_current_mtp);

plot(macro_current_mtp_sum);
plot(macro_current_mtp_sum, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_current_mtp_sum, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

plot(macro_current_mean_projection > mean(thresholds_macro$MTP)); # alternative method;
plot(macro_current_mean_projection > mean(thresholds_macro$MTP), xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_current_mean_projection > mean(thresholds_macro$MTP), xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

macro_stack_current_mtss <- macro_stack_current;
for (i in 1:20) {
	macro_stack_current_mtss[[i]] <- macro_stack_current_mtss[[i]] >= thresholds_macro$MTSS[i];
}
macro_current_mtss_sum <- sum(macro_stack_current_mtss);

plot(macro_current_mtss_sum);
plot(macro_current_mtss_sum, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_current_mtss_sum, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

plot(macro_current_mean_projection > mean(thresholds_macro$MTSS)); # alternative method;
plot(macro_current_mean_projection > mean(thresholds_macro$MTSS), xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_current_mean_projection > mean(thresholds_macro$MTSS), xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

macro_stack_current_q05 <- macro_stack_current;
for (i in 1:20) {
	macro_stack_current_q05[[i]] <- macro_stack_current_q05[[i]] >= quantile(unlist(extract(macro_stack_current[[i]], unique(macro_occs$cell))), 0.05);
}
macro_current_q05_sum <- sum(macro_stack_current_q05);

plot(macro_current_q05_sum);
plot(macro_current_q05_sum, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(macro_current_q05_sum, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

############################################################
############################################################
############################################################

########## NEREOCYSTIS ##########

##### Pre-process data #####

# prepare 20 different datasets where 80% of unique occurrences and 80% of unique backgrounds are retained;
# these can be used for having 20 different replicates of the model, and averaging the results, to generate some uncertainty and to ensure that the results do not depend on the exact set of presence / background records;
# note that we can NOT simply use the randomFolds() function from SDMtune, because this function would divide the entire set of records *including duplicates*, which we do not want to include (more than one observation per raster cell);
# however, we cannot simply remove duplicates first and then use randomFolds() to address that issue, because then we lose information about sampling effot in certain pixels that are extremely overrepresented, distoring the sampling effort;
# thus, I need to customize the datasets myself;

# logic is to randomly resort the dataframe then remove duplicates for cell number, then take the first 100 samples;

nereo_pool <- nereo_occs;
nereo_pool$ID <- 1:nrow(nereo_occs);

nereo_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(nereo_pool), ncol = 20));

for (i in 1:20) {
	
	nereo_randomized <- nereo_pool[sample(nrow(nereo_pool)), ];
	nereo_randomized_noDuplicates <- nereo_randomized[!(duplicated(nereo_randomized$cell)), ];
	nereo_randomized_IDs <- nereo_randomized_noDuplicates$ID[1:round(nrow(nereo_randomized_noDuplicates) * 0.8)];
	
	rownames(nereo_randomized_noDuplicates[1:round(nrow(nereo_randomized_noDuplicates) * 0.8), ]);
	
	nereo_sampling_20reps[nereo_pool$ID %in% nereo_randomized_IDs, i] <- T;

}
 
#

phaeo_pool <- phaeo_occs;
phaeo_pool$ID <- 1:nrow(phaeo_occs);

phaeo_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(phaeo_pool), ncol = 20));

for (i in 1:20) {
	
	phaeo_randomized <- phaeo_pool[sample(nrow(phaeo_pool)), ];
	phaeo_randomized_noDuplicates <- phaeo_randomized[!(duplicated(phaeo_randomized$cell)), ];
	phaeo_randomized_IDs <- phaeo_randomized_noDuplicates$ID[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8)];
	
	rownames(phaeo_randomized_noDuplicates[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8), ]);
	
	phaeo_sampling_20reps[phaeo_pool$ID %in% phaeo_randomized_IDs, i] <- T;

}

#####
 
##### Select top correlated variables to retain #####

# FOR THIS VERSION HERE, USE CORRELATED VARIABLE RETENTION BASED ON BIOLOGICAL PRINCIPLES;

#

# start with all variables;

retained_var <- names(current_buffer12km);

# inspect correlations - only within our 12km buffer where likely to be relevant;

bg_sample <- spatSample(current_buffer12km, size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg <- prepareSWD(species = "Background points", a = bg_sample, env = current_buffer12km);

plotCor(bg, method = "spearman", cor_th = 0.8);

# retain biogeo08_5m (mean annual SSS) - remove all variables that are correlated;

env_cor <- corVar(bg, method = "spearman", cor_th = 0.8);

env_cor$Var1 <- as.character(env_cor$Var1); # don't want factor;
env_cor$Var2 <- as.character(env_cor$Var2); # don't want factor;

to_remove <- c(env_cor$Var2[env_cor$Var1 == "biogeo08_5m"], env_cor$Var1[env_cor$Var2 == "biogeo08_5m"]);

retained_var <- retained_var[!(retained_var %in% to_remove)];

# re-examine correlations;

bg_sample_retained <- spatSample(current_buffer12km[[retained_var]], size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg_retained <- prepareSWD(species = "Background points", a = bg_sample_retained, env = current_buffer12km[[retained_var]]);

plotCor(bg_retained, method = "spearman", cor_th = 0.8);

# next, retain biogeo15_5m (SST of the warmest ice-free month) - remove all variables that are correlated;

to_remove <- c(env_cor$Var2[env_cor$Var1 == "biogeo15_5m"], env_cor$Var1[env_cor$Var2 == "biogeo15_5m"]);

retained_var <- retained_var[!(retained_var %in% to_remove)];

# re-examine correlations;

bg_sample_retained <- spatSample(current_buffer12km[[retained_var]], size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg_retained <- prepareSWD(species = "Background points", a = bg_sample_retained, env = current_buffer12km[[retained_var]]);

plotCor(bg_retained, method = "spearman", cor_th = 0.8);

# we have biogeo16_5m and biogeo17_5m at 0.99 correlation - it likely does not matter which we remove with this level of correlation - arbitrarily retain biogeo16_5m (range in SST) rather than biogeo17_5m (variance in SST) as range is a bit more intuitive;

retained_var <- retained_var[!(retained_var == "biogeo17_5m")];

# we also have bathy_5m and biogeo06_5m correlated at -0.85 - bathy_5m is bathymetry whereas biogeo06_5m is bathymetric slope - the raw bathymetry is a priori much more biologically relevant to kelp distribution than the bathymetric slope (especially at large pixel size where slope over an enormous area is likely irrelevant) - remove bathymetric slope and retain raw bathymetry;

retained_var <- retained_var[!(retained_var == "biogeo06_5m")];

# confirm that no retained variables are correlated above the 0.8 threshold;

bg_sample_retained <- spatSample(current_buffer12km[[retained_var]], size = 10000, method = "random", na.rm = T, xy = T, values = F);
bg_retained <- prepareSWD(species = "Background points", a = bg_sample_retained, env = current_buffer12km[[retained_var]]);

plotCor(bg_retained, method = "spearman", cor_th = 0.8);

# final list of retained uncorrelated variables from this step;

uncor_var_nereo <- retained_var;

# save for later to not have to repeat;
#write.table(as.data.frame(uncor_var_nereo), file = "Model_output_241113/Nereocystis/env_var_uncorrelated.txt", sep = "\t", row.names = F, col.names = F, quote = F);

uncor_var_nereo <- unlist(as.vector(read.table("Model_output_241113/Nereocystis/env_var_uncorrelated.txt", header = F, sep = "\t", stringsAsFactors = F)));

#####

##### Loop through the 20 sampling replicates to select top correlated variables to retain #####

# note that we will use the same logic as reduceVar(), except that we must use a custom function because we want to average across 20 sampling replicates, instead of using only the full dataset;

minImp = 0;
retained_var_nereo <- uncor_var_nereo;

while(minImp < 5) {
	
	for (i in 1:20) {
				
		dat <- prepareSWD(species = "Nereocystis", env = current_buffer12km[[retained_var_nereo]], p = nereo_occs[nereo_sampling_20reps[ , i], 1:2], a = phaeo_occs[phaeo_sampling_20reps[ , i], 1:2]);
		
		# again, we need to find which fc doesn't result in errors;
		model <- train(method = "Maxnet", data = dat, fc = "lqph");
		#model <- train(method = "Maxnet", data = dat, fc = "lqp");
		
		res_imp <- varImp(model, permut = 20);
		
		if (i == 1) {
			imp_20reps <- res_imp[ , c("Variable", "Permutation_importance")];
		} else {
			colnames(res_imp)[2] <- paste0("Permutation_importance_", i);
			imp_20reps <- merge(imp_20reps, res_imp[ , 1:2]);
		}
		
	}
	
	imp_20reps$mean <- apply(imp_20reps[ , 2:21], 1, mean);
	
	# sort so that the first row is the least important variable;
	imp_20reps <- imp_20reps[order(imp_20reps$mean), ];

	minImp <- imp_20reps$mean[1];
			
	if (minImp < 5) {
		retained_var_nereo <- retained_var_nereo[!(retained_var_nereo == imp_20reps$Variable[1])];
		cat(paste0("Retained variables: ", retained_var_nereo, collapse = ""));
	}	
	
}

# final list of (previously filtered to be uncorrelated and) important variables;

final_var_nereo <- retained_var_nereo;

# save for later to not have to repeat;
#write.table(as.data.frame(final_var_nereo), file = "Model_output_241113/Nereocystis/env_var_final.txt", sep = "\t", row.names = F, col.names = F, quote = F);

final_var_nereo <- unlist(as.vector(read.table("Model_output_241113/Nereocystis/env_var_final.txt", header = F, sep = "\t", stringsAsFactors = F)));

#####

##### Tune model hyperparameters on the full dataset #####

# We are using the full dataset instead of 20 replicates because it is conceptually difficult to combine the results of 20 different hyperparameter tunings - instead, just use the full data;

# Use cross-validation rather than testing and training data;

# reload variables if necessary - continuing from where left off previously;

final_var_nereo <- unlist(as.vector(read.table("Model_output_241113/Nereocystis/env_var_final.txt", header = F, sep = "\t", stringsAsFactors = F)));

#

dat <- prepareSWD(species = "Nereocystis", env = current_buffer12km[[final_var_nereo]], p = nereo_occs[!(duplicated(nereo_occs$cell)), 1:2], a = phaeo_occs[!(duplicated(phaeo_occs$cell)), 1:2]);

folds <- randomFolds(dat, k = 5, only_presence = F); # use only_presence = F to also permute the background observations;

cv_model <- train("Maxnet", data = dat, folds = folds, fc = "lqpht");

h_grid <- list(reg = seq(0.2, 3, 0.2), fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"));

exp_grid <- gridSearch(cv_model, hypers = h_grid, metric = "auc"); # do not specify test data, as it is contained within the folds of cv_model_test;

plot(exp_grid, title = "Grid search results");

exp_grid@results;

# get the best model - the one that has the highest AUC in testing data;

head(exp_grid@results[order(-exp_grid@results$test_AUC), ]);

index <- which.max(exp_grid@results$test_AUC);

best_hyperparam <- exp_grid@results[index, ];

# save for later to not have to repeat;
#write.table(best_hyperparam, file = "Model_output_241113/Nereocystis/hyperparam_final.txt", sep = "\t", row.names = F, col.names = T, quote = F);

best_hyperparam <- read.table("Model_output_241113/Nereocystis/hyperparam_final.txt", sep = "\t", header = T, stringsAsFactors = F);

#########################
#########################
#########################

##### Ready to run the final models, evaluate, and make predictions #####

# We want to get predictions for the 20 resampled datasets, so that the final results will incorporate uncertainty and not rely too strongly on the exact sets of presences and absences;

# reload variables and re-create 20 sampling reps, as needed;

final_var_nereo <- unlist(as.vector(read.table("Model_output_241113/Nereocystis/env_var_final.txt", header = F, sep = "\t", stringsAsFactors = F)));

best_hyperparam <- read.table("Model_output_241113/Nereocystis/hyperparam_final.txt", sep = "\t", header = T, stringsAsFactors = F);

nereo_pool <- nereo_occs;
nereo_pool$ID <- 1:nrow(nereo_occs);

nereo_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(nereo_pool), ncol = 20));

for (i in 1:20) {
	
	nereo_randomized <- nereo_pool[sample(nrow(nereo_pool)), ];
	nereo_randomized_noDuplicates <- nereo_randomized[!(duplicated(nereo_randomized$cell)), ];
	nereo_randomized_IDs <- nereo_randomized_noDuplicates$ID[1:round(nrow(nereo_randomized_noDuplicates) * 0.8)];
	
	rownames(nereo_randomized_noDuplicates[1:round(nrow(nereo_randomized_noDuplicates) * 0.8), ]);
	
	nereo_sampling_20reps[nereo_pool$ID %in% nereo_randomized_IDs, i] <- T;

}
 
#

phaeo_pool <- phaeo_occs;
phaeo_pool$ID <- 1:nrow(phaeo_occs);

phaeo_sampling_20reps <- as.data.frame(matrix(F, nrow = nrow(phaeo_pool), ncol = 20));

for (i in 1:20) {
	
	phaeo_randomized <- phaeo_pool[sample(nrow(phaeo_pool)), ];
	phaeo_randomized_noDuplicates <- phaeo_randomized[!(duplicated(phaeo_randomized$cell)), ];
	phaeo_randomized_IDs <- phaeo_randomized_noDuplicates$ID[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8)];
	
	rownames(phaeo_randomized_noDuplicates[1:round(nrow(phaeo_randomized_noDuplicates) * 0.8), ]);
	
	phaeo_sampling_20reps[phaeo_pool$ID %in% phaeo_randomized_IDs, i] <- T;

}

# get the environmental datasets;

current_buffer12km_nereo <- current_buffer12km[[final_var_nereo]];
LGM_ensemble_buffer12km_nereo <- LGM_ensemble_buffer12km[[final_var_nereo]];
LGM_ccsm_buffer12km_nereo <- LGM_ccsm_buffer12km[[final_var_nereo]];

# set up a dataframe for the validation results;

val_nereo <- data.frame(rep = 1:20, AUC_training = NA, AUC_testing = NA, TSS_training = NA, TSS_testing = NA, AICc_training = NA, AICc_testing = NA);

# ready to run;

for (i in 1:20) {
		
	##### Prepare the data for each replicate;
	
	dat <- prepareSWD(species = "Nereocystis", env = current_buffer12km_nereo, p = nereo_occs[nereo_sampling_20reps[ , i], 1:2], a = phaeo_occs[phaeo_sampling_20reps[ , i], 1:2]);

	##### Split the data into folds, to use with cross-validation;
	
	folds <- randomFolds(dat, k = 5, only_presence = F); # use only_presence = F to also permute the background observations;

	##### Train the model (i.e., run the model) with maxnet;

	model <- train(method = "Maxnet", data = dat, folds = folds, reg = best_hyperparam$reg, fc = best_hyperparam$fc);

	##### get the projections on a map;

	map_current_buffer12km <- predict(model, data = current_buffer12km_nereo, type = "cloglog");
	map_LGM_ensemble_buffer12km <- predict(model, data = LGM_ensemble_buffer12km_nereo, type = "cloglog");
	map_LGM_ccsm_buffer12km <- predict(model, data = LGM_ccsm_buffer12km_nereo, type = "cloglog");

	##### save the projections;
	
	writeRaster(map_current_buffer12km, paste0("Model_output_241113/Nereocystis/projections/current_", i, ".tif"), overwrite = T);
	writeRaster(map_LGM_ensemble_buffer12km, paste0("Model_output_241113/Nereocystis/projections/LGM_ensemble_", i, ".tif"), overwrite = T);
	writeRaster(map_LGM_ccsm_buffer12km, paste0("Model_output_241113/Nereocystis/projections/LGM_ccsm_", i, ".tif"), overwrite = T);

	##### get the model evaluation metrics on the training and test data - note that AICc doesn't run on cross-validation models (and doesn't seem relevant in this situation);
	
	val_nereo$AUC_training[i] <- auc(model);
	val_nereo$AUC_testing[i] <- auc(model, test = T);
	
	val_nereo$TSS_training[i] <- tss(model);
	val_nereo$TSS_testing[i] <- tss(model, test = T);
	
	##### evaluate variable importance and variable response curves;
	
	if (i == 1) {
		pi_nereo <- varImp(model, permut = 20)[ , 1:2];
		colnames(pi_nereo) <- c("Variable", "PI_1");
		#
		jk_nereo <- doJk(model, metric = "auc"); # no test data are specified with a cross-validation model;
		colnames(jk_nereo) <- c("Variable", "Train_AUC_without_1", "Train_AUC_withonly_1");
	} else {
		new_pi_nereo <- varImp(model, permut = 20)[ , 1:2];
		colnames(new_pi_nereo) <- c("Variable", paste0("PI_", i));
		pi_nereo <- merge(pi_nereo, new_pi_nereo);
		#
		new_jk_nereo <- doJk(model, metric = "auc"); # no test data are specified with a cross-validation model;
		colnames(new_jk_nereo) <- c("Variable", paste0("Train_AUC_without_", i), paste0("Train_AUC_withonly_", i));
		jk_nereo <- merge(jk_nereo, new_jk_nereo);
	}
	
	for (j in 1:length(final_var_nereo)) {
		outname <- paste0("Model_output_241113/Nereocystis/response_curves/current_", final_var_nereo[j], "_", i, ".pdf")
		pdf(outname);
		print(plotResponse(model, var = final_var_nereo[j], type = "cloglog", only_presence = T, marginal = F, rug = T));
		dev.off();
	}

	##### get common thresholds for presence-absence;
	
	# note that the thresholds() function requires a single model, which does not work with a SDMmodelCV object - instead, we can get the five individual models of SDMmodel objects and calculate thresholds individually on each of them, then get the mean (the mean will correspond to and be compatible with the mean predicted habitat suitability from the SDMmodelCV, which is simply the mean of the five model predictions);
	
	ths_1 <- thresholds(model@models[[1]], type = "cloglog");
	ths_2 <- thresholds(model@models[[2]], type = "cloglog");
	ths_3 <- thresholds(model@models[[3]], type = "cloglog");
	ths_4 <- thresholds(model@models[[4]], type = "cloglog");
	ths_5 <- thresholds(model@models[[5]], type = "cloglog");
	
	ths_combined <- cbind(ths_1[1:2], ths_2[2], ths_3[2], ths_4[2], ths_5[2]);
				
	if (i == 1) {
		thresholds_nereo <- data.frame(rep = 1, MTP = mean(as.numeric(ths_combined[1, 2:6])), ETSS = mean(as.numeric(ths_combined[2, 2:6])), MTSS = mean(as.numeric(ths_combined[3, 2:6])));
	} else {
		thresholds_nereo <- rbind(thresholds_nereo, c(i, mean(as.numeric(ths_combined[1, 2:6])), mean(as.numeric(ths_combined[2, 2:6])), mean(as.numeric(ths_combined[3, 2:6]))));
	}

}

write.table(thresholds_nereo, "Model_output_241113/Nereocystis/results_presenceAbsence_thresholds.txt", sep = "\t", row.names = F, col.names = T, quote = F);

write.table(val_nereo, "Model_output_241113/Nereocystis/results_modelMetrics.txt", sep = "\t", row.names = F, col.names = T, quote = F);

#

pi_nereo$PI_mean <- apply(pi_nereo[ , 2:21], 1, mean);
pi_nereo$PI_sd <- apply(pi_nereo[ , 2:21], 1, sd);

jk_nereo$Train_AUC_without_mean <- apply(jk_nereo[ , paste0("Train_AUC_without_", 1:20)], 1, mean);
jk_nereo$Train_AUC_without_sd <- apply(jk_nereo[ , paste0("Train_AUC_without_", 1:20)], 1, sd);

jk_nereo$Train_AUC_withonly_mean <- apply(jk_nereo[ , paste0("Train_AUC_withonly_", 1:20)], 1, mean);
jk_nereo$Train_AUC_withonly_sd <- apply(jk_nereo[ , paste0("Train_AUC_withonly_", 1:20)], 1, sd);

write.table(pi_nereo, "Model_output_241113/Nereocystis/results_permutation_importance.txt", sep = "\t", row.names = F, col.names = T, quote = F);
write.table(jk_nereo, "Model_output_241113/Nereocystis/results_jackknife.txt", sep = "\t", row.names = F, col.names = T, quote = F);

# get a single final raster, plus uncertainty, and make a map;

nereo_stack_current_files <- paste0("Model_output_241113/Nereocystis/projections/", list.files("Model_output_241113/Nereocystis/projections/", pattern = "current"));
nereo_stack_current <- rast(nereo_stack_current_files);

nereo_current_mean_projection <- mean(nereo_stack_current);

plot(nereo_current_mean_projection);
plot(nereo_current_mean_projection, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_current_mean_projection, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

# 

nereo_stack_LGM_ensemble_files <- paste0("Model_output_241113/Nereocystis/projections/", list.files("Model_output_241113/Nereocystis/projections/", pattern = "LGM_ensemble"));
nereo_stack_LGM_ensemble <- rast(nereo_stack_LGM_ensemble_files);

nereo_LGM_ensemble_mean_projection <- mean(nereo_stack_LGM_ensemble);

plot(nereo_LGM_ensemble_mean_projection);
plot(nereo_LGM_ensemble_mean_projection, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_LGM_ensemble_mean_projection, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

nereo_stack_LGM_ccsm_files <- paste0("Model_output_241113/Nereocystis/projections/", list.files("Model_output_241113/Nereocystis/projections/", pattern = "LGM_ccsm"));
nereo_stack_LGM_ccsm <- rast(nereo_stack_LGM_ccsm_files);

nereo_LGM_ccsm_mean_projection <- mean(nereo_stack_LGM_ccsm);

plot(nereo_LGM_ccsm_mean_projection);
plot(nereo_LGM_ccsm_mean_projection, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_LGM_ccsm_mean_projection, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

# plot presence/absence using:
	# MTP - minimum training presence;
	# ETSS - equal training sensitivity and specificity;
	# MTSS - maximum training sensitivity plus specificity;
	# Q05 - 5 percent quantile (my own method);

thresholds_nereo <- read.table("Model_output_241113/Nereocystis/results_presenceAbsence_thresholds.txt", sep = "\t", header = T, stringsAsFactors = F);

#

nereo_stack_current_mtp <- nereo_stack_current;
for (i in 1:20) {
	nereo_stack_current_mtp[[i]] <- nereo_stack_current_mtp[[i]] >= thresholds_nereo$MTP[i];
}
nereo_current_mtp_sum <- sum(nereo_stack_current_mtp);

plot(nereo_current_mtp_sum);
plot(nereo_current_mtp_sum, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_current_mtp_sum, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

plot(nereo_current_mean_projection > mean(thresholds_nereo$MTP)); # alternative method;
plot(nereo_current_mean_projection > mean(thresholds_nereo$MTP), xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_current_mean_projection > mean(thresholds_nereo$MTP), xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

nereo_stack_current_mtss <- nereo_stack_current;
for (i in 1:20) {
	nereo_stack_current_mtss[[i]] <- nereo_stack_current_mtss[[i]] >= thresholds_nereo$MTSS[i];
}
nereo_current_mtss_sum <- sum(nereo_stack_current_mtss);

plot(nereo_current_mtss_sum);
plot(nereo_current_mtss_sum, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_current_mtss_sum, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

plot(nereo_current_mean_projection > mean(thresholds_nereo$MTSS)); # alternative method;
plot(nereo_current_mean_projection > mean(thresholds_nereo$MTSS), xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_current_mean_projection > mean(thresholds_nereo$MTSS), xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));

#

nereo_stack_current_q05 <- nereo_stack_current;
for (i in 1:20) {
	nereo_stack_current_q05[[i]] <- nereo_stack_current_q05[[i]] >= quantile(unlist(extract(nereo_stack_current[[i]], unique(nereo_occs$cell))), 0.05);
}
nereo_current_q05_sum <- sum(nereo_stack_current_q05);

plot(nereo_current_q05_sum);
plot(nereo_current_q05_sum, xlim = c(-3e+06, 2e+06), ylim = c(-3e+06, 3e+06));
plot(nereo_current_q05_sum, xlim = c(-1e+06, 0.5e+06), ylim = c(0.5e+06, 2.5e+06));
