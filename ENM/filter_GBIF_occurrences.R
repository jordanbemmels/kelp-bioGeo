##### Code illustrates how raw GBIF occurrence records were filtered down to the set of high-quality filtered records subsequently used for construction ecological niche models;

# Author use only: GitHub version derived from custom ENMs script filter_GBIF.R, 2025/12/03;

#####

### Prerequisites:

# User has already downloaded raw occurrence records for Macrocystis, Nereocystis, and brown algae (Phaeophyceae) from GBIF and saved them as three separate tab-separated .txt files;

# Citations for the raw occurrence record downloads:

# GBIF.org. 2024a. ‘Macrocystis’ GBIF Occurrence Download, 10 October 2024. https://doi.org/10.15468/dl.pvf237.

# GBIF.org. 2024b. ‘Nereocystis’ GBIF Occurrence Download, 10 October 2024. https://doi.org/10.15468/dl.rebpy7.

# GBIF.org. 2024c. ‘Phaeophyceae’ GBIF Occurrence Download, 10 October 2024. https://doi.org/10.15468/dl.n8evsc.

#####

macro <- read.table("GBIF/Macrocystis_global_241010.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", comment = "");
nereo <- read.table("GBIF/Nereocystis_global_241010.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", comment = "");
phaeo <- read.table("GBIF/Phaeophyceae_NEPacific_241010.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", comment = "");

nrow(macro); # 93793;
nrow(nereo); # 41622;
nrow(phaeo); # 568876;

##############################
##############################

##### Reduce macro and nereo to the same geographic subset as phaeo;

# from my notes about occurrence record downloads, for Phaeophyceae the range of coordinates was from 13.5 to 72.5º latitude and -180.0º to -91.0º longitude;

macro <- macro[!(is.na(macro$decimalLatitude)) & !(is.na(macro$decimalLongitude)) & macro$decimalLatitude >= 13.5 & macro$decimalLatitude <= 72.5 & macro$decimalLongitude >= -180.0 & macro$decimalLongitude <= -91.0, ];
nereo <- nereo[!(is.na(nereo$decimalLatitude)) & !(is.na(nereo$decimalLongitude)) & nereo$decimalLatitude >= 13.5 & nereo$decimalLatitude <= 72.5 & nereo$decimalLongitude >= -180.0 & nereo$decimalLongitude <= -91.0, ];
phaeo <- phaeo[!(is.na(phaeo$decimalLatitude)) & !(is.na(phaeo$decimalLongitude)) & phaeo$decimalLatitude >= 13.5 & phaeo$decimalLatitude <= 72.5 & phaeo$decimalLongitude >= -180.0 & phaeo$decimalLongitude <= -91.0, ];

nrow(macro); # 79008;
nrow(nereo); # 41242;
nrow(phaeo); # 568876 - as expected, no change;

plot(macro$decimalLongitude, macro$decimalLatitude, pch = ".");
plot(nereo$decimalLongitude, nereo$decimalLatitude, pch = ".");
plot(phaeo$decimalLongitude, phaeo$decimalLatitude, pch = ".");

##############################
##############################

##### Perform additional filtering;

colnames(macro);

unique(macro$genus);
unique(nereo$genus);

unique(macro$species);
unique(nereo$species);

unique(macro$scientificName);
unique(nereo$scientificName);

unique(macro$verbatimScientificName);
unique(nereo$verbatimScientificName);

unique(c(macro$occurrenceStatus, nereo$occurrenceStatus, phaeo$occurrenceStatus));

#####

unique(c(macro$basisOfRecord, nereo$basisOfRecord, phaeo$basisOfRecord));

tapply(macro$basisOfRecord, macro$basisOfRecord, length);
tapply(nereo$basisOfRecord, nereo$basisOfRecord, length);

# descriptions of GBIF basisOfRecord categories: https://docs.gbif.org/course-data-use/en/basis-of-record.html;

# most of our records are HUMAN_OBSERVATION, OCCURRENCE, or PRSERVED_SPECIMEN;

# HUMAN_OBSERVATION are mostly records such as iNaturalist, and don't have evidence recorded;
unique(macro$institutionCode[macro$basisOfRecord == "HUMAN_OBSERVATION"]);
unique(nereo$institutionCode[nereo$basisOfRecord == "HUMAN_OBSERVATION"]);

# OCCURRENCE mostly do not have any institution code whatsoever - note that this typically means there is no further info available: https://discourse.gbif.org/t/basis-of-record-ocurrence/3269;
unique(macro$institutionCode[macro$basisOfRecord == "OCCURRENCE"]);
unique(nereo$institutionCode[nereo$basisOfRecord == "OCCURRENCE"]);

# the only trustworthy basis here is PRESERVED_SPECIMEN, although note that it may also be (rarely) possible for a preserved specimen to have arisen itself from a zoo or garden, etc.;

macro <- macro[macro$basisOfRecord == "PRESERVED_SPECIMEN", ];
nereo <- nereo[nereo$basisOfRecord == "PRESERVED_SPECIMEN", ];
phaeo <- phaeo[phaeo$basisOfRecord == "PRESERVED_SPECIMEN", ];

nrow(macro);
nrow(nereo);
nrow(phaeo);

#####

# filter the data to remove those with recognized severe issues, as well as low coordinate precision;

# for kelp here, I am not considering "TAXON_MATCH_FUZZY" to be a problem, because there are multiple versions of Macrocystis's name due to taxonomic uncertainty, and I *do* want to include anything Macrocystis, including other scientific names or misspellings;
# for example, we have "Macrocystis integrifolia" that is interpreted as "Macrocystis integrifolius Bory" which is a fuzzy taxonomic match, but clearly not a problem;
# instead, I have checked that the verbatimScientiticName category above is reasonable for all macro and nereo individuals;

# include "CONTINENT_COORDINATE_MISMATCH" as I discovered there were some kelp sample in the ocean that made not sense and which had this designation - and I had already included COUNTRY_COORDINATE_MISTMATCH so this makes logical sense to include;

# require at least two decimals of precision to keep a record - I don't think one is enough as the record could still be off by ~10 km (111 km per degree at the equator * 0.1 degrees = 11.1 km);

# also use the GBIF-provided coordinateUncertaintyInMeters, and set this to a maximum of 1,000 m;

filterData <- function(occData) {
	
	### remove records with severe issues - issues that would indicate something has gone terribly wrong with identifying the latitude and longitude of the specimen in particular - other issues are minor or not related to the lat/long so they can be kept for now;
	
	severeIssues = c("CONTINENT_COORDINATE_MISMATCH", "CONTINENT_COUNTRY_MISMATCH", "COORDINATE_INVALID", "COORDINATE_REPROJECTION_FAILED", "COORDINATE_REPROJECTION_SUSPICIOUS", "COUNTRY_COORDINATE_MISMATCH", "GEODETIC_DATUM_INVALID", "TAXON_MATCH_HIGHERRANK", "TAXON_MATCH_NONE", "ZERO_COORDINATE");
	
	occData$severeIssue <- 0;
	for (i in 1:nrow(occData)) {
		if (sum(strsplit(occData[i, "issue"], ";")[[1]] %in% severeIssues) >= 1) {
			occData$severeIssue[i] <- 1;
		} 
	}
	
	occData <- occData[occData$severeIssue == 0, ];
				
	### remove low-precision records;

	# first identify records with few decimal places;
	# see https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r;
	decimalplaces <- function(x) {
		stopifnot(class(x)=="numeric")
		x <- sub("0+$","",x)
		x <- sub("^.+[.]","",x)
		nchar(x)
	}
	
	occData$latPrecision <- decimalplaces(occData$decimalLatitude);
	occData$longPrecision <- decimalplaces(occData$decimalLongitude);

	# also identify records that round to 0, 10, 20, 30, 40, 50, 60 minutes, or to 0, 15, 30, 45, 60 minutes, as these are likely to be imprecise coordinates recorded in a format like 12º10, 17º50, etc., that are likely approximations rather than actual coordinates;

	occData$latTenRemainder <- occData$decimalLatitude %% (10/60);
	occData$longTenRemainder <- occData$decimalLongitude %% (10/60);
	
	occData$latQuaterRemainder <- occData$decimalLatitude %% (15/60);
	occData$longQuaterRemainder <- occData$decimalLongitude %% (15/60);

	# need to have a tolerance for the 10/60 remainders because 0.66 and 0.67 are both bad rounding, so if the remainder is LESS THAN OR EQUAL TO  0.01 then it is a problem, but > 0.01 indicates the coordinate is fine - don't need a tolerance for the 15/60 remainders because these are 0.25, 0.50, 0.75 and should be recorded exactly;

	tenRemainderTolerance <- 0.01; 

	# now, check which records have problems;

	occData$precisionError <- rep(0, nrow(occData));
	occData$tenRemainderError <- rep(0, nrow(occData));
	occData$quarterRemainderError <- rep(0, nrow(occData));

	for (i in 1:nrow(occData)) {
	
		# note that we use & rather than | because if one of the lat or long does not return an precision or remainder error, this suggests that the coordinates in general are fine and the potential error was simply due to chance;
	
		if ((occData$latPrecision[i] < 2) & (occData$longPrecision[i]) < 2) {
			occData$precisionError[i] <- 1;
		}
		
		if ((occData$latTenRemainder[i] <= tenRemainderTolerance) & (occData$longTenRemainder[i]) <= tenRemainderTolerance) {
			occData$tenRemainderError[i] <- 1;
		}
	
		if ((occData$latQuaterRemainder[i] == 0) & (occData$longQuaterRemainder[i] == 0)) {
			occData$quarterRemainderError[i] <- 1;
		}
	
	}

	# check that it worked;	
	occData[ , c("decimalLatitude", "decimalLongitude", "latPrecision", "longPrecision", "precisionError")];	
	occData[ , c("decimalLatitude", "decimalLongitude", "latTenRemainder", "longTenRemainder", "tenRemainderError")];	
	occData[ , c("decimalLatitude", "decimalLongitude", "latQuaterRemainder", "longQuaterRemainder", "quarterRemainderError")];
	
	# finally, remove the offending records;

	occData$toRemove <- rep(0, nrow(occData));
	
	for (i in 1:nrow(occData)) {
		
		if ((occData$precisionError[i] == 1) | (occData$tenRemainderError[i] == 1) | (occData$quarterRemainderError[i] == 1) | (!(is.na(occData$coordinateUncertaintyInMeters[i])) & occData$coordinateUncertaintyInMeters[i] >= 1000)) {
			occData$toRemove[i] <- 1;
		}
		
	}
	
	occData <- occData[occData$toRemove != 1, ];
		
}

macro_filtered <- filterData(macro);
nereo_filtered <- filterData(nereo);
phaeo_filtered <- filterData(phaeo);

nrow(macro_filtered); #1269;
nrow(nereo_filtered); # 669;
nrow(phaeo_filtered); # 45,118;

plot(macro$decimalLongitude, macro$decimalLatitude, pch = 21, col = "red");
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "black");

plot(nereo$decimalLongitude, nereo$decimalLatitude, pch = 21, col = "red");
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "black");

plot(phaeo$decimalLongitude, phaeo$decimalLatitude, pch = 21, col = "red");
points(phaeo_filtered$decimalLongitude, phaeo_filtered$decimalLatitude, pch = 21, col = "black");

# optional - view with map;

require(raster);
NAmerica <- shapefile("NA_PoliticalDivisions/data/bound_p/boundary_p_v2.shp");

basic_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";

NAmerica_basic <- spTransform(NAmerica, basic_crs);

plot(NAmerica_basic[!(NAmerica_basic$COUNTRY == "water/agua/d'eau"), ], col = "grey", xlim = c(-180, -91));
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[!(NAmerica_basic$COUNTRY == "water/agua/d'eau"), ], col = "grey", xlim = c(-180, -91));
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[!(NAmerica_basic$COUNTRY == "water/agua/d'eau"), ], col = "grey", xlim = c(-180, -91));
points(phaeo_filtered$decimalLongitude, phaeo_filtered$decimalLatitude, pch = 21, col = "red");

#########################
#########################

##### there are still some occurrences that don't make any sense for each species, so inspect or remove these;

##### ALASKA - NON-PANHANDLE #####

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-AK" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey", xlim = c(-180, -150));
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-AK" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey", xlim = c(-180, -150));
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

### samples in the Pribilof Islands, Alaska - outside the known range of either species;

macro_filtered[macro_filtered$decimalLongitude < -165, ];

# the Macrocystis sample is marked on the original herbarium image as "Macrocystis?" - dubious - best to remove;
# image: https://sweetgum.nybg.org/science/vh/specimen-details/?irn=3436110

nereo_filtered[nereo_filtered$decimalLongitude < -165, ];

# there are two samples, one marked on the original herbarium image with a "?" but the other does not - both samples are just a piece of a blade with no distinguishing characteristics - dubious - best to remove;
# image: https://sweetgum.nybg.org/science/vh/specimen-details/?irn=3444946
# image: https://sweetgum.nybg.org/science/vh/specimen-details/?irn=3444947

macro_filtered <- macro_filtered[!(macro_filtered$decimalLongitude < -165), ];
nereo_filtered <- nereo_filtered[!(nereo_filtered$decimalLongitude < -165), ];

##### ALASKA - PANHANDLE #####

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-AK" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey", xlim = c(-150, -130), ylim = c(54, 62));
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-AK" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey", xlim = c(-150, -130), ylim = c(54, 62));
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

# no problems;

##### BC #####

plot(NAmerica_basic[NAmerica_basic$STATEABB == "CA-BC" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[NAmerica_basic$STATEABB == "CA-BC" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

### for Macrocystis, note that the one off NW Vancouver Island is likely valid - location of Triangle Island;

# there is one of each both species located in the middle of Vancouver Island;

macro_filtered[macro_filtered$decimalLongitude > -124.884942 & macro_filtered$decimalLongitude < -123.85093 & macro_filtered$decimalLatitude > 48.70 & macro_filtered$decimalLatitude < 49.1, ];

nereo_filtered[nereo_filtered$decimalLongitude > -124.884942 & nereo_filtered$decimalLongitude < -123.85093 & nereo_filtered$decimalLatitude > 48.70 & nereo_filtered$decimalLatitude < 49.1, ];

# this Macrocystis one had its coordinates transcribed incorrectly on GBIF! they do not match the herbarium sheet - best to remove, as I don't have the ability to hand-check coordinates of all other samples;
# https://collections.beatymuseum.ubc.ca/specimen/search?catalogNumber=A042169&entity=1987713433

# for Nereocystis, there is no herbarium image available, but it's the same site description and appears to be the same coordinates rounding problem as for Macrocystis - also remove;
# https://collections.beatymuseum.ubc.ca/specimen/search?q=A042168

macro_filtered <- macro_filtered[!(macro_filtered$decimalLongitude > -124.884942 & macro_filtered$decimalLongitude < -123.85093 & macro_filtered$decimalLatitude > 48.70 & macro_filtered$decimalLatitude < 49.1), ];

nereo_filtered <- nereo_filtered[!(nereo_filtered$decimalLongitude > -124.884942 & nereo_filtered$decimalLongitude < -123.85093 & nereo_filtered$decimalLatitude > 48.70 & nereo_filtered$decimalLatitude < 49.1), ];

### there is a Macrocystis is Nanimo, which is outside its known range;

macro_filtered[macro_filtered$decimalLongitude > -125 & macro_filtered$decimalLatitude > 49, ];

# indeed, the specimen apperas to be identified by the collectors as a Macrocystis and is substantial with morphological details (not just a rectangular blade cutting) and is from "Departure Bay, BC" (i.e., Nanaimo), and is from 1887 - however, this is totally unreasonable given the known distribution of Macrocystis for it to have occurred here;
# see also Rigg 1913. THE DISTRIBUTION OF MACROCYSTIS PYRIFERA ALONG THE AMERICAN SHORE OF THE STRAIT OF JUAN DE FUCA. (https://www.jstor.org/stable/40595390?seq=2) who clearly state that extensive surveys have NOT found any Macrocystis in the Strait of Georgia or the San Juan Islands or the innermost reaches of Strait of Juan de Fuca (this is all referring to US locations, though) - it would seem totally unreasonable for this individual to be alive and growing in Departure Bay as it's so far outside the range of the species - remove this specimen;

macro_filtered <- macro_filtered[!(macro_filtered$decimalLongitude > -125 & macro_filtered$decimalLatitude > 49), ];

### there are several Macrocystis in Victoria and even farther east in the San Juans - may be historical populations - investigate;

macro_filtered[macro_filtered$decimalLongitude > -123.5067 & macro_filtered$decimalLatitude > 48, ];

# they are three historical samples with dates 1908, 1913, 1925;
# see comments above, this is not compatible with Rigg 1913 as cited above;
# this also explicitly states that a sample from Whidby Island (not any of the ones currently in question) was previously used to argue that it occurred in the San Juans, but this is believed erroneous as the collector confirmed that it was found FLOATING and not growing - so, there is at least some precedent for things from this time period appearing in herbaria that are floating;
# plus, Rigg 1913 described distribution is not compatible with Macrocystis occurring in the San Juans (one of the three specimens here), and unlikely to be compatible with Macrocystis occurring in Victoria;
	# Rigg 1913 says the distribution extends to Low Point, Washington, which is across from Sooke, and matches the known distribution well - it is unlikely that in early 1900s the distribution matched the known distribution from today on the American side, but extended much farther eastward into the strait on the Canadian side up to Victoria - remove these specimens too from Victoria;

macro_filtered <- macro_filtered[!(macro_filtered$decimalLongitude > -123.5067 & macro_filtered$decimalLatitude > 48), ];

##### WASHINGTON #####

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-WA" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-WA" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

# no problems;

##### OREGON #####

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-OR" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-OR" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

### does Macrocystis grow in Oregon naturally?;

macro_filtered[macro_filtered$decimalLatitude > 42 & macro_filtered$decimalLatitude <46, ];

# https://hdcgcx1.deq.state.or.us/Geocortex/Essentials/REST/sites/GRP_Strategies/map/mapservices/10/layers/330
# "In Oregon, the giant kelp, Macrocystis sp. is found only at Cape Arago."

# most of the samples are clustered around Cape Arago, there is one clustered just south a bit (43.11900,  -124.4084) that appears to be near "Cat and Kitten Rocks" which does look like giant kelp habitat on Google Maps - retain this record too;

# this means that most of our locations are likely correct but one may by suspcicious - inspect;

macro_filtered[macro_filtered$decimalLatitude > 44 & macro_filtered$decimalLatitude <46, ];

# record is from 1937 from "US Marine Gardens, Otter Rock", but there is a Marine Reserve at Otter Rock - it seems that locally these areas are referred to as gardens even though they are nature reserves - not likely a human garden - this is plausible as a historical record, as Googling seems to indicate Otter Rock is a very rocky area, which doesn't currently (in 21st century) have kelp but has substrate specifically identified as possible for kelp, but that it used to be an area with a ton of otters until they were extirpated in the early 20th century - thus it seems very plausible that there was actually kelp growing here at one point;
# also this definitely looks like a real specimen, is high quality, clearly labelled as from Otter Rock: https://oregonflora.org/imglib/OSU_A/OSC-A-015/OSC-A-015504.JPG;
# retain this specimen;

##### CALIFORNIA #####

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-CA" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[NAmerica_basic$STATEABB == "US-CA" & !(is.na(NAmerica_basic$STATEABB)), ], col = "grey");
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

### a Macrocystis inside of San Francisco Bay seems dubious - investigate;

macro_filtered[macro_filtered$decimalLongitude > -122.4528 & macro_filtered$decimalLatitude > 37.5294, ];

# search here https://webapps.cspace.berkeley.edu/ucjeps/publicsearch/publicsearch/
# Specimen IDs UC1983242, UC1716049, UC1716050;
# these records do not have any unusual info to trigger me to be suspicious, other than that they are collected at Pier 39, San Francisco, in a highly urban place - though we also have urban kelp in Vancouver so perhaps not cause for concern;

### several Macrocystis in the middle of the ocean;

macro_filtered[macro_filtered$decimalLongitude < -122.3644 & macro_filtered$decimalLatitude < 36.8809, ];

# for all of these three, the locality descriptions do not match the latitude and longitude - clearly some sort of transcription or recording error - remove;

macro_filtered <- macro_filtered[!(macro_filtered$decimalLongitude < -122.3644 & macro_filtered$decimalLatitude < 36.8809), ];

### a Macrocystis in the middle of land;

macro_filtered[macro_filtered$decimalLongitude > -120.2741 & macro_filtered$decimalLatitude > 34.8862, ];

# no locality, date is 1895, clearly something wrong - remove;

macro_filtered <- macro_filtered[!(macro_filtered$decimalLongitude > -120.2741 & macro_filtered$decimalLatitude > 34.8862), ];

### there are a few more Macrocystis that appear to be in the middle of the ocean, but these are likely Santa Barbara Island (very tiny and doesn't show up on map);

### multiple Nereocystis in southern California that may(?) be beyond its southermost distributional limits;

# shapefile from surveys of distribution: https://catalog.data.gov/dataset/kelp-distribution-off-california1;
# annoying, file does not exist!;

# several internet sites seem to say the southern limit is either Point Conception, or San Luis Obispo County (just a bit north of Point Conception), but I can't find any verifiable data or primary literature;

# aha! this one seems fairly reliable and would be written by bull kelp experts:
# https://bullkelp.info/regions
# indeed, shows the distribution ending at San Luis Obispo County;

nereo_filtered[nereo_filtered$decimalLatitude < 34.5, ];

# aha! a great many of them are described as "cast ashore" or "beach" or "cultured" or "on the beach";
# these appear to be NOT NATIVELY GROWING IN AREA - remove all;

nereo_filtered <- nereo_filtered[!(nereo_filtered$decimalLatitude < 34.5), ];

##### MEXICO #####

plot(NAmerica_basic[NAmerica_basic$COUNTRY == "MEX" & !(is.na(NAmerica_basic$COUNTRY)), ], col = "grey");
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[NAmerica_basic$COUNTRY == "MEX" & !(is.na(NAmerica_basic$COUNTRY)), ], col = "grey");
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

### one Macrocystis sample far inland;

macro_filtered[macro_filtered$decimalLongitude > -115.2738 & macro_filtered$decimalLatitude > 29.8121, ];

# locality simply given as Baja California and seems to be georeferenced to the middle of the state - remove;

macro_filtered <- macro_filtered[!(macro_filtered$decimalLongitude > -115.2738 & macro_filtered$decimalLatitude > 29.8121), ];

### one Macrocystis in the middle of the ocean;

macro_filtered[macro_filtered$decimalLatitude > 23 & macro_filtered$decimalLatitude < 27, ];

# WOW this is  REAL RECORD! locality is "Alijos rocks, off Baja California", which is a SUPER TINY little series of rocks that really do exist! - retain record at 24.95, -115.7333!;

### one Macrocystis very far south on Socorro Island;

# this is likely a real record, too:
# https://www.nature.com/articles/s41598-023-38944-7 "still occurs in elusive habitats (e.g., Socorro Island in Mexico)92"
# 92 is Taylor, W. R. Pacific marine algae of the Allan Hancock Expeditions to the Galapagos Islands. (Allan Hancock Pacific Expeditions, 1945).;

#########################
#########################

##### Remove areas from Phaeophyceae background that we don't want to include;
##### UPDATE 2024/11/13 - NEW FILTERING: previously I had excluded Hawaii and the Arctic Coast of Alaska - now I am retaining Hawaii and the Arctic Coast of Alaska to help better constrain the models in very cold and very warm waters;

### our study area should be the Northeast Pacific, in Canada, USA, and Mexico;

plot(NAmerica_basic[!(NAmerica_basic$COUNTRY == "water/agua/d'eau"), ], col = "grey", xlim = c(-180, -91));
points(phaeo_filtered$decimalLongitude, phaeo_filtered$decimalLatitude, pch = 21, col = "red");

### remove samples from the Canadian Arctic only (biogeographically too distant);

phaeo_filtered <- phaeo_filtered[!(phaeo_filtered$decimalLatitude > 65.60998 & phaeo_filtered$decimalLongitude > -141), ]

### remove samples from the Atlantic Ocean;

phaeo_filtered <- phaeo_filtered[!(phaeo_filtered$decimalLatitude > 17.73320 & phaeo_filtered$decimalLongitude > -98.4590), ];

### UPDATE: do NOT remove samples from Hawaii - even though this is a different biogeographic region from the main coast of North America, we want to help constrain that it is not found in warm tropical waters, and this falls within the lat/long of the study area;

#phaeo_filtered <- phaeo_filtered[!(phaeo_filtered$decimalLatitude < 34 & phaeo_filtered$decimalLongitude < -147), ];

#########################
#########################

##### RECONFIRM RANGE-WIDE SAMPLES #####

plot(NAmerica_basic[!(NAmerica_basic$COUNTRY == "water/agua/d'eau"), ], col = "grey", xlim = c(-180, -91));
points(macro_filtered$decimalLongitude, macro_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[!(NAmerica_basic$COUNTRY == "water/agua/d'eau"), ], col = "grey", xlim = c(-180, -91));
points(nereo_filtered$decimalLongitude, nereo_filtered$decimalLatitude, pch = 21, col = "red");

plot(NAmerica_basic[!(NAmerica_basic$COUNTRY == "water/agua/d'eau"), ], col = "grey", xlim = c(-180, -91));
points(phaeo_filtered$decimalLongitude, phaeo_filtered$decimalLatitude, pch = 21, col = "red");

##################################################
##################################################
##################################################

##### Pre-emptively remove samples from land #####

### There are many samples where the coordinates are on land - some of these are clearly mistakes that I haven't dealt with yet and can simply be removed, others might be from areas that are very close to shore and even slight rounding errors push them onto land;
### We can use the same environmental rasters we plan to use for ENMs to remove samples that fall on land;
### Initially I thought we would also check if coordinates can be "rescued" by picking a surrounding cell, to account for samples that only appear on land due to rounding - however, the raster map favours water over land VERY STRONGLY, so most of the samples appearing on land might be true errors or else quite far from the coast (if coordinates correct) such that the neighbouring pixels may not accurately represent the environmental conditions - also the number of points on land is likely to be very small - probably better practice to simple exclude these points;

### load an example raster;
# raster is originally from the MARSPEC dataset but has been reprojected to an Albers Equal Area projection specially formulated to focus on Vancouver Island: "+proj=aea +lat_0=40 +lon_0=-125 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";

bathy <- raster("current_5m_vi_aea/bathy_5m.tif");

vi_aea_crs <- crs(bathy);

### convert sample points to the projection of bathy;

cds_macro_filtered <- macro_filtered[ , c("decimalLongitude", "decimalLatitude")];
coordinates(cds_macro_filtered) <- c("decimalLongitude", "decimalLatitude");
proj4string(cds_macro_filtered) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";

cds_macro_filtered_vi_aea <- spTransform(cds_macro_filtered, vi_aea_crs);
macro_filtered$Long_vi_aea <- cds_macro_filtered_vi_aea@coords[ , 1];
macro_filtered$Lat_vi_aea <- cds_macro_filtered_vi_aea@coords[ , 2];

plot(bathy);
points(macro_filtered$Long_vi_aea, macro_filtered$Lat_vi_aea);

#

cds_nereo_filtered <- nereo_filtered[ , c("decimalLongitude", "decimalLatitude")];
coordinates(cds_nereo_filtered) <- c("decimalLongitude", "decimalLatitude");
proj4string(cds_nereo_filtered) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";

cds_nereo_filtered_vi_aea <- spTransform(cds_nereo_filtered, vi_aea_crs);
nereo_filtered$Long_vi_aea <- cds_nereo_filtered_vi_aea@coords[ , 1];
nereo_filtered$Lat_vi_aea <- cds_nereo_filtered_vi_aea@coords[ , 2];

plot(bathy);
points(nereo_filtered$Long_vi_aea, nereo_filtered$Lat_vi_aea);

#

cds_phaeo_filtered <- phaeo_filtered[ , c("decimalLongitude", "decimalLatitude")];
coordinates(cds_phaeo_filtered) <- c("decimalLongitude", "decimalLatitude");
proj4string(cds_phaeo_filtered) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";

cds_phaeo_filtered_vi_aea <- spTransform(cds_phaeo_filtered, vi_aea_crs);
phaeo_filtered$Long_vi_aea <- cds_phaeo_filtered_vi_aea@coords[ , 1];
phaeo_filtered$Lat_vi_aea <- cds_phaeo_filtered_vi_aea@coords[ , 2];

plot(bathy);
points(phaeo_filtered$Long_vi_aea, phaeo_filtered$Lat_vi_aea);

### remove points on land;

macro_filtered$isOceanPixel <- extract(bathy, macro_filtered[ , c("Long_vi_aea", "Lat_vi_aea")]);
macro_filtered$isOceanPixel[!(is.na(macro_filtered$isOceanPixel))] <- 1;
macro_filtered$isOceanPixel[is.na(macro_filtered$isOceanPixel)] <- 0;

plot(macro_filtered$Long_vi_aea[macro_filtered$isOceanPixel == 0], macro_filtered$Lat_vi_aea[macro_filtered$isOceanPixel == 0]);
plot(bathy, add = T);
points(macro_filtered$Long_vi_aea[macro_filtered$isOceanPixel == 0], macro_filtered$Lat_vi_aea[macro_filtered$isOceanPixel == 0]);

sum(macro_filtered$isOceanPixel == 0)  / nrow(macro_filtered); # 1.0% of records;

#

nereo_filtered$isOceanPixel <- extract(bathy, nereo_filtered[ , c("Long_vi_aea", "Lat_vi_aea")]);
nereo_filtered$isOceanPixel[!(is.na(nereo_filtered$isOceanPixel))] <- 1;
nereo_filtered$isOceanPixel[is.na(nereo_filtered$isOceanPixel)] <- 0;

plot(nereo_filtered$Long_vi_aea[nereo_filtered$isOceanPixel == 0], nereo_filtered$Lat_vi_aea[nereo_filtered$isOceanPixel == 0]);
plot(bathy, add = T);
points(nereo_filtered$Long_vi_aea[nereo_filtered$isOceanPixel == 0], nereo_filtered$Lat_vi_aea[nereo_filtered$isOceanPixel == 0]);

sum(nereo_filtered$isOceanPixel == 0)  / nrow(nereo_filtered); # 0.3% of records;

#

phaeo_filtered$isOceanPixel <- extract(bathy, phaeo_filtered[ , c("Long_vi_aea", "Lat_vi_aea")]);
phaeo_filtered$isOceanPixel[!(is.na(phaeo_filtered$isOceanPixel))] <- 1;
phaeo_filtered$isOceanPixel[is.na(phaeo_filtered$isOceanPixel)] <- 0;

plot(phaeo_filtered$Long_vi_aea[phaeo_filtered$isOceanPixel == 0], phaeo_filtered$Lat_vi_aea[phaeo_filtered$isOceanPixel == 0]);
plot(bathy, add = T);
points(phaeo_filtered$Long_vi_aea[phaeo_filtered$isOceanPixel == 0], phaeo_filtered$Lat_vi_aea[phaeo_filtered$isOceanPixel == 0]);

sum(phaeo_filtered$isOceanPixel == 0)  / nrow(phaeo_filtered); # 1.7% of records;

#

macro_filtered <- macro_filtered[macro_filtered$isOceanPixel == 1, ];
nereo_filtered <- nereo_filtered[nereo_filtered$isOceanPixel == 1, ];
phaeo_filtered <- phaeo_filtered[phaeo_filtered$isOceanPixel == 1, ];

nrow(macro_filtered); #1,245;
nrow(nereo_filtered); #647
nrow(phaeo_filtered); #43,792;

##################################################
##################################################
##################################################

##### View the final product and save;

pdf("GBIF_filtered/pdfs/Macrocystis_global_241010_filtered.pdf");
plot(bathy, col = "skyblue", legend = F);
points(macro_filtered$Long_vi_aea, macro_filtered$Lat_vi_aea, pch = ".", col = "red");
dev.off();

pdf("GBIF_filtered/pdfs/Nereocystis_global_241010_filtered.pdf");
plot(bathy, col = "skyblue", legend = F);
points(nereo_filtered$Long_vi_aea, nereo_filtered$Lat_vi_aea, pch = ".", col = "red");
dev.off();

pdf("GBIF_filtered/pdfs/Phaeophyceae_global_241010_filtered.pdf");
plot(bathy, col = "skyblue", legend = F);
points(phaeo_filtered$Long_vi_aea, phaeo_filtered$Lat_vi_aea, pch = ".", col = "red");
dev.off();

write.table(macro_filtered, "GBIF_filtered/Macrocystis_global_241010_filtered.txt", sep = "\t", row.names = F, col.names = T, quote = F);
write.table(nereo_filtered, "GBIF_filtered/Nereocystis_global_241010_filtered.txt", sep = "\t", row.names = F, col.names = T, quote = F);
write.table(phaeo_filtered, "GBIF_filtered/Phaeophyceae_global_241010_filtered.txt", sep = "\t", row.names = F, col.names = T, quote = F);
