# Ecological niche models

### filter_GBIF_occurrences.R

This _R_ script illustrates how the raw occurrence records from GBIF (gbif.org) were filtered prior to constructing ecological niche models (ENMs). The code subsets all occurrences to an appropriate extent and background area, removes occurrences flagged with severe issues, removes occurrences with low precision or extreme rounding of latitude/longitude, removes individual records that are suspicious or outside the known range of either species, and removes pixels that fall on land.
