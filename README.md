### MOCHBreedingPlasticity
Data and code used for analyses and figures from Whitenack et al. 2025 "Individual repeatability and plasticity of   reproductive phenology in a resident montane bird", published in Behavioral Ecology and Sociobiology.

File preparation date: 17 July 2025

Authors: Lauren E. Whitenack, Benjamin R. Sonnenberg, Carrie L. Branch, Angela M. Pitera, Joseph F. Welklin, Virginia K. Heinen, Vladimir V. Pravosudov

Contact Lauren Whitenack at lauren.whitenack@gmail.com with questions

## List of R scripts 
# Script: plasanalyses.R
	Summary: Use this script to run all plasticity analyses and produce figures for publication. 

## List of data files

# Data file: mochhigh.csv
	Summary: 11 years of breeding data (2013-2023) from high elevation mountain chickadees at Sagehen Experimental Forest, CA, USA.
	Each row contains data from a single nest associated with one season.
	Column descriptions:
		YEAR - Field season year
		SEASON - Winter season preceding that year's breeding season
		BOX - Nest box ID
		ELEVATION - 'H' = high elevation; 'L' = low elevation
		NEST.TYPE - all nests are initial nests ('INITIAL'), no renests or second nesting attemps included
		FIRST.EGG - date of first egg laid; format = DD/MM/YYYY
		J.FIRST.EGG - julian date of first egg laid
		CLUTCH - number of eggs at initiation of incubation
		BROOD - number of eggs at day 16 in the nest
		MEANMASS - mean nestling mass (g)
		CV - coefficient of variation in nestling mass ((standard deviation÷mean)×100) 
		year - duplicate of YEAR, to match column names with snotel_data

# Data file: mochlow.csv
	Summary: 11 years of breeding data (2013-2023) from low elevation mountain chickadees at Sagehen Experimental Forest, CA, USA.
	Each row contains data from a single nest associated with one season.
	Column descriptions:
		YEAR - Field season year
		SEASON - Winter season preceding that year's breeding season
		BOX - Nest box ID
		ELEVATION - 'H' = high elevation; 'L' = low elevation
		NEST.TYPE - all nests are initial nests ('INITIAL'), no renests or second nesting attemps included
		FIRST.EGG - date of first egg laid; format = DD/MM/YYYY
		J.FIRST.EGG - julian date of first egg laid
		CLUTCH - number of eggs at initiation of incubation
		BROOD - number of eggs at day 16 in the nest
		MEANMASS - mean nestling mass (g)
		CV - coefficient of variation in nestling mass ((standard deviation÷mean)×100) 
		year - duplicate of YEAR, to match column names with snotel_data

# Data file: snotel_data.csv
	Summary: Daily weather data downloaded from USDA SNOTEL stations (https://wcc.sc.egov.usda.gov/nwcc/). High elevation 
	data come from SNOTEL Site 541 and low elevation data come from taking the average of SNOTEL Sites 539 and 540. Both the
	raw data for Sites 539 and 540 and the averaged data for low elevation are included in this data file.
	Column descriptions:
		Columns with prefix X539 - Data from SNOTEL Site 539
		Columns with prefix X540 - Data from SNOTEL Site 540
		Columns with prefix X541 - Data from SNOTEL Site 541; used for high elevation analyses
		Columns with prefix Low - Average of Site 539 and 540; used for low elevation analyses

