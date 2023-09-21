# living-up-to-thermal-potentials

## Introduction
This repository contains the data and code used to produce results presented in Moore et al. 2023: _Temperate species underfill their tropical thermal potentials on land_. This study tests hypotheses about how the role of temperature in limiting species ranges might be overriden by antagonistic species interactions and other range-limiting factors. The study was designed to be reproducible, and running all scripts in order should successfully reproduce analyses and figures.

![Pretty picture](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/figures/main/inkscape-files/fig2_niches_final.png)

## Data 
A minimum dataset needed to reproduce all main analyses and large files that cannot be hosted on Github can be downloaded from the paper's Figshare repository, which is linked in the article Data Availability statement. For scripts to successfully run, the large files should be donwloaded into a folder titled 'large-files' in the repository's working directory.
  
Additional temperature, elevation, and depth data files used in analyses can be downloaded here:
 - [Berkeley Earth Daily Land (Experimental; 1880 – Recent)](http://berkeleyearth.org/data/) - Average High Temperature (TMAX) and Average High Temperature (TMAX)]
 - [NOAA OI SST V2 High Resolution Dataset](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html) - Mean Daily Sea Surface Temperature
 - [General Bathymetric Chart of the Oceans](https://gebco.net/data_and_products/gridded_bathymetry_data/) - TID Grid
 - [EarthEnv Global 10-km topography](http://www.earthenv.org/topography) - Maximum and Minimum, GMTED2010
 
Realized range polygons needed to run analysis can be downloaded here:
 - [IUCN Spaial Files](https://www.iucnredlist.org/resources/spatial-data-download) - AMPHIBIANS, BLENNIES, CLUPEIFORMES, FILTERED, FW_FISH, MAMMALS, MARINEFISH, PUFFERFISH, REPTILES, SEABREAMS_PORGIES_PICARELS, SHARKS_RAYS_CHIMAERAS, WRASSES_PARROTFISHES
 - [GARD Version 1.1](http://www.gardinitiative.org/data.html)

## Code
The following are all of the R scripts contained in this repository and a short description of what each accomplishes:
#### [01_creating-temperature-depth-elev-data.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/01_creating-temperature-depth-elev-data.R)
 - prepares global air and sea surface temperature, depth, and elevation data for use in analyses 
#### [02_collating-realized-range-maps.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/02_collating-realized-range-maps.R)
 - extracts and organizes realized range polygons  
#### [03_trait-wrangling.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/03_trait-wrangling.R)
- collates and cleans trait data for species incldued in the analysis
#### [04_acclimatisation-analysis.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/04_acclimatisation-analysis.R)
- collates and cleans pre-existing acclimation response ratio databases for use in the acclimatisation analysis 
#### [05_modelling-operative-temperatures.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/05_modelling-operative-temperatures.R)
- models operative temperatures in different microhabitat for terrestrial species 
#### [06_creating-potential-ranges.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/06_creating-potential-ranges.R)
- creates potential thermal ranges
#### [07_quantifying-niche-filling.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/07_quantifying-niche-filling.R)
- measures potential thermal niche filling in thermal space
#### [08_quantifying-range-filling.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/08_quantifying-range-filling.R)
- measures potential thermal range filling in geographic space
#### [09_model-selection_niche.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/09_model-selection_niche.R)
- fits models to warm and cool niche filling and plots predictions
#### [10_model-selection_range.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/10_model-selection_range.R)
- fits models to range filling and asymmetry in underfilling and plots predictions
#### [11_sensitivity_model-selection_niche.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/11_sensitivity_model-selection_niche.R)
- checks sensitivity of warm and cool niche filling results to behaviour and acclimatisation
#### [12_sensitivity_model-selection_range.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/12_sensitivity_model-selection_range.R)
- checks sensitivity of range filling and asymmetry in underfilling results to acclimatisation
#### [13_sensitivity_model-selection_niche_no-dormancy.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/13_sensitivity_model-selection_niche_no-dormancy.R)
- checks sensitivity of warma and cool niche filling results to exclusion of dormant species 
#### [14_sensitivity_model-selection_range_no-dormancy.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/14_sensitivity_model-selection_range_no-dormancy.R)
- checks sensitivity of range filling and asymmetry in underfilling results to exclusion of dormant species 
#### [15_sensitivity_range-source.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/15_sensitivity_range-source.R)
- checks sensitivity of warm and cool niche filling results to realized range source 
#### [16_niche-filling-figures.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/16_niche-filling-figures.R)
- creates warm and cool niche filling figures
#### [17_range-filling-figures.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/17_range-filling-figures.R)
- creates range filling figures
#### [18_sensitivity_NicheMapR.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/18_sensitivity_NicheMapR.R)
- tests sensitivity of results on land to parameters used in operative temperature models
#### [19_phylo-gls-sensitivity-analysis.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/19_phylo-gls-sensitivity-analysis.R)
- tests sensitivity of results to choice of method used to control for evolutionary relatedness
#### [20_create-minimum-dataset.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/20_create-minimum-dataset.R)
- creates a minimum dataset needed to reproduce main analyses
  
## Ownership
This reposity is owned by Nikki A. Moore and Jennifer M. Sunday. All data, scripts, and figures may be used by others with proper acknowledgement.
