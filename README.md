# living-up-to-thermal-potentials

## Introduction
This repository contains the data and code used to assess how ectothermic animal species fill their potential thermal niches. This study was designed to be easily reproducible, and running all scripts in order should successfully reproduce all analyzes and figures.

Large files needed to re-run this analysis cannot be hosted on github and so must be downloaded from here:

For scripts to successfully run, large files should remain in a folder titled 'large-files' in the repository's working directory.

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
#### [05_modelling-operative-temperatures.R](https://github.com/nicole-a-moore/living-up-to-thermal-potentials/blob/main/R/05_modelling-operative-temperatures.R) - models operative temperatures in different microhabitat for terrestrial species 
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

## Ownership
This reposity is owned by Nikki Moore and Jennifer Sunday. All data, scripts, and figures may be used by others with proper acknowledgement.
