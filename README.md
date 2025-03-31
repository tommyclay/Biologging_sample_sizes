# Modest sample sizes are often sufficient for biologging to inform conservation

Thomas A. Clay, Richard A. Phillips Ana P.B. Carneiro, Janet R. D. Silk, Joan Navarro, Jacob González-Solís, Andrea Manica

## Overview and data

The following repository contains codes to 1) predict minimum sample sizes required based on spatial biologging data, and to reproduce 2) analyses and 3) figures from the paper, submitted to Conservation Biology. 

## Codes

The codes are laid out into three sections: 

### 1) Predicting minimum sample sizes required based on spatial biologging data

**1a_resampling_predicting_sample_sizes.r**: This code...

**1b_individual_overlap_sample_sizes.r**: This comprises code....


### 2) Reproducing analyses from the paper

**2a_resampling_area_occupied_predicting_sample_sizes.r**: This code to...

**2b_resampling_birdlife_predicting_sample_sizes.r**: This comprises code to...

**2c_glmms_comparing_predicted_sample_sizes.r**: This comprises code to...

**2d_testing_og_sample_size_on_predictions.r**: This comprises code to...

**2e_calculate_among_individual_overlap.r**: This comprises code to...

**2f_linking_individual_overlap_sample_sizes.r**: This comprises code to...

**2g_resampling_fisheries_overlap.r**: This comprises code to...


### 3) Reproducing figures from the paper

**3a_fig2_asymptotic_plots.r**: load in resampled datasets and models to predict asymptotic relationships between the number of individuals tracked and home-range area (Fig. 2).

**3b_fig3_sensitivity_analysis_og_sample_size.r**: load in predictions from sensitivity analysis (2d, above) to plot the relationship between original sample sizes and minimum sample sizes required 
(Fig. 3). 

**3c_fig4_individual_overlap_maps.r**: load in GLS datasets, creates kernel UDs for each individual and the population and plots maps for each species and year combination (Fig. 4).

**3d_fig5_among-individual_overlap_sample_size.r**: This comprises code...

**3e_fig6_sample_size_fisheries_overlap.r**: load in resampled fisheries overlap datasets and plots the relationships between the number of individuals tracked and the total and change in fisheries overlap scores (Fig. 6). 


## Data use statement (please read before using)

The attached archived file(s) contain data derived from long-term tracking projects on albatrosses and petrels at Bird Island (South Georgia). This is a request to please let us know if you use them.  Several people have made extensive efforts overseeing the data collection for the last 25+ years. If you plan to analyse the data, there are a number of reasons why it would be very helpful if you could contact Richard Phillips (raphil@bas.ac.uk) before doing so: 

1. Occasionally we discover and correct errors in the data.
2. The data are complex and workers who do not know the study system are likely to benefit from advice on interpretation.
3. At any one time, several people within the existing project collaboration are analysing data from this project. Someone else may already be conducting the analysis you have in mind and it is desirable both to prevent duplication of effort and the risk of compromising the work of a PhD student, for example, whose future career may depend on publications arising from analyses of these data.
4. In order to maintain funding for the project and for further analyses, every few years we submit proposals to funding agencies. It is therefore very helpful for those running the project to know which data analyses are in progress elsewhere.

If you are interested in analysing the detailed project data in any depth, you may find it helpful to have access to the full database rather than the files available here.  If so, we are always open to further collaboration.
