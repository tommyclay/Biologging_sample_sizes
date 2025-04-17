# Modest sample sizes are often sufficient for biologging to inform marine conservation

Thomas A. Clay, Richard A. Phillips Ana P.B. Carneiro, Janet R. D. Silk, Joan Navarro, Jacob González-Solís, Andrea Manica

## Overview

The following repository contains codes and datasets to 1) predict minimum sample sizes required based on spatial biologging data, and to reproduce 2) analyses and 3) figures from the paper, submitted to Conservation Biology. 

## Codes

The codes are laid out in three sections: 

### 1) Predicting minimum sample sizes required (apply to own spatial biologging data)

**1a_resampling_predicting_sample_sizes.r**: for example datasets, resample home-range area occupied with increasing sample size, run asymptotic models, choose the most appropriate and the predict the minimum sample size required. The resampling code can take a long time to run depending on the size of the dataset. 

**1b_individual_overlap_sample_sizes.r**: for example datasets, calculate the average among-individual spatial overlap (Bhattacharyya's affinity, BA) score and use to predict a coarse value of the mum sample size required. 


### 2) Reproducing analyses from the paper

**2a_resampling_area_occupied_predicting_sample_sizes.r**: loads in GLS datasets and resamples home-range area occupied with an increasing sample size for the 50% and 95% utilization distributions (UDs), runs asymptotic models on resampled outputs and the predicts minimum sample size required.

**2b_resampling_birdlife_predicting_sample_sizes.r**: loads in GLS datasets and calculates kernel UDs and then runs code from the track2KBA R package to estimate representativeness of datasets based on the inclusion approach. 

**2c_glmms_comparing_predicted_sample_sizes.r**: loads in minimum sample sizes based on area occupied and inclusion approaches (above) and runs GLMMs too compare predictions according to UD level (50% versus 95%), resampling approach, asymptotic model used, and the population size models are projected to.

**2d_testing_og_sample_size_on_predictions.r**: loads in GLS datasets and randomly subsets at intervals of n = 5 and then resamples the home-range area occupied. Asymptotic models are run on resampled outputs to predict minimum sample sizes required, and a GLMM run to examine the effect of original sample size on predictions. 

**2e_calculate_among_individual_overlap.r**: calculates spatial overlap among all pairwise combinations of individuals for each dataset. 

**2f_linking_individual_overlap_sample_sizes.r**: loads in minimum sample sizes and among-individual overlap scores based on GLS, GPS and PTT datasets from South Georgia and runs a model to establish a relationship between the two variables. Overlap scores for 52 species from Oppel et al. (2018) are loaded in and minimum sample size are generated based on this model. 

**2g_resampling_fisheries_overlap.r**: loads in GLS datasets for four species and iteratively calculates spatial overlap with pelagic longline fisheries (based on effort data from Global Fishing Watch) for an increasing number of individuals. 


### 3) Reproducing figures from the paper

**3a_fig2_asymptotic_plots.r**: loads in resampled datasets and models to predict the asymptotic relationships between sample size and home-range area (Fig. 2).

**3b_fig3_sensitivity_analysis_og_sample_size.r**: loads in predictions from the sensitivity analysis (2d, above) to plot the relationship between original sample sizes and minimum sample sizes required 
(Fig. 3). 

**3c_fig4_individual_overlap_maps.r**: loads in GLS datasets, creates kernel UDs for each individual and the whole population, and plots maps for each species and year combination (Fig. 4).

**3d_fig5_among-individual_overlap_sample_size.r**: loads in minimum sample sizes and among-individual overlap scores from South Georgia datasets and plots the relationship (Fig. 5a). Also loads in Oppel et al. (2018) values and plots BA overlap scores and minimum sample sizes at the Family level (Fig. 5b, c). 

**3e_fig6_sample_size_fisheries_overlap.r**: loads in resampled fisheries overlap datasets and plots the relationships between sample size and the total value of and change in fisheries-overlap score (Fig. 6). 


## Data use statement (please read before using)

The attached archived file(s) contain data derived from long-term tracking projects on albatrosses and petrels at Bird Island (South Georgia). This is a request to please let us know if you use them.  Several people have made extensive efforts overseeing the data collection for the last 25+ years. If you plan to analyse the data, there are a number of reasons why it would be very helpful if you could contact Richard Phillips (raphil@bas.ac.uk) before doing so: 

1. Occasionally we discover and correct errors in the data.
2. The data are complex and workers who do not know the study system are likely to benefit from advice on interpretation.
3. At any one time, several people within the existing project collaboration are analysing data from this project. Someone else may already be conducting the analysis you have in mind and it is desirable both to prevent duplication of effort and the risk of compromising the work of a PhD student, for example, whose future career may depend on publications arising from analyses of these data.
4. In order to maintain funding for the project and for further analyses, every few years we submit proposals to funding agencies. It is therefore very helpful for those running the project to know which data analyses are in progress elsewhere.

If you are interested in analysing the detailed project data in any depth, you may find it helpful to have access to the full database rather than the files available here.  If so, we are always open to further collaboration.
