This folder contains all the main script used to reproduce simulations and analyses.

The folder "functions and secondary scripts" contains some functions that are called in the main scripts or alternative scripts that have been used while building up the project, for exploration.

# Codes

### Initial conditions

Scripts 1 & 1bis defined initial conditions for theoretical and empirical communities, respectively, to allows the reproducibility of the simulations.

### Simulations under Julia

Scripts 2.1.0 & 2.1.1 are Julia scripts that were used to perform the coevolution simulations, using theoretical initial conditions. The script ran on a HPC plateform.

Script 2.2 are Julia scripts that were used to perform the coevolution simulations, using empirical initial conditions. The script ran on a HPC plateform.

### Extracts indices from theoretical and empirical simulations 

Script 3.1 extracts community measures from theoretical simulations. The middle part of script 3.1, was also ran on a HPC plateform to gain time.

Script 3.2 extracts community measures from empirical simulations. The middle part of script 3.2, was also ran on a HPC plateform to gain time.

### Figures for the theoretical simulations 

Script 4 plots the results, from figure 2 to figure 4. Figure 1 was done in a sepate scripts labelled *Figure1.r*

### Analyse empirical data and simulations with empirical initial conditions (=with empirical network dimensions) 

Scripts 5 extracts communitiy measures from empirical networks. 

Script 6 merges community measures from empirical networks and simulations to do figure 5.

# Outputs of simulations and empirical data

Codes, outputs of simulations and empirical data can be found here: https://doi.org/10.5281/zenodo.13270452

These outputs of simulations and empirical data allow to start analyses from anywhere in the worklow. For example, if one wants to reproduce the figures and analyses, they can start analyses from script 4, without running the simulations again, because all intermediate outputs are included in the data folder.

### Project structure
```
.

└── data_zenodo

    ├── data
    
              ├── empirical
              
                     ├── initial_conditions_simulations
                     
                     ├── interactions
                     
                     └── outputs_simulations
                     
              ├── simulated
              
                     ├── initial_conditions_simulations
                     
                     └── outputs_simulations
                     
    ├── Figures
    
    ├── scripts
    
              ├── functions and secondary scripts
              
              | numbered scripts
              
    ├── README.md
```
The folder "data_zenodo/data/empirical/interactions" contains the files with empirical mutualistic interactions and species phenologies. *flow_pheno_empirical.csv* and *poll_pheno_empirical.csv* contain the empirical phenological parameters for plant and pollinator species, respectively: the mean activity day (mu) and its standard deviation (sde) representing the duration of the activity period. *matrices_empirical_networks.RData* contains an R object with the 17 networks used. Plants are in rows and pollinators in columns, with each cell representing the average interaction value across sampling rounds, corrected by abundances. You can access it in R via:

```R
#load data
load("matrices_empirical_networks.RData")

#see the structure (a list of 17 networks)
str(networks)

#access the first network
networks[[1]]
```
