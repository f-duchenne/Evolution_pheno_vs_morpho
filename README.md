This folder contains all the main script used to reproduce simulations and analyses.

# Codes

Scripts 1 & 1bis define initial conditions for theoretical and empirical communities, respectively, to allows the reproducibility of the simulations.

Script 2 is a Julia script, which is used to perform the coevolution simulations, and ran on a HPC plateform.

Script 3 extracts communitiy measures from *in silico* simulations. The middle part of script 3, was also ran on a HPC plateform to gain time.

Script 4 plots the results, from figure 2 to figure 4. Figure 1 was done in a sepate scripts labelled *Figure1.r*

Scripts 5 extracts communitiy measures from empirical networks and from simulations performed with initial empirical network sizes. 

Script 6 merges community measures from empirical networks and simulations to do figure 5.

The folder "functions and secondary scripts" contains some functions that are called in the main scripts or alternative scripts that have been used while building up the project, for exploration.

# Empirical data

The empirical data can be found here: *link to come*

In this folder *flow_pheno_empirical.csv* and *poll_pheno_empirical.csv* contain the empirical phenological parameters for plant and pollinator species, respectively: the mean activity day (mu) and its standard deviation (sde) representing the duration of the activity period.

*matrices_empirical_networks.RData* contains an R object with the 17 networks used. Plants are in rows and pollinators in columns, with each cell representing the average interaction value across sampling rounds, corrected by abundances.

You can access it in R via:

```R
#load data
load("matrices_empirical_networks.RData")

#see the structure (a list of 17 networks)
str(networks)

#access the first network
networks[[1]]
```
