# Using Bayesian Hierarchical Modeling to Track Distributed Industrial Equipment Performance

This GitHub repository has the code and files required to replicate and analyze the model implemented in the article.

It includes four scripts in R Markdown files (.Rmd), ordered from 1 to 4, ready to be run, a Utils R file with accessory functions and a data folder with three data/results files stored as RData. All Rmd files have their corresponding md files rendered for ease of readibility. 

To run the code, an extra RData file including the model posterior sampling is required. This model file is stored as a [bugs](https://cran.r-project.org/web/packages/R2WinBUGS/index.html) object (**model.RData**) and kept in the [Zenodo](https://doi.org/10.5281/zenodo.10413762) repository [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10413762.svg)](https://doi.org/10.5281/zenodo.10413762).

## 1. Code files
Code files are given a number according to the order they should be run. Scripts 02 and 04 can be run without script 01 and 02 being previously executed. Intermediate results are already stored in the data folder to speed execution.

### 01ModelExecution.Rmd
This file has the code required to sample the posterior of the model parameters using the [JAGS](https://mcmc-jags.sourceforge.io/) Gibbs-Sampler. To run the script, one needs to have [JAGS](https://mcmc-jags.sourceforge.io/) installed and the [R2JAGS](https://cran.r-project.org/web/packages/R2jags/index.html) library available to link R to JAGS. Due to confidentiality reasons the model is built using only information from four machines out of the 357 available. This lets the user see the model execute. However, the obtained results are not relevant due to data incompleteness. The actual article model with the whole installed printer base is stored in the **model.RData** file. 

### 02ModelResults.Rmd
Running this script lets the user analyze the model convergence, its posterior and posterior predictive distributions as well as generate the graphs shown in the article (article figures 1, 2 and 3) based on the complete model. The script requires  the complete BUGS model stored in the **model.RData** file in Zenodo.

### 03ModelPerformance.Rmd
This script runs the simulation needed to assess the model performance and validity as explained in section 5 in the article. It generates the data-frame needed to plot the OC curves shown in the paper (figure 4). The script takes a while to run: 3 minutes in a Intel Core i9, 64GB Windows 11 workstation. It runs 25 times each of the 9 scenarios for the 357 printers included, sampling each of them 5,000 times. One can run the next script (04OCcurves.Rmd) without the need to run this one, 03ModelPerformance, as the simulation results data is stored in the data folder, file **results.RData**.

### 04OCcurves.Rmd
This file has the code needed to generate the model Operating Characteristic (OC) curves shown in the article's figure 4. It needs to load the results.RData file that contains the simulation results used int the article

### 00_Utils.R
This script includes accessory functions to draws control charts.

## 2. Data files
Data files are stored in the data folder. They are three:

### data.RData

Includes the data-frame with the printer-week information for printers 7, 116, 207 and 216. It has 96 rows (printers 7, 116, 207 and 216 have 25, 24, 25 and 22 recorded weeks each) and 8 columns:

 1. **rowId**, the consecutive row number from the original 8,113 printer-weeks dataset for all 357 printers.
 2. **id**, printer ID
 3. **jagsId**, the ID needed for the JAGS model to run (01ModelExecution.Rmd) as it expects consecutively numbered
    printers. It goes from 1 to 4.
 4. **week**, observed week ISO number.
 5. **nTotal**, total number of errors logged by the printer on this week.
 6. **timeImHour**, total weekly printing time.
 7. **PVCsecPerc**, total percent time spent printing on PVC on this week.
 8. **jobsPerTime**, average number of jobs per printing hour.

### printWeeks.RData

Two column data-frame listing the observed weeks for each of the 357 printers analyzed. It has 8,113 rows (one per printer-week) and two columns:

 1. **id**, printer id.  
 2. **week**, observed week ISO number.

### results.RData

This data-frame records the aggregated results of each of the different model performance simulations as obtained from script 03ModelPerformance.Rmd. It has 225 rows (9 scenarios times 25 different delta values each) and 8 columns:

 1. **scenario**, scenario number, from 1 to 9 as three different usage conditions are combined with three different print times.
 2. **delta**, out of control multiplier used in the current simulation.
 3. **tp**, true positives observed in the simulation.
 4. **fn**, false negatives observed in the simulation. The sum of true positives and false negatives gives the total 1,785,000
    printer-weeks sampled per simulation (357 printers, 5,000 times each).
 5. **type**, scenario usage conditions. It is a combination of PVC printed and jobs per hour done. It has three levels: low (0% PVC and 1 job every 2 hours), medium (50% PVC and 1 job change per hour) and high (100% PVC and 3 job changes per hour).
 6.   **timeImHour**, total weekly printing time.
 7. **PVCsecPerc**, total percent time spent printing on PVC on this week.
 8. **jobsPerTime**, average number of jobs per printing hour.

## 3. Model files

Stored in [Zenodo](https://doi.org/10.5281/zenodo.10413762)

### model.RData
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10413762.svg)](https://doi.org/10.5281/zenodo.10413762)

It includes a [bugs](https://cran.r-project.org/web/packages/R2WinBUGS/index.html) model object with the results of sampling the proposed model with the total available 8,113 printer weeks from the 357 printers using script 01ModelExecution.Rmd. It contains 10,000 samples from 3 chains of the posterior distribution most relevant model parameters with a previous 1,000 burnin phase.
The model is stored in [Zenodo](https://doi.org/10.5281/zenodo.10413762) as its size, 350MB, and the absence of code versioning tracking makes it unsuitable for GitHub storage. 
