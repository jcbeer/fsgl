# Fused Sparse Group Lasso Simulation Study

As reported in the manuscript *Incorporating Prior Information with Fused Sparse Group Lasso: Application to Prediction of Clinical Measures from Neuroimages*

## Workflow
### 1. Run simulation
Edit and run verions of the script s01_FSGLassoSimulation.R for each of the 9 simulation scenarios. Parallelize by splitting the 100 random seeds per scenario into 4 or more groups. This may take a few weeks to run. 

**Input:** None. \
**Output:** Files compagg.csv, partagg.csv, compdist.csv, spcompagg.csv, sppartagg.csv, spcompdist.csv, exsp.csv, missp.csv, misspsp.csv. Each row represents one simulation for the given alpha, gamma combination. 100 simulations were done for each alpha, gamma combination, so there are 2100 rows. \
Variables in these files include: \
**seed:** random seed \
**alpha:** alpha value \
**gamma:** gamma value \
**opt.lambda:** value of lambda yielding the minimum 5-fold cross-validation error \
**mean.cve:** the average mean squared error of predicted y across 5 cross-validation folds at the optimal lambda \
**mean.cve:** the standard deviation of mean squared error of predicted y across 5 cross-validation folds at the optimal lambda \
**mse.betas:** mean squared error of estimated betas \
**mse.y:** mean squared error of predicted y for the training set \
beta1 to beta400: estimated beta values \


### 2. Calculate test set error and bias-variance decomposition
Run the script s02_FSGLassoSimulation_TestDataPredictionError.R. This script will generate new test set data and calculate test set prediction error for each simulation, and will decompose the mean squared error into bias squared and variance. 

**Input:** Files compagg.csv, partagg.csv, compdist.csv, spcompagg.csv, sppartagg.csv, spcompdist.csv, exsp.csv, missp,csv, misspsp.csv. \
**Output:** Files compaggplot.csv, partaggplot.csv, compdistplot.csv, spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, exspplot.csv, misspplot.csv, misspspplot.csv, and compaggbv.csv, partaggbv.csv, compdistbv.csv, spcompaggbv.csv, sppartaggbv.csv, spcompdistplot.csv, exspplot.csv, misspplot.csv, misspspplot.csv. \
Variables in the \*plot.csv files include variables defined above and: \
**mse.y.test:** mean squared error of predicted y for the test set \
Variables in the \*bv.csv files include: \
**id:** observationID \
**X.beta:** the true outcome y values  \
**yhat.var\*:** variance of all 100 simulation predicted y's for that observation for the given alpha and gamma values \
**yhat.bias\*:** estimated bias of all 100 simulation predicted y's for that observation for the given alpha and gamma values \
**yhat.bias2\*:** estimated bias squared of all 100 simulation predicted y's for that observation for the given alpha and gamma values \
**yhat.mse\*:** estimated mean squared error of all 100 simulation predicted y's for that observation for the given alpha and gamma values \


### 3. Make tables
#### Table 2
Use the script FSGLassoSimulation_Table2.R. \
**Input:** Files compaggplot.csv, partaggplot.csv, compdistplot.csv, spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, exspplot.csv, misspplot.csv, misspspplot.csv. \
**Output:** LaTeX code printed to console. \


#### Supplemental Tables S2 to S10
Use the script FSGLassoSimulation_SupplementalTable_S2toS10.R. \
**Input:** Files compaggplot.csv, partaggplot.csv, compdistplot.csv, spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, exspplot.csv, misspplot.csv, misspspplot.csv, and compaggbv.csv, partaggbv.csv, compdistbv.csv, spcompaggbv.csv, sppartaggbv.csv, spcompdistplot.csv, exspplot.csv, misspplot.csv, misspspplot.csv. \
**Output:** LaTeX code printed to console. \


### 3. Make figures
#### Figure 1. Group Structure and True Coefficients
Use the script FSGLassoSimulation_Figure1_GroupStructureFig.R. \
**Input:** None. \
**Output:** 1 pdf figure. \


#### Figures 2 to 4. Box Plots
Use the script FSGLassoSimulation_Figure1_GroupStructureFig.R. \
**Input:** Files compaggplot.csv, partaggplot.csv, compdistplot.csv, spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, exspplot.csv, misspplot.csv, misspspplot.csv. \
**Output:** 3 pdf figues. \


#### Supplementary Figure S1.
Use the script FSGLassoSimulation_FigureS1_BiasVariancePlots.R. \
**Input:** Files compaggbv.csv, partaggbv.csv, compdistbv.csv, spcompaggbv.csv, sppartaggbv.csv, spcompdistplot.csv, exspplot.csv, misspplot.csv, misspspplot.csv. \
**Output:** 1 pdf figues. \