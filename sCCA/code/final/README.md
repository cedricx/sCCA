### go1sample_split 
* Use this file to apply exclusion criteria and create the discovery and replication samples 

### compile_network_pwr 
* Use this file to grab connectivity matrices and clinical data from CfN cluster

### Select features_pwr 
* Apply median absolute deviation to select the top 10% edges
* Regress out covariates (age, sex, race, motion)

### sCCA_pwr 

* Grid Search for regularization parameter
* Covariance explained
* CCA correlation
* Permutation testing
* Visulization of clinical items
* Bootstrapping (resampling)
* Calculate within/between module connectivity

### ccafunctions 
* custom functions that sCCA_pwr calls


