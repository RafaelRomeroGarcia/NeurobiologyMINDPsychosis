# NeurobiologyMINDPsychosis

This repository contains code and data created in support of the project *"TITLE"*. All code was written in Matlab, R, and Python. Folders, files, and first steps are described below.

## **Data**

The `Data` folder contains all the data required for running the analyses. Here are the files that need to be downloaded and stored in a specific location. The remaining files will be automatically generated:

-	The code used to compute MIND networks (`MIND_networks_PAFIP` folder) is available at https://github.com/isebenius/MIND and corresponds to [MIND_01_MIND.py](Code/MIND_01_MIND.py)
  
-	The `all_microsc_DesikanKilliany68.csv` file is available at https://github.com/netneurolab/netneurotools

-	The code used for the PCA-CCA analyses (`cca_pls_toolkit_dev-grcca` folder) is available at https://github.com/RafaelRomeroGarcia/cca_pls_toolkit

## **First steps**

1.	Download `Code` folder, which contains the scripts and functions used for the analyses.

2.	Download `Data` folder, which contains data used to run the analyses, and will include data derived from them.

3.	Create `volumes`, `centiles`, `degree`, `edges`, and `epicenters` folders in `Data`.


## **Code**

The `Code` folder contains all the code required for running the analyses and generate data and figures. All scripts are designed to be sequentially executed. Don't forget to change the location variable regularly. 

-	[MIND_01_MIND.py](Code/MIND_01_MIND.py) – computes MIND networks from FreeSurfer (by default stored in the surf/ folder). It returns a .csv file for each individual.

-	[MIND_02_degrees_and_edges.m](Code/MIND_02_degrees_and_edges.m) – calculates for each HC and FEP individual the edges, degrees, and their effect sizes from MIND networks (`MIND_networks_PAFIP` folder). The results are stored in `degree` and `edges` folders. It returns a .csv file for each type (degree, effsizes_degree, edges, effsizes_edges) and for each clinical outcome (cognition, BPRS, SAPS, SANS) in the effect sizes cases.

-	[MIND_03_brain_mapping.R](Code/MIND_03_brain_mapping.R) – generates the regional brain maps of: (1) MIND degree of HC and FEP, and (2) effect sizes of MIND degree from the .csv files previously generated.

-	[MIND_04_brain_mapping_3D.m](Code/MIND_04_brain_mapping_3D.m) – generates the 3D spatial representation of the effect sizes of MIND edges and degrees from the .csv files previously generated.

-	[MIND_05_maturational_features.m](Code/MIND_05_maturational_features.m) – computes the associations between (1) MIND and centiles, (2) MIND and psychosis co-vulnerability, (3) hierarchical neurodevelopment and MIND and centiles, (4) peaks of volume/velocity and MIND and centiles.

-	[MIND_06_epicenters.m](Code/MIND_06_epicenters.m) – computes the epicenters of the disease. It returns a .csv file for each dignosis (HC/FEP) and clinical outcome.

-	[MIND_07_epicenter_mapping.R](Code/MIND_07_epicenter_mapping.R) – generates the regional brain maps of HC and FEP epicenters from the .csv files previously generated.

-	[MIND_08_parfor_FEP_degrees.m](Code/MIND_08_parfor_FEP_degrees.m) – runs PCA-CCA analyses. It returns the neurobiological loadings of the HC network and the effect sizes of degree for each clinical outcome. See [MIND_09_CCA_cent_var_FEP_degrees.m](Code/MIND_09_CCA_cent_var_FEP_degrees.m) and [MIND_10_run_CCA.m](Code/MIND_10_run_CCA.m) scripts before execution.
 
-	[MIND_09_CCA_cent_var_FEP_degrees.m](Code/MIND_09_CCA_cent_var_FEP_degrees.m) – it is not a script to be executed, but for setting parameters by default for [MIND_08_parfor_FEP_degrees.m](Code/MIND_08_parfor_FEP_degrees.m) script. Default parameters include (see https://mlnl.github.io/cca_pls_toolkit/cfg/ for detailed information):

```matlab
% Machine settings

 cfg.machine.name = 'cca';
 
 cfg.machine.metric = {'correl' 'trexvarx' 'trexvary'}; 
 
 cfg.machine.param.name = {'VARx', 'VARy'}; % explained variance by the PCA components

 cfg.machine.param.VARx = 0.6:0.1:0.9; % variance of data kept in the principal components during the SVD step of PCA-CCA  
 
 cfg.machine.param.VARy = 1;   
 
 cfg.machine.svd.varx = 1; % variance of X kept during the SVD step of PCA-CCA 

 cfg.machine.svd.vary = 1; % variance of Y kept during the SVD step of PCA-CCA

 cfg.machine.alignw = 'wX';

% Framework settings

cfg.frwork.name = 'permutation';     
 
 cfg.frwork.split.nout % number of outer splits/folds
 
 cfg.frwork.nlevel = 1;
    
% Deflation settings

 cfg.defl.name = 'generalized'; 
    
% Environment settings

 cfg.env.comp = 'local'; %  ['local', 'cluster']

 cfg.env.save.tableHeading = {'set' 'varx' 'correl' 'pval' 'npcax'};
    
% Number of permutations

 cfg.stat.nperm = 1000;

 cfg.stat.nboot = 1000;

```

-	[MIND_10_run_CCA.m](Code/MIND_10_run_CCA.m) – it is not a script to be executed, but for setting parameters by default for [MIND_08_parfor_FEP_degrees.m](Code/MIND_08_parfor_FEP_degrees.m) script. Default parameters include (see https://mlnl.github.io/cca_pls_toolkit/res/ for detailed information):

```matlab
% Weight or loading visualization

res.gen.weight.type = 'correlation'; % for plotting Loadings instead of weights
%     res.behav.weight.norm = 'minmax'; % Normalize weight {'std','zscore'}
res.gen.weight.subset.type = 'significant'; % Filter weights by keeping only significant ones 
% {'top' (by keeping top ones based on absolute value) 'minmax' (by keeping top positive and negative ones)}
res.gen.weight.stat.ci = 'sd'; % {'pi'} Confidence interval
res.behav.weight.errorbar.side = 'one';

```

