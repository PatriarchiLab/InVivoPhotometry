# InVivoPhotometry
Folder "doric_photometry_analysis" contains the functions for analysis of the recorded photometry signals. We used similar algorithms for calculating the dFF (see Methods), but structured these functions by the figures it relates to. The pre-processed data can be found in the Zenodo repository (see the data availability section). The data is sorted according to the experiment/figure and must be unzipped before analysis.
•	“average_open_field_data_AlloLiteCtrl” was used to calculate dFF signal for mice injected with AlloLite-ctr sensor and performing open filed test, average it across mice, and quantify the AUC (Figure5 i,j,k, Extended Fig6 h)

•	“average_open_field_data” was focused on calculating dFF signal for mice injected with dLight1.3b sensor and performing an open filed test, averaging it across mice and quantifying the AUC (Figure5 l,m,n, Extended Fig6 d)

•	“average_pavlovian_expert_days” was used to calculate dFF signal during the Pavlovian Conditioning experiment, average it across mice and quantify the peak response as well as Licking response parameters (Figure7 d, e, f, g)

•	“average_unpredicted_reward” was used to calculate dFF signal during the Unpredicted Reward experiment, calculate the percent of the detected response and average it across mice and quantify the peak response as well as Licking duration (Figure7 j,k,l,m,n)

•	” average_open_field_behavior” was used to quantify the behavior metrics of mice performing Open Field test (Extended Fig 5 f-i)

Folder “utilities” contains the functions used within photometry data analysis (for producing and saving the plots, calculating z-score etc)
