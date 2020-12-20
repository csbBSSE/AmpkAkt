# AmpkAkt

This repository holds all the codes used for simulation and analysis of the AMPK-AKT network in the manuscript.

The data for figure 1 can be generated using the Random_parameter_sims.m script. Further figures consist of analysis done on bistable parameter sets, which can be found in the new folder generated upon running the Random_parameter_sims.m script.

The nullcline analysis is done by editing the Nullcline_script.m file in the Nullcline analysis folder (edit file_name variable to the path to the bistable parameter set file). The corresponding results are saved in the parameter set folder.

Sensitivity analysis is done similarly as well.

Bifurcation analysis is done by running the "AMPK_AKT_BD.m" script with similar edit to the file_name variable as before.

Clinical data analysis is carried out in R and the corresponding script is given.
