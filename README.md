# DRO-FDED
Thanks for your interest in this project. 
This project includes the following files:
- Simulink file for the primary frequency response model: PFR_IEEE_39_two_DSR.mdl. 
- Monte_Carlo_Simulation.m.  This file calls the PFR_IEEE_39_two_DSR.mdl and launches the linear regression and cross-validation. The fitted model is stored in: Coefficient.mat, residue_variance.mat, and max_regression_residue.mat.
- IEEE39_Parameter.m. This file stores the data for the IEEE 39-bus system.
- DRO_OPF.m. This file solves the distributional robust economic dispatch model.
- Robust_OPF.m. This file solves the robust and deterministic economic dispatch model.
- DRO_reformulation.pdf. This file explains the reformulation process for the DRO problem.
