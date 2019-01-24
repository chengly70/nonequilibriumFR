# nonequilibriumFR
Code to implement moment closure method on coupling firing rate (Wilson-Cowan) heterogeneous network receiving correlated 
background inputs.  All functions/scripts written in MATLAB

Driver Files:
scptAn_3.m -- driver file to simulate 3 cell network (Figs 1, 2), calls on rhs_momentsFP.m, getFstats.m, 
          iter_method.m, saves results in dAn_n3Tw.mat

scptAn_3sin.m -- almost the exact same as scptAn_3 but has sinusoidal mu_t input, saves results in dAn_n3Tw_sin.mat
scptAn_50.m -- driver file to simulate 40 cell network, dense coupling; allow varying coupling strength, 
          calls on rhs_momentsFP.m, getFstats.m, iter_method.m, saves results in dAn_n50[a/b/c/d].mat, a/b/c/d denote different coupling 
          strengths
scptAn_50sin.m -- almost the exact same as scptAn_50 but has sinusoidal mu_t input, saves results in dAn_n50_sin.mat

scptMC.m -- driver to run the Monte Carlo simulations AFTER the scptAn_[].m script(s) have been run; same file for pulse input 
or sinusoidal input.  Calls on mc_WCpstrtT.m and LOADS either dAn_n3Tw.mat or dAn_n3Tw_sin.mat; 
!!CHANGE file name save to MATCH which dAn_ file loading!!


Main RHS file
rhs_momentsFP.m -- RHS of our method; nonlinear ODEs solved with ode45
mc_WCpstrT.m -- function to run Monte Carlo simulation
getFstats.m -- returns firing rate stats F(x) given x and F parameters (assuming sigmoidal transfer function)
iter_method.m -- our old steady-state PRE method just to setup Init. Cond. quickly

DISPLAY files:
plots_net3.m -- shows Nc=3 networks, either fast pulse or sinusoidal input (change which_set variable in file)
plots_50.m -- shows Nc=50 networks with pulse input, 
           showing only the error (can easily show individual comparisons like Figs1-2, all data is there)
plots_50sin.m -- shows Nc=50 networks with sinusoidal input (similar to plots_50.m)

MAT files are called by plots_[].m, generated by the scripts/functions in this directory.
