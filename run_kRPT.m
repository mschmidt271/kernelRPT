% This script is the driver for kRPT_fxn.m

clear variables

% length of 1-D domain
omega = 1;
% reaction rate coefficient, as derived by law of mass action
k = 5;
% molecular diffusion coefficient
D = 1e-5;
% number of particles of respective species
NA = 1000;
NB = NA;
% initial concentration of reactants (assumed to be the same for both species)
C0 = 1;
% particle kernel variance (set to 0 for delta kernels)
Avar = 0.12;
Bvar = Avar;
% total run time for individual simulation
maxtime = 1e2;
% number of loops in ensemble of simulations
loops = 10;

[totmass, runtimes]  =  kRPT_fxn(omega, k, D, NA, NB, C0, Avar, Bvar, maxtime, loops);