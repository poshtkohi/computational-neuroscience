maxNumCompThreads(1);
%parpool(2);
%delete(gcp('nocreate'));
%----- Prepare the Environment -------%
close all; clc; clear; clearvars;
format long g;
c = Constants; % Gets constants for the model
%------Read data from file -------%
CAi = 45e-9;
NAi = 8e-3;
opt.NAi = NAi;
opt.CAi = CAi;

%%g = 0.614215234326583;
%%I_PMCA_BAR = 0.740341573296545;
%%I_NCX_BAR = 1.97477417445101;

g_na = 4.51528407234581e-05;
g_ca = 8.51528407234581e-05;
I_PMCA_BAR = 3.67847491339453e-06;
I_NCX_BAR = 1.03163149888025;

%i_cai_tot = -(NCX_Current_Fit(CAi, I_NCX_BAR) + PMCA_Current_Fit(CAi, I_PMCA_BAR) + Leak_Current_Fit(CAi, g))

%return

%y0(1) = 0.4; % g s−1 % NA
y0(1) = 0.6; % g s−1 % CA
y0(2) = 0.9;   %  I_PMCA_Bar A/m^2
y0(3) = 1; % I_NCX_BAR A/m^2 

ops.LBounds = 0.0;
ops.UBounds = Inf;
ops.PopSize = 72;
ops.StopFitness = 1e-60;
ops.DiagonalOnly = 1;
ops.MaxIter = 1000000;    % maximum number of iterations needed for estimation of parameters
ops.DispModulo = 1; % disp messages after every i-th iteration
ops.Resume = 'no';
ops.SaveVariables = 'on';
ops.StopOnStagnation = 'off';
ops.TolUpX = '1e30*max(insigma)';
ops.TolFun = '1e-20 % stop if fun-changes smaller TolFun';
ops.tolhistfun = '1e-20';
ops.SaveFilename = 'variablescmaes-p2x-channel-parameters.mat';
ops.LogFilenamePrefix = 'outcmaes-p2x-channel-parameters  % files for output data';

%w = i_shape;

%sigma(1) = 1e-2;
%sigma(2) = 1e-6;
%sigma(3) = 0.5;

sigma = 0.5;

[X, F, COUNT, STOP, OUT, BESTEVER]= cmaes('fhngen_p2x_fit_channel_parameters', y0, sigma', ops, opt);
X
%--------------------%