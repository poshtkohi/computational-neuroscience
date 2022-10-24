%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

maxNumCompThreads(1);
%parpool(2);
%delete(gcp('nocreate'));
%----- Prepare the Environment -------%
close all; clc; clear; clearvars;
format long g;
c = Constants; % Gets constants for the model
%------Read data from file -------%
CAx = 2e-3;   % M
CAi = 45e-9;  % M
NAi = 8e-3;   %M
NAx = 130e-3; % M

%%g = 0.614215234326583;
%%I_PMCA_BAR = 0.740341573296545;
%%I_NCX_BAR = 1.97477417445101;

g = 4.51528407234581e-05;
I_PMCA_BAR = 3.67847491339453e-06;
I_NCX_BAR = 1.03163149888025;


R = 8.314; % Ideal gas constant; unit: J/K.mol;
T = 310; % Absolute temperature; unit: K;
F = 96485.33212; % Faraday's constant; unit: C/mol;
Vm = -60e-3;
Z_Ca = 2;
Z_Na = 1;
Ena =  R * T * Z_Na^-1 * F^-1 * log(NAx/NAi);
Eca =  R * T * Z_Ca^-1 * F^-1 * log(CAx/CAi);
%g = 0.05;
%g = 0.013954; % s−1
%%g = 0.614215234326583;
g = 0.253465502377278;
%I_CA_Leak = g * (Vm - E);
%I_Ca_Leak = 0;

I_NA_NCX = NA_NCX_Current(NAi, CAi, Vm)
I_CA_NCX = -2 * 3^-1 * I_NA_NCX;
%I_NA_NCX = 3 * NCX_Current(NAi, CAi);
%I_CA_NCX = -2 * NCX_Current(NAi, CAi);

%I_CA_NCX*317.83*1e-12 * 1e12
%I_CA_NCX
%return
pmca = PMCA_Current(CAi)

i_na_efflux = 0;%300*NAi;
g_NA_Leak = (-I_NA_NCX) / (Vm - Ena)
%I_CA_NCX = 0;
g_CA_Leak = (-I_CA_NCX - PMCA_Current(CAi)) /  (Vm - Eca)
%--------------------%