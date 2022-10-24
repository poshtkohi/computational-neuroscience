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
resume = 'no'; % yes or no
pop_size = 10;
%------Read data from file -------%
d = load('data/rNCX.dat');

v = d(:, 1) * 1e-3;
i_tot = d(:, 2) * 1e-12; % pA
%cai_shape = cumtrapz(t, i_shape);
%plot(t, i_shape, t, cai_shape);
%%figure;
%%plot(v, i_tot);

%----- Variables -------%
%---- Preparing for cureve fitting ------%
k(1) = 1; % I_NCX_BAR A/m2
k(2) = 0.5; % gama NCX partition parameter
k(3) = 3; % n
k(4) = 8e-3;  % NAi M
k(5) = 45e-9;   % CAi M
opt.v = v;

%return;
global rt_fig;
rt_fig = figure('Name', 'rNCX fitted curve');
rt_fig.WindowState = 'maximized';
rt_fig.Units = 'normalized';
rt_fig.OuterPosition = [0 0 1 1];
rt_fig.WindowStyle = 'docked';
%----- Curve Fitting -------%
% Parameter estimation using GA (genetics algorithm) 
y0 = k;          % Initial guess for k
%ops.LBounds = LB;      % Sets of parameters to be estimated
%opt.UBounds = UB; 
%for i=1:1:length(y0)
%    if i > length(opt.K) - 2
%        UB(i) = 1.0;
%    else
%        UB(i) = Inf;
%    end
%end
ops.LBounds = 0.0;
ops.UBounds = Inf;
ops.PopSize = pop_size;
opt.PopSize = ops.PopSize;
ops.StopFitness = 1e-100;
ops.DiagonalOnly = 1;
ops.MaxIter = 1000000;    % maximum number of iterations needed for estimation of parameters
ops.DispModulo = 1; % disp messages after every i-th iteration
ops.Resume = resume;
ops.SaveVariables = 'on';
ops.StopOnStagnation = 'off';
ops.TolUpX = '1e9*max(insigma)';
ops.TolFun = '1e-30 % stop if fun-changes smaller TolFun';
ops.tolhistfun = '1e-30';
ops.SaveFilename = 'variablescmaes-rncx.mat';
ops.LogFilenamePrefix = 'outcmaes-rncx % files for output data';

w = i_tot;

sigma = 0.5;%[1.0 0.5 0.5 1e-3]';

CAi = 45e-9;  % M
NAi = 8e-3;    % M
for i=1:1:length(v)
    I_NA_NCX = NA_NCX_Current_Fit(k, v(i));
    I_CA_NCX = -2 * 3^-1 * I_NA_NCX;
    s_p = 243.53 * 1e-12; % the area of processes
    s_b = 74.3 * 1e-12; % the area of cell body
    s_tot = s_p + s_b;
    ww(i) = (I_NA_NCX + I_CA_NCX) * s_tot;
end
figure;
plot(v, ww * 1e12);

I_NA_NCX = NA_NCX_Current_Fit(k, v(i));
I_CA_NCX = -2 * 3^-1 * I_NA_NCX;
s_p = 243.53 * 1e-12; % the area of processes
s_b = 74.3 * 1e-12; % the area of cell body
s_tot = s_p + s_b;
I_Bar = 17.26e-12 / (I_CA_NCX * s_tot)
return;

[X, F, COUNT, STOP, OUT, BESTEVER]= cmaes('fhngen_rncx_fit', y0, sigma, ops, opt, w);

printk(X);
%------Functions -------%
function printk(k)
    m = length(k);
    s = "";
    for i=1:1:m
        s = sprintf("%s %g", s, k(i));
    end
    disp(s);
end
%--------------------%