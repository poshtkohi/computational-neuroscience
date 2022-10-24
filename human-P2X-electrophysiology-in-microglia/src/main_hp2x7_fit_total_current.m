%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

maxNumCompThreads(1);
%parpool(2);
%delete(gcp('nocreate'));
%----- Prepare the Environment -------%
close all; clc; clear; clearvars;
format long;
c = Constants; % Gets constants for the model
resume = c.resume; % yes or no
pop_size = c.PopSize;
%------Read data from file -------%
d = load('data/hP2X7-ATP-5mM-total-current.dat');

ATP = 5e-3; % unit in millimollars
ATP_t1 = 2; ATP_t2 = 4.6; t_max = 8.5;% unit in seconds
t = d(:, 1);
i_tot = d(:, 2);% * 1e-3; % pA
i_shape = normalize(-i_tot, 'range');% * 2.416 / 2.5;
%i_shape = normalize(-i_tot, 'range');
%normalize(-i_tot, 'range') * 2.416 / 2.61;
cai_shape = cumtrapz(t, i_shape);
%plot(t, i_shape, t, cai_shape);
%figure;
%plot(t, i_shape);

global nT_P2X T;
T = 0; % Period in seconds
nT_P2X = 0;

% Smoothing the data sets
mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
f = fit(t, i_shape, mode);
t1 = linspace(0.0, t(length(t)), 500);
t = t1;
i_shape = f(t1);
%int = integrate(f, t1, 0);
%figure;
%plot(t, int);
%return
%----- Constants -------%
%nt = length(t);     % number of time points to simulate
%tol = 1e-20;
odeopt = odeset;%('RelTol', tol, 'AbsTol', tol);
%----- Variables -------%
tspan = [0.0 t(length(t))];
tint = t;  % all time points
%S = readmatrix('STM-P2X.xlsx'); % Loads the stoichiometry matrix. small model
%[NROWS, NCOLS] = size(S);
%---- Preparing for cureve fitting ------%
K_init = load('data/initial-guess-k-hp2x7-total-current.dat');
%previous_experiment = '1000';
[K_new, num_added] = load_initial_conditions_from_file_hp2x7_total_current(K_init);
opt.K = K_new;
opt.initial_condition_index = length(K_init) + 1;
opt.x0 = zeros(4, 1);
opt.x0(1, 1) = 1;
opt.x0 = load_initial_conditions_from_k_hp2x7_total_current(opt.K, opt);
%opt.S = sparse(S);
opt.odeopt = odeopt;
opt.tspan = tspan;
opt.tint = tint;
opt.ATP = ATP;
opt.ATP_t1 = ATP_t1;
opt.ATP_t2 = ATP_t2;
opt.observable_index = 4; % Cytosolic concentration of calcium

%return;
global rt_fig fig_states_rt err_prev counter total_elapsed_time;
rt_fig = figure('Name', 'Real-time fitted current curve');
fig_states_rt = figure('Name', 'hP2X7 states before fitting');
fig_states_rt.WindowStyle = 'docked';
err_prev = Inf;
counter = 1;
total_elapsed_time = 0.0;
rt_fig.WindowState = 'maximized';
rt_fig.Units = 'normalized';
rt_fig.OuterPosition = [0 0 1 1];
rt_fig.WindowStyle = 'docked';
%----- Curve Fitting -------%
% Parameter estimation using GA (genetics algorithm) 
y0 = opt.K;          % Initial guess for k
%ops.LBounds = LB;      % Sets of parameters to be estimated
%opt.UBounds = UB; 
%for i=1:1:length(y0)
%    if i > length(opt.K) - 2
%        UB(i) = 1.0;
%    else
%        UB(i) = Inf;
%    end
%end
%ops.LBounds = 0.0;
%ops.UBounds = Inf;
ops.LBounds = zeros(length(y0), 1);
ops.UBounds = Inf * ones(length(y0), 1);

%{
ops.LBounds(25) = 0.1;%
ops.LBounds(26) = 0.1;%
ops.LBounds(27) = 0.1;%

ops.UBounds(25) = 1.1;%
ops.UBounds(26) = 1.1;%
ops.UBounds(27) = 1.1;%
%}
%{
ops.LBounds(7) = 1;%
ops.LBounds(10) = 1;%
ops.LBounds(13) = 1;%
ops.LBounds(16) = 1;%
ops.LBounds(19) = 1;%
ops.LBounds(22) = 1;%

ops.UBounds(7) = 2;%
ops.UBounds(10) = 2;%
ops.UBounds(13) = 2;%
ops.UBounds(16) = 2;%
ops.UBounds(19) = 2;%
ops.UBounds(22) = 2;%
%}
ops.PopSize = pop_size;
opt.PopSize = ops.PopSize;
ops.StopFitness = 1e-16;
ops.DiagonalOnly = 1;
ops.MaxIter = 1000000;    % maximum number of iterations needed for estimation of parameters
ops.DispModulo = 1; % disp messages after every i-th iteration
ops.Resume = resume;
ops.SaveVariables = 'on';
ops.StopOnStagnation = 'off';
ops.TolUpX = '1e13*max(insigma)';
ops.TolFun = '1e-40 % stop if fun-changes smaller TolFun';
ops.tolhistfun = '1e-20';
ops.SaveFilename = sprintf('variablescmaes-hp2x7-ATP-%d-popsize-%d-total-current.mat', ATP, pop_size);
ops.LogFilenamePrefix = sprintf('outcmaes-hp2x7-ATP-%d-popsize-%d-total-current  % files for output data', ATP, pop_size); 

w = i_shape;

sigma = 0.5 * ones(length(y0), 1);
%sigma(7:24) = 0.1;

[X, F, COUNT, STOP, OUT, BESTEVER]= cmaes('fhngen_hp2x7_fit_total_current', y0, sigma, ops, opt, w);

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