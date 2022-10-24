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
BzATP = 300e-6; % unit in milimollars
ATP = 11 * BzATP; % BzATP is 10/12-fold more potent than ATP for P2X7 receptors in human microglia
%ATP_t1 = 2; ATP_t2 = 4.6; t_max = 8.5;% unit in seconds - Orginal
ATP_t1 = 30.0/60; ATP_t2 = (116.0 + 30.0)/60; t_max = 180/60;% unit in seconds

d = load('data/hP2X7-BzATP-300uM-CAi.dat');
t = d(:, 1) / 60;
cai = d(:, 2); % AU
%cai = normalize(d(:, 2), 'range');
%figure; plot(t, cai);

% Smoothing the data sets
mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
f = fit(t, cai, mode);
t1 = linspace(0.0, t(length(t)), 500);
t = t1;
cai = f(t1);
figure; plot(t, cai);
%return;
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
K_init = load('data/initial-guess-k-hp2x7-cai.dat');
%previous_experiment = '1000';
[K_new, num_added] = load_initial_conditions_from_file_hp2x7_cai(K_init);
opt.K = K_new;
opt.initial_condition_index = length(K_init) + 1;
opt.x0 = zeros(9, 1);
%x0_new(1) = K(opt.initial_condition_index + 0);		% s1.
opt.x0(1, 1) = 1.0; % C
opt.x0 = load_initial_conditions_from_k_hp2x7_cai(opt.K, opt);
%opt.S = sparse(S);
opt.odeopt = odeopt;
opt.tspan = tspan;
opt.tint = tint;
opt.ATP = ATP;
opt.ATP_t1 = ATP_t1;
opt.ATP_t2 = ATP_t2;
opt.observable_index = 5; % Cytosolic concentration of calcium

%return;
global rt_fig fig_states_rt err_prev counter total_elapsed_time;
rt_fig = figure('Name', 'Real-time fitted calcium curve');
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
ops.LBounds = 0.0;
ops.UBounds = Inf;
ops.PopSize = pop_size;
opt.PopSize = ops.PopSize;
ops.StopFitness = 1e-16;
ops.DiagonalOnly = 1;
ops.MaxIter = 1000000;    % maximum number of iterations needed for estimation of parameters
ops.DispModulo = 1; % disp messages after every i-th iteration
ops.Resume = resume;
ops.SaveVariables = 'on';
ops.StopOnStagnation = 'off';
ops.TolUpX = '1e9*max(insigma)';
ops.TolFun = '1e-20 % stop if fun-changes smaller TolFun';
ops.tolhistfun = '-1e-20';
ops.SaveFilename = sprintf('variablescmaes-hp2x7-ATP-%d-popsize-%d-cai.mat', ATP, pop_size);
ops.LogFilenamePrefix = sprintf('outcmaes-hp2x7-ATP-%d-popsize-%d-cai  % files for output data', ATP, pop_size); 

w = cai;

sigma = 0.5;

[X, F, COUNT, STOP, OUT, BESTEVER]= cmaes('fhngen_hp2x7_fit_cai', y0, sigma, ops, opt, w);

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