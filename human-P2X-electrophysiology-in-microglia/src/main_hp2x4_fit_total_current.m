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
resume = c.resume; % yes or no
pop_size = c.PopSize;
%------Read data from file -------%
d = load('data/hP2X4-ATP-0.1mM-total-current.dat');

ATP = 0.1e-3; % unit in mollars
%ATP_t1 = 19.56; ATP_t2 = 19.56 + 30.43; t_max = 70.0;% unit in seconds
ATP_t1 = 15.68; ATP_t2 = 25.255; t_max = 53;% unit in seconds
t = d(:, 1);
i_tot = d(:, 2);% * 1e-3; % pA
i_shape = normalize(-i_tot, 'range') * 372;% * 6.2 / 6.8;
%cai_shape = cumtrapz(t, i_shape);
%plot(t, i_shape, t, cai_shape);
%figure;
%plot(t, i_shape);
%return

global nT_P2X T;
T = 0; % Period in seconds
nT_P2X = 0;

% Smoothing the data sets
mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
f = fit(t, i_shape, mode);
t1 = linspace(0.0, t(length(t)), 200);
t = t1;
i_shape = f(t1);
%int = integrate(f, t1, 0);
%figure;
%plot(t, i_shape);
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
K_init = load('data/initial-guess-k-hp2x4-total-current.dat');
%previous_experiment = '1000';
[K_new, num_added] = load_initial_conditions_from_file_hp2x4_total_current(K_init);
opt.K = K_new;
opt.initial_condition_index = length(K_init) + 1;
opt.x0 = zeros(4, 1);
opt.x0(1, 1) = 1;
opt.x0 = load_initial_conditions_from_k_hp2x4_total_current(opt.K, opt);
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
fig_states_rt = figure('Name', 'rP2X4 states before fitting');
fig_states_rt.WindowStyle = 'docked';
err_prev = Inf;
counter = 1;
total_elapsed_time = 0.0;
rt_fig.WindowState = 'maximized';
rt_fig.Units = 'normalized';
rt_fig.OuterPosition = [0 0 1 1];
rt_fig.WindowStyle = 'docked';
%----- Curve Fitting -------%
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

%ops.LBounds(25) = 1;%
%ops.UBounds(25) = 40;%

%{
ops.LBounds(25) = 49;%
ops.LBounds(26) = 49;%
ops.LBounds(27) = 49;%

ops.UBounds(25) = 51;%
ops.UBounds(26) = 51;%
ops.UBounds(27) = 51;%
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

%ops.LBounds(3:6) = 1;%open rp2
%ops.UBounds(3:6) = 20;

y0 = [1.870575e+05 1.111369e+03 6.223432e+00 7.297772e+00 1.936948e+01 1.340619e-02 4.804276e+04 1.783911e+04 2.769054e+03 2.884102e+02 5.062157e+02 2.453336e+04 3.434921e+03 6.176162e+03 3.122662e+04 2.851495e+05 8.622153e+03 8.525302e+02 5.387433e+02 1.751705e+04 1.803310e+02 6.515703e+03 6.890038e+03 2.073169e+02 4.714344e+04 5.768429e+03 1.540326e+03 5.280772e+02];

ops.PopSize = pop_size;
opt.PopSize = ops.PopSize;
ops.StopFitness = 1e-4;
ops.DiagonalOnly = 1;
ops.MaxIter = 1000000;    % maximum number of iterations needed for estimation of parameters
ops.DispModulo = 1; % disp messages after every i-th iteration
ops.Resume = resume;
ops.SaveVariables = 'on';
ops.StopOnStagnation = 'off';
ops.TolUpX = '1e9*max(insigma)';
%ops.TolFun = '1e-20 % stop if fun-changes smaller TolFun';
%ops.tolhistfun = '1e-20';
ops.SaveFilename = sprintf('variablescmaes-hp2x4-ATP-%d-popsize-%d-total-current.mat', ATP, pop_size);
ops.LogFilenamePrefix = sprintf('outcmaes-hp2x4-ATP-%d-popsize-%d-total-current  % files for output data', ATP, pop_size); 

w = i_shape;

sigma = 1.0 * ones(length(y0), 1);
%sigma(1:2) = 1;
%sigma = 0.5;
%sigma(25) = 0.01;

[X, F, COUNT, STOP, OUT, BESTEVER]= cmaes('fhngen_hp2x4_fit_total_current', y0, sigma, ops, opt, w);

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