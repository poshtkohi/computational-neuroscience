%----- Prepare the Environment -------%
close all; clc; clear; clearvars;
format long;
fclose('all');
c = Constants; % Gets constants for the model
%------Read data from file -------%
offset = 0.0;
ADP_t1 = offset + 13.88; ADP_t2 = offset + 13.88 + 29.16; max_t_p2y12 = offset + 66.67;
ADP = 50; %uM
resume = c.resume; % yes or no
pop_size = c.PopSize;
%global nT_P2Y12 T;
%T = 0; % Period in minutes
%nT_P2Y12 = 0;

% hP2Y12 data
d = load('data/hP2Y12-ADP-50uM-CAi.dat');
t = d(:, 1) + offset;

w = normalize(d(:, 2), 'range') * 0.5 + 0.1; % d in uM;
% Smoothing the data sets
mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
global f1 w_experimental first;
f1 = fit(t, w, mode);
if c.should_use_matlab_solver == 1
    tw = linspace(0.0, t(length(t)), 1000);
else
    tw = 0:0.5:t(length(t));
end
w_experimental = f1(tw);
first = true;

%{
figure('Name', 'Experimental CAi for the human P2Y12 receptor');
[d1, d2] = differentiate(f1, tw);
hold on
plot(tw, -d1); xlabel('Time (s)'); ylabel('CAi (uM)');
plot(tw, w); xlabel('Time (s)'); ylabel('CAi (uM)');
hold off
%xlim([0 3]);
legend('hP2Y12');
%}

%return


%figure('Name', 'Experimental CAi');
%plot(tw, w); xlabel('Time (s)'); ylabel('CAi (uM)');
%----- Constants -------%
c = Constants; % Gets constants for model
nt = length(t);     % number of time points to simulate
%nt
%tol = 1e-9;
odeopt = odeset;%('InitialStep', 1e-3, 'MaxStep', 1e-3);%('RelTol', tol, 'AbsTol', tol);
%----- Variables -------%
x0 = zeros(8, 1);
tspan = [tw(1) tw(length(tw))];
tint = tw;  % all time points
%tint = linspace(0, c.max_t_P2Y12, nt)';  % all time points
%S = readmatrix('STM-hP2Y12.xlsx'); % Loads the stoichiometry matrix
%[NROWS,NCOLS] = size(S);
%---- Boundary Conditions ------%
%x0 = c.X0_init;
global K_init;
K_init = load('data/initial_guess_k_hp2y12.dat');
%K_init = load('data/initial_guess_k_p2y12-manipulated1.dat');
%opt.S = sparse(S);
opt.K = K_init;%%
opt.odeopt = odeopt;
opt.x0 = x0;
opt.tspan = tspan;
opt.tint = tint;
opt.ADP = ADP;
opt.ADP_t1 = ADP_t1;
opt.ADP_t2 = ADP_t2;
opt.observable_index = 8; % Cytosolic concentration of calcium
opt.nn = length(K_init);
[K_new, num_added] = load_initial_conditions_from_file_hp2y12(opt.K);
opt.K = K_new;
%global is;
%is = true;
opt.initial_condition_index = length(K_init) + 1;
opt.x0 = load_initial_conditions_from_k_hp2y12(opt.K, opt);
is = false;
opt.mode = c.curefitting_mode;
global rt_fig fig_states_rt err_prev counter total_elapsed_time;
rt_fig = figure('Name', 'Real-time fitted calcium curve');
err_prev = Inf;
counter = 1;
total_elapsed_time = 0.0;
rt_fig.WindowState = 'maximized';
rt_fig.Units = 'normalized';
rt_fig.OuterPosition = [0 0 1 1];
rt_fig.WindowStyle = 'docked';

s = sprintf('hP2Y12 states for ADP %g uM', ADP);
fig_states_rt = figure('Name', s);
fig_states_rt.WindowState = 'maximized';
fig_states_rt.Units = 'normalized';
fig_states_rt.OuterPosition = [0 0 1 1];
fig_states_rt.WindowStyle = 'docked';
%----- Curve Fitting -------%
% Parameter estimation using GA (genetics algorithm) 
y0 = opt.K;          % Initial guess for k
ops.LBounds = zeros(length(y0), 1);

ops.UBounds = Inf * ones(length(y0), 1);

ops.PopSize = pop_size;
opt.PopSize = ops.PopSize;
ops.StopFitness = 1e-20;
ops.DiagonalOnly = 1;
ops.MaxIter = 1000000;    % maximum number of iterations needed for estimation of parameters
ops.DispModulo = 1; % disp messages after every i-th iteration
ops.Resume = resume;
ops.SaveVariables = 'on';
ops.StopOnStagnation = 'off';
ops.StopOnWarnings = 'off';
ops.TolUpX = '1e20*max(insigma)';
%ops.TolX = '1e-12*max(insigma) % stop if x-change smaller TolX';
ops.SaveFilename = sprintf('variablescmaes-hp2y12-ADP-%d-popsize-%d.mat', ADP, pop_size);
ops.LogFilenamePrefix = sprintf('outcmaes-hp2y12-ADP-%d_popsize-%d- ', ADP, pop_size);

sigma = 1.0 * ones(length(y0), 1);
sigma(13) = 0.01;
sigma(26) = 0.01;


[X, F, COUNT, STOP, OUT, BESTEVER]= cmaes('fhngen_hp2y12_fit', y0, sigma, ops, opt);%, w_experimental);
printk(X);
kest = X;
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