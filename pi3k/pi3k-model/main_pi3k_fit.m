maxNumCompThreads(1);
%parpool(2);
%delete(gcp('nocreate'));
%----- Prepare the Environment -------%
close all; clc; clear; clearvars;
format long g;
TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
c = Constants; % Gets constants for the model
resume = c.resume; % yes or no
pop_size = c.PopSize;
baseline = 100; %45
opt_p2x_cai.baseline = baseline;
opt.baseline = baseline;
%------Read data from file -------%
d = load('data/pAkt-ADP-50uM.dat');

offset = 20 / 60;
ADP = 50; % unit in micromollars
ADP_t1 = 0 + offset; ADP_t2 = 30 + offset; t_max = 2000/60;%t_max = 30 * 60 + offset;% unit in seconds
t =  d(:, 1) / 60;
pAkt = 0.2 * normalize(d(:, 2), 'range'); % pA 700
w = pAkt;
%--------------- Simulates the P2Y model ---------------%
es_hp2y12_cai = load(sprintf('variablescmaes-hp2y12-ADP-50-popsize-%d.mat', c.PopSize), '-mat');
K_hp2y12_cai = es_hp2y12_cai.out.solutions.bestever.x;
K_hp2y12_cai = read_last_line(c);
k = K_hp2y12_cai;
%k(46:50)
%return
ff = find(K_hp2y12_cai < 0);
if length(ff) > 1
    disp('A negative value found in K. So we set it to positive');
    %K_hp2y12_cai = abs(K_hp2y12_cai);
end

K_init = load('data/initial_guess_k_hp2y12.dat');
opt_p2y_cai.nn = length(K_init);

k_hp2y12_cai = K_hp2y12_cai(1:opt_p2y_cai.nn);
opt_p2y_cai.x0 = zeros(8, 1);
opt_p2y_cai.initial_condition_index = length(k_hp2y12_cai) + 1;
opt_p2y_cai.x0 = load_initial_conditions_from_k_hp2y12(K_hp2y12_cai, opt_p2y_cai);
%opt_p2y_cai.x0(14) = 0.045;
opt_p2y_cai.K = k_hp2y12_cai;
opt_p2y_cai.ADP = ADP;
opt_p2y_cai.ADP_t1 = ADP_t1 * 60; % Convert mintues to seconds
opt_p2y_cai.ADP_t2 = ADP_t2 * 60; % Convert mintues to seconds
opt_p2y_cai.tspan = [0.0 t_max * 60]; % Convert mintues to seconds
opt_p2y_cai.observable_index = 8;

odeopt = odeset;%('InitialStep', 1e-2, 'MaxStep', 1e-2);%('RelTol', 1e-9, 'AbsTol', 1e-9);
sol_p2y_cai = ode15s(@(t, x)(reaction_network_hp2y12(t, x, opt_p2y_cai)), opt_p2y_cai.tspan, opt_p2y_cai.x0, odeopt);
%--------------- simulate the P2X model ---------------%
% Sensitivity analysis (SA)
opt_p2x_cai.SA = zeros(2 + 2 * 8, 1);
opt_p2x_cai.SAncx = 1;
opt_p2x_cai.SApmca = 2;
opt_p2x_cai.SAk1f_p2x7 = 3;
opt_p2x_cai.SAk2f_p2x7 = 4;
opt_p2x_cai.SAk3f_p2x7 = 5;
opt_p2x_cai.SAk1b_p2x7 = 6;
opt_p2x_cai.SAk2b_p2x7 = 7;
opt_p2x_cai.SAk3b_p2x7 = 8;
opt_p2x_cai.SAk4b_p2x7 = 9;
opt_p2x_cai.SAg_p2x7 = 10;
opt_p2x_cai.SAk1f_p2x4 = 11;
opt_p2x_cai.SAk2f_p2x4 = 12;
opt_p2x_cai.SAk3f_p2x4 = 13;
opt_p2x_cai.SAk1b_p2x4 = 14;
opt_p2x_cai.SAk2b_p2x4 = 15;
opt_p2x_cai.SAk3b_p2x4 = 16;
opt_p2x_cai.SAk4b_p2x4 = 17;
opt_p2x_cai.SAg_p2x4 = 18;

[opt_p2x_cai.g_CA_Leak, opt_p2x_cai.g_NA_Leak] = compute_leak_conductances(opt_p2x_cai);

es_hp2x7_cai = load('variablescmaes-hp2x7-ATP-5.000000e-03-popsize-18-total-current.mat');
K_hp2x7_cai = es_hp2x7_cai.out.solutions.bestever.x;
es_hp2x4_cai = load('variablescmaes-hp2x4-ATP-1.000000e-04-popsize-18-total-current.mat');
K_hp2x4_cai = es_hp2x4_cai.out.solutions.bestever.x;
%es_rncx = load('variablescmaes-rncx.mat');
%K_rncx = es_rncx.out.solutions.bestever.x;
k_hp2x7_cai = K_hp2x7_cai;%(1:10);
k_hp2x4_cai = K_hp2x4_cai;%(1:7);
opt_p2x_cai.x0 = zeros(11, 1);
opt_p2x_cai.x0(1, 1) = 1; % [C]_0
opt_p2x_cai.x0(5, 1) = 1; % [C]_0
opt_p2x_cai.x0(9, 1) = 8e-3; % Initial condtion for intracellular/cytosolic Na+ concentration 8mM
opt_p2x_cai.x0(10, 1) = baseline * 1e-9; % Initial condtion for intracellular/cytosolic Ca+2 concentration 45nM
opt_p2x_cai.x0(11, 1) = -60e-3; % Initial condtion for membrane voltage
opt_p2x_cai.K_p2x7 = k_hp2x7_cai;
opt_p2x_cai.K_p2x4 = k_hp2x4_cai;
opt_p2x_cai.ATP = ADP * 1e-6; % uM to M
opt_p2x_cai.ATP_t1 = ADP_t1 * 60; % Convert mintues to seconds
opt_p2x_cai.ATP_t2 = ADP_t2 * 60; % Convert mintues to seconds
opt_p2x_cai.tspan = [0.0 t_max * 60]; % Convert mintues to seconds
opt_p2x_cai.observable_index = 10;
opt_p2x_cai.test = false;
opt_p2x_cai.CAiB = opt_p2x_cai.x0(10, 1);


odeopt = odeset;%('InitialStep', 1e-3, 'MaxStep', 1e-3);%('RelTol', 1e-9, 'AbsTol', 1e-9);
sol_p2x_cai = ode15s(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt_p2x_cai)), opt_p2x_cai.tspan, opt_p2x_cai.x0, odeopt);
%return;

%%{
%t_d = 2.7703 * 60; % minutes to seconds
%t_d = 10 * 60;
tint = linspace(opt_p2x_cai.tspan(1), opt_p2x_cai.tspan(2), 1000);  % all time points.
y_p2x_cai = deval(sol_p2x_cai, tint);
cai = y_p2x_cai(opt_p2x_cai.observable_index, :);
opt_p2x_cai.cai = cai;
opt_p2x_cai.tint_cai = tint;

%{
tint = linspace(opt_p2y_cai.tspan(1), opt_p2y_cai.tspan(2), 10000);  % all time points.
y_p2y_cai = deval(sol_p2y_cai, tint);
fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2Y12 receptor in human microglia');
fig1.WindowStyle = 'docked';
cai = y_p2y_cai(opt_p2y_cai.observable_index, :);
hold on;
plot(tint/60, cai, '-ko', 'MarkerIndices', 1:steps:length(cai));
%plot(tint/60, ATP * 0.005 , '-.r*', 'MarkerIndices', 1:steps:length(ATP));
%plot(tint, cai2, '-.r*', 'MarkerIndices', 1:steps:length(cai2));
%plot(tint, cai, '-ko', 'MarkerIndices', 1:100:length(cai));
%plot(tint, i_tot/1, '--b', 'MarkerIndices', 1:100:length(i_tot));
%plot(tint, cumtrapz(tint, i_tot) , '--b', 'MarkerIndices', 1:100:length(i_tot));
%plot(tint, y_hp2y12_cai(1, :), '-ko', 'MarkerIndices', 1:40:length(O));
plot([ADP_t1/1 ADP_t2/1], [0.2 0.2], '-k' , 'LineWidth', 2);
text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.22, 'ADP');
%text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.22, 'ATP');
hold off;
legend('[Ca^2+_i]', 'Agonist', 'FontWeight', 'bold');
xlabel('Time (m)');
ylabel('[Ca^2+_i] (\muM)');
%}
% -------Fitting -------------%
%t = t + 20;
%t
%[t pAkt]
%t = [0 0.5 1 3 5 10 30]' * 60;% + offset * ones(7, 1);
%d = [0 873 10812 8637 10917 7269 0]';
%d = [0 0 10812 8637 10917 7269 0]';
%t(1) = 0;
%shape = normalize(d, 'range') * 0.1;
%t

% Smoothing the data sets
%{
mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
f = fit(t, pAkt, mode);
t1 = linspace(0.0, t(length(t)), 1000);
%t1 = [linspace(0.0, 6, 1000) linspace(6.1, t(length(t)), 10)];
t = t1;
w = f(t1);
%int = integrate(f, t1, 0);
%figure;
%plot(t, w);
%xlabel('Time (s)');
%%ylabel('pAkt (\muM)');
%%xlim([0 t_max]);
%return
%}


mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
global f1 w_experimental first;
f1 = fit(t, w, mode);
if c.should_use_matlab_solver == 1
    tw = linspace(0.0, t(length(t)), 1000);
else
    tw = 0:0.5:t(length(t));
end
differentiate(f1, 17.55)
differentiate(f1, 40)
%return
w_experimental = f1(tw);
first = true;
%----- Constants -------%
%nt = length(t);     % number of time points to simulate
%tol = 1e-20;
odeopt = odeset;%('RelTol', tol, 'AbsTol', tol);
%----- Variables -------%
tspan = [0.0 t(length(t))];
tint = tw;  % all time points
%S = readmatrix('STM-P2X.xlsx'); % Loads the stoichiometry matrix. small model
%[NROWS, NCOLS] = size(S);
%---- Preparing for cureve fitting ------%
K_init = load('data/initial-guess-k-pi3k.dat');
%previous_experiment = '1000';
[K_new, num_added] = load_initial_conditions_from_file_pi3k(K_init);
opt.K = K_new;
opt.initial_condition_index = length(K_init) + 1;
opt.x0 = zeros(15, 1);
%opt.x0(1, 1) = 1;
opt.x0 = load_initial_conditions_from_k_pi3k(opt.K, opt);
%opt.S = sparse(S);
opt.odeopt = odeopt;
opt.tspan = tspan;
opt.tint = tint;
opt.ADP = ADP;
opt.ADP_t1 = ADP_t1;
opt.ADP_t2 = ADP_t2;
opt.observable_index = 10; % Cytosolic concentration of pAKT
opt.sol_p2y_cai = sol_p2y_cai;
%opt.mode = c.pi3k_mode;
opt.CAiB_p2y = opt_p2y_cai.x0(opt_p2y_cai.observable_index); % Steady-state baseline of CAi: units in uM
opt.CAiB_p2x = opt_p2x_cai.x0(10, 1); % Steady-state baseline of CAi: units in M
opt.cai_p2x = opt_p2x_cai.cai;
opt.tint_p2x = opt_p2x_cai.tint_cai;


%return;
global rt_fig fig_states_rt err_prev counter total_elapsed_time;
rt_fig = figure('Name', 'Real-time fitted current curve');
fig_states_rt = figure('Name', 'PI3K states before fitting');
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
y0 = [2.224115e+03 7.834268e-01 2.101810e+01 2.045164e+01 6.053654e+01 4.000991e+00 4.889444e+00 3.134684e+01 4.936186e-02 6.852647e-02 3.158896e+00 2.584806e-01 0.5 2.071029e-02 1.429151e+01 1.395122e+01 7.487444e+00 1.367922e+00 1.085065e+02 8.557009e-01 8.830612e+00 1.050761e+02 1.491477e+02 6.663650e+01 3.209582e+00 2.810433e+01 1.527106e+01 1.025415e+03 1.031701e+02 2.592506e+02 6.940814e+02 1.583591e+02 9.976268e+02];
%y0 = opt.K;          % Initial guess for k
%ops.LBounds = LB;      % Sets of parameters to be estimated
%opt.UBounds = UB; 
%for i=1:1:length(y0)
%    if i > length(opt.K) - 2
%        UB(i) = 1.0;
%    else
%        UB(i) = Inf;
%    end
%end
ops.LBounds = zeros(length(y0), 1);
%ops.LBounds(17) = 1;
%ops.LBounds(18) = 1;
ops.LBounds(25) = 2;
ops.LBounds(26) = 1;
ops.LBounds(27) = 10;
ops.UBounds = Inf * ones(length(y0), 1);
%ops.UBounds(17) = 8;
%ops.LBounds = 0.0 * ones(length(y0), 1);
%ops.UBounds = 1000 * ones(length(y0), 1);
ops.PopSize = pop_size;
opt.PopSize = ops.PopSize;
ops.StopFitness = 1e-16;
ops.DiagonalOnly = 1;
ops.MaxIter = 1000000;    % maximum number of iterations needed for estimation of parameters
ops.DispModulo = 1; % disp messages after every i-th iteration
ops.Resume = resume;
ops.SaveVariables = 'on';
ops.StopOnStagnation = 'off';
%ops.TolUpX = '1e9*max(insigma)';
%ops.TolFun = '1e-20 % stop if fun-changes smaller TolFun';
%ops.tolhistfun = '1e-20';
ops.SaveFilename = sprintf('variablescmaes-pi3k-ADP-%d-popsize-%d.mat', ADP, pop_size);
ops.LogFilenamePrefix = sprintf('outcmaes-pi3k-ADP-%d-popsize-%d-  % files for output data', ADP, pop_size); 

sigma = c.sigma * ones(length(y0), 1);
[X, F, COUNT, STOP, OUT, BESTEVER]= cmaes('fhngen_pi3k_fit', y0, sigma, ops, opt);%, w);

printk(X);
%------Functions -------%
function [g_CA_Leak, g_NA_Leak] = compute_leak_conductances(opt_p2x_cai)
    CAx = 2e-3;   % M
    %CAi = 45e-9;  % M
	CAi = opt_p2x_cai.baseline * 1e-9;  % M
    NAi = 8e-3;   %M
    NAx = 130e-3; % M
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    Vm = -60e-3;
    Z_Ca = 2;
    Z_Na = 1;
    Ena =  R * T * Z_Na^-1 * F^-1 * log(NAx/NAi);
    Eca =  R * T * Z_Ca^-1 * F^-1 * log(CAx/CAi);
    
    I_NA_NCX = NA_NCX_Current(NAi, CAi, Vm, opt_p2x_cai);
    I_CA_NCX = -2 * 3^-1 * I_NA_NCX;

    g_NA_Leak = (-I_NA_NCX) / (Vm - Ena);
    g_CA_Leak = (-I_CA_NCX - PMCA_Current(CAi, opt_p2x_cai)) /  (Vm - Eca);
end
%--------------------%
function [K] = read_last_line(c)
    fid = fopen(sprintf('outcmaes-hp2y12-ADP-50_popsize-%d-xmean.dat', c.PopSize),'r');     %# Open the file as a binary
    lastLine = '';                   %# Initialize to empty
    offset = 1;                      %# Offset from the end of file
    fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
    newChar = fread(fid,1,'*char');  %# Read one character
    while (~strcmp(newChar,char(10))) || (offset == 1)
      lastLine = [newChar lastLine];   %# Add the character to a string
      offset = offset+1;
      fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
      newChar = fread(fid,1,'*char');  %# Read one character
    end
    fclose(fid);  %# Close the file
    array = split(lastLine, " ");
    K = zeros(length(array) - 6, 1);
    for i=6:1:(length(array)-1)
        K(i - 5) = str2double(array(i));
    end
end
%--------------------%
function printk(k)
    m = length(k);
    s = "";
    for i=1:1:m
        s = sprintf("%s %g", s, k(i));
    end
    disp(s);
end
%--------------------%