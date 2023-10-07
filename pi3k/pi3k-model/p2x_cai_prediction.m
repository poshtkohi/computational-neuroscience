%------Functions -------%
%function p2x_cai_prediction()
    close all; clc; clear; clearvars;
    format long g;
    c = Constants; % Gets constants for the model
    %%ATP = [0.05e-3 0.1e-3 0.5e-3 1.0e-3 3.0e-3];
    %ATP = [0 0.03e-3 0.04e-3 0.05e-3 0.10e-3];
    ATP = [0 0.05e-3 0.05e-3 0.05e-3 0.05e-3];
    %ATP = [0.05e-3];
    %ATP = [0.01e-3 0.1e-3 0.2e-3 0.3e-3 0.4e-3];
    %ATP = [1e-3 1e-3 1e-3 1e-3 1e-3];
    %ATP = [0.01e-3 3e-3 1e-3 3.0e-3 10e-3];
    %ATP = [0.01e-3  0.1e-3 0.2e-3 0.3e-3 1e-3];
    %ATP = [0.01e-3 0.1e-3 1e-3 3.0e-3 10e-3];
    %ATP_t1 = 1; ATP_t2 = 3.6; t_max = 6;% unit in seconds
    %ATP_t1 = 0.1; ATP_t2 = 0.6; t_max = 1;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 32; t_max = 40;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 122; t_max = 200;% unit in seconds
    %ATP_t1 = 10.0; ATP_t2 = (116.0 + 10.0); t_max = 120;% unit in seconds
    %ATP_t1 = 80.0; ATP_t2 = 480; t_max = 540;% unit in seconds
    %ATP_t1 = 20; ATP_t2 = 350; t_max = 360;% unit in seconds
    %ATP_t1 = 20; ATP_t2 = 140; t_max = 160;% unit in seconds
    %ATP_t1 = 20; ATP_t2 = 80.0; t_max = 100;% unit in seconds
    %ATP_t1 = 10; ATP_t2 = 40.0; t_max = 50;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 17; t_max = 300;% unit in seconds
    %ATP_t1 = 5; ATP_t2 = 20; t_max = 25;% unit in seconds
    %ATP_t1 = 5; ATP_t2 = 30; t_max = 35;% unit in seconds
    ATP_t1 = 2; ATP_t2 = 17; t_max = 200;% unit in seconds
    offset = 20;
    ATP_t1 = 0 + offset; ATP_t2 = 30 * 60 + offset; t_max = 2000 + 9000;%t_max = 30 * 60 + offset;% unit in seconds
    
    baseline = 100;
    opt.baseline = baseline;
    %----------------------
    global t1 t2 g enabled;
    g = 0.5;
    enabled = 1;
    t1 = ATP_t1;
    t2 = ATP_t2;
    t_max = t_max;
    opt.delay = c.delay * 60;
    opt.width = t2 - t1;
    opt.amplitude = ATP;
    %----------------------
    % Sensitivity analysis (SA)
    opt.SA = zeros(2 + 2 * 8, 1);
    opt.SAncx = 1;
    opt.SApmca = 2;
    opt.SAk1f_p2x7 = 3;
    opt.SAk2f_p2x7 = 4;
    opt.SAk3f_p2x7 = 5;
    opt.SAk1b_p2x7 = 6;
    opt.SAk2b_p2x7 = 7;
    opt.SAk3b_p2x7 = 8;
    opt.SAk4b_p2x7 = 9;
    opt.SAg_p2x7 = 10;
    opt.SAk1f_p2x4 = 11;
    opt.SAk2f_p2x4 = 12;
    opt.SAk3f_p2x4 = 13;
    opt.SAk1b_p2x4 = 14;
    opt.SAk2b_p2x4 = 15;
    opt.SAk3b_p2x4 = 16;
    opt.SAk4b_p2x4 = 17;
    opt.SAg_p2x4 = 18;
    
    [opt.g_CA_Leak, opt.g_NA_Leak] = compute_leak_conductances(opt);

    es_hp2x7_cai = load('variablescmaes-hp2x7-ATP-5.000000e-03-popsize-18-total-current.mat');
	K_hp2x7_cai = es_hp2x7_cai.out.solutions.bestever.x;
    es_hp2x4_cai = load('variablescmaes-hp2x4-ATP-1.000000e-04-popsize-18-total-current.mat');
	K_hp2x4_cai = es_hp2x4_cai.out.solutions.bestever.x;
    %es_rncx = load('variablescmaes-rncx.mat');
	%K_rncx = es_rncx.out.solutions.bestever.x;
    k_hp2x7_cai = K_hp2x7_cai;%(1:10);
	k_hp2x4_cai = K_hp2x4_cai;%(1:7);
	opt.x0 = zeros(11, 1);
    opt.x0(1, 1) = 1; % [C]_0
    opt.x0(5, 1) = 1; % [C]_0
    opt.x0(9, 1) = 8e-3; % Initial condtion for intracellular/cytosolic Na+ concentration 8mM
    opt.x0(10, 1) = baseline * 1e-9; % Initial condtion for intracellular/cytosolic Ca+2 concentration 45nM
    opt.x0(11, 1) = -60e-3; % Initial condtion for membrane voltage
    %opt.initial_condition_index = length(k_hp2x7_cai) + 1;
	%opt.x0 = load_initial_conditions_from_k_hp2x7_total_current(K_hp2x7_cai, opt);
	%opt.initial_condition_index = length(k_hp2x4_cai) + 1;
	%opt.x0 = load_initial_conditions_from_k_hp2x4_total_current(K_hp2x4_cai, opt);
	%opt.x0
	%K_hp2x4_cai
    opt.K_p2x7 = k_hp2x7_cai;
	opt.K_p2x4 = k_hp2x4_cai;
	opt.ATP = ATP;
	opt.ATP_t1 = ATP_t1;
	opt.ATP_t2 = ATP_t2;
	opt.tspan = [0.0 t_max];
    opt.observable_index = 10;
    opt.test = false;
    %opt.K_rncx = K_rncx;
    %opt.x0
    
    %odeopt = odeset;%('RelTol', 1e-12, 'AbsTol', 1e-12);
	%odeopt = odeset;%('Iacobian', @compute_Iacobian_matrix_hp2x4_total_current);
    %sol_p2x_cai = ode23s(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
    
    global my_k_p2x7 my_k_p2x4 I my_ATP;
	my_k_p2x7 = opt.K_p2x7;
	my_k_p2x4 = opt.K_p2x4;
    I = zeros(9, 9);
	%odeopt = odeset;%('RelTol', 1e-12, 'AbsTol', 1e-12);
    odeopt = odeset;%('InitialStep', 1e-1, 'MaxStep', 1e-1);%('RelTol', 1e-16, 'AbsTol', 1e-16);
	%odeopt = odeset('Iacobian', @compute_Iacobian_matrix_p2x);%, 'RelTol', 1e-12, 'AbsTol', 1e-12);
    for i=1:1:length(ATP)
        opt.ATP = ATP(i);
        opt.amplitude = ATP(i);
        g = 0.5;
        enabled = 1;
        t1 = ATP_t1;
        t2 = ATP_t2;
        my_ATP = ATP(i);
        sol_p2x_cai(i) = ode23tb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
        %sol_p2x_cai(i) = ode23ttb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
    end

	%tint = linspace(0, t_max, 1000);  % all time points.
	%y_p2x_cai = deval(sol_p2x_cai, tint);
    

    %steps = 100; tint = linspace(0, t_max, 1000);  % all time points. % For short duration ode23t
    steps = 10000; tint = linspace(0, t_max, 100000);  % all time points. % For long duration ode23tb
    %tint1 = linspace(0, ATP_t2, 2000);  % all time points.
    %tint2 = linspace(ATP_t2, t_max, 50);  % all time points.
    %tint = [tint1 tint2]';
    
	y_p2x_cai1 = deval(sol_p2x_cai(1), tint);
    %%y_p2x_cai2 = deval(sol_p2x_cai(2), tint);
    %%y_p2x_cai3 = deval(sol_p2x_cai(3), tint);
    %%y_p2x_cai4 = deval(sol_p2x_cai(4), tint);
    %%y_p2x_cai5 = deval(sol_p2x_cai(5), tint);

	fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2X4/7 receptor in human microglia');
	fig1.WindowStyle = 'docked';
    cai1 = y_p2x_cai1(10, :) * 1e9;
    %%cai2 = y_p2x_cai2(10, :) * 1e9;
    %%cai3 = y_p2x_cai3(10, :) * 1e9;
    %%cai4 = y_p2x_cai4(10, :) * 1e9;
    %%cai5 = y_p2x_cai5(10, :) * 1e9;
    %%caii = 0.7 * y_p2x_cai5(3, :) + 1 * y_p2x_cai5(6, :);
    %i_tot_hp2x7 = y_p2x_cai(3, :) * (2.5 / 2.416) * 145;
    %i_tot_hp2x4 = y_p2x_cai(6, :) * (7.0 / 6.2) * 372;
    %i_tot = i_tot_hp2x7 + i_tot_hp2x4;
	hold on;

    %tint = tint / 60.0;
    gg = 1;
    alpha = 0;%- gg * 45 + 45;
    cai1x = cai1 * gg + alpha;
    %%cai2x = cai2 * gg + alpha;
    %%cai3x = cai3 * gg + alpha;
    %%cai4x = cai4 * gg + alpha;
    %%cai5x = cai5 * gg + alpha;
    %plot(tint, cai2x, 'g--', 'MarkerIndices', 1:steps:length(cai2));
    %plot(tint, cai3x, '-kx', 'MarkerIndices', 1:steps:length(cai3));
    %plot(tint, cai4x, '-.r*', 'MarkerIndices', 1:steps:length(cai4));
    %plot(tint, cai5x, 'b--*', 'MarkerIndices', 1:steps:length(cai5));
    if c.periodic == 1
        g = 0.5;
        enabled = 1;
        t1 = ATP_t1;
        t2 = ATP_t2;
        A = zeros(length(tint), 1);
        G = zeros(length(tint), 1);
        for i=1:1:length(tint)
            %[A(i), G(i)] = periodic_function(tint(i), opt.delay, opt.width, 50);
            [A(i), G(i)] = periodic_function(tint(i), opt.ATP_t1, opt.ATP_t2, opt.delay, opt.amplitude);
        end
        %%plot(tint, normalize(cai1x, 'range'), '-ko', 'MarkerIndices', 1:steps:length(cai1));
        %%plot(tint, 0.5 * normalize(A, 'range'), 'b--*', 'MarkerIndices', 1:steps:length(cai1));
        plot(tint, G, 'g--', 'MarkerIndices', 1:steps:length(cai1));
        %plot(tint, normalize(G, 'range'), 'g--', 'MarkerIndices', 1:steps:length(cai1));
        legend('[Ca^2+_i]', 'ADP', 'g', 'FontWeight', 'bold');
    else
        mm = max(cai1x) - baseline;
        plot(tint, cai1x, '-ko', 'MarkerIndices', 1:steps:length(cai1));
        plot([ATP_t1 ATP_t2], [baseline + mm/2 baseline + mm/2], '-k' , 'LineWidth', 2);
        text(ATP_t1 + (ATP_t2 - ATP_t1)/2, baseline + mm/2 + 0.05 * mm/2, 'ADP');
    end
	hold off;
	%text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 51.5, 'ATP');
    ATP = ATP * 1e3;
    %legend(sprintf('[Ca^2+_i] at ADP=%gmM', ATP(1)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(2)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(3)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(4)), sprintf('[Ca^2+_i] at ADP=%gmM', ATP(5)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    %ylabel('\Delta[Ca^2+_i] (\muM)');
    ylabel('[Ca^2+_i] (nM)');
    %ylabel('[Ca^2+_i] (M)');
    %ylabel('[Q] (M)');
    xlim([0 t_max]);
%end
%--------------------%
function [g_CA_Leak, g_NA_Leak] = compute_leak_conductances(opt)
    CAx = 2e-3;   % M
    CAi = opt.baseline * 1e-9;  % M
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
    
    I_NA_NCX = NA_NCX_Current(NAi, CAi, Vm, opt);
    I_CA_NCX = -2 * 3^-1 * I_NA_NCX;

    g_NA_Leak = (-I_NA_NCX) / (Vm - Ena)
    g_CA_Leak = (-I_CA_NCX - PMCA_Current(CAi, opt)) /  (Vm - Eca)
end
%--------------------%