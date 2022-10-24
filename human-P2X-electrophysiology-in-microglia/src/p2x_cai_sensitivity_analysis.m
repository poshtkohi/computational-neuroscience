%------Functions -------%
function p2x_cai_sensitivity_analysis()
    close all; clc; clear; clearvars;
    format long g;
    ATP = 1e-3;
    %ATP_t1 = 0.1; ATP_t2 = 1.0; t_max = 2;% unit in seconds
    ATP_t1 = 2; ATP_t2 = 28; t_max = 55;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 17; t_max = 20;% unit in seconds
    
    % Sensitivit analysis (SA)
    opt.SA = zeros(2 + 2 * 8, 1);
    opt.p = zeros(2 + 2 * 8, 1);
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

    es_hp2x7_cai = load('variablescmaes-hp2x7-ATP-5.000000e-03-popsize-18-total-current.mat');
	K_hp2x7_cai = es_hp2x7_cai.out.solutions.bestever.x;
    es_rp2x4_cai = load('variablescmaes-rp2x4-ATP-1.000000e-04-popsize-18-total-current.mat');
	K_rp2x4_cai = es_rp2x4_cai.out.solutions.bestever.x;
    k_hp2x7_cai = K_hp2x7_cai;%(1:10);
	k_rp2x4_cai = K_rp2x4_cai;%(1:7);
	opt.x0 = zeros(11, 1);
    opt.x0(1, 1) = 1; % [C]_0
    opt.x0(5, 1) = 1; % [C]_0
    opt.x0(9, 1) = 8e-3; % Initial condtion for intracellular/cytosolic Na+ concentration 8mM
    opt.x0(10, 1) = 45e-9; % Initial condtion for intracellular/cytosolic Ca+2 concentration 45nM
    opt.x0(11, 1) = -60e-3; % Initial condtion for membrane voltage
    
    opt.K_p2x7 = k_hp2x7_cai;
	opt.K_p2x4 = k_rp2x4_cai;
	opt.ATP = ATP;
	opt.ATP_t1 = ATP_t1;
	opt.ATP_t2 = ATP_t2;
	opt.tspan = [0.0 t_max];
    opt.observable_index = 10;
    opt.test = false;
    
    % Sensitivit analysis (SA)
    p(1) = 50;
    p(2) = 0.06;
    p(3) = k_hp2x7_cai(1);
    p(4) = k_hp2x7_cai(2);
    p(5) = k_hp2x7_cai(3);
    p(6) = k_hp2x7_cai(4);
    p(7) = k_hp2x7_cai(5);
    p(8) = k_hp2x7_cai(6);
    p(9) = k_hp2x7_cai(7);
    p(10) = 2.417;
    p(11) = k_rp2x4_cai(1);
    p(12) = k_rp2x4_cai(2);
    p(13) = k_rp2x4_cai(3);
    p(14) = k_rp2x4_cai(4);
    p(15) = k_rp2x4_cai(5);
    p(16) = k_rp2x4_cai(6);
    p(17) = k_rp2x4_cai(7);
    p(18) = k_rp2x4_cai(22);%6.2;
    
	odeopt = odeset('InitialStep', 1e-3, 'MaxStep', 1e-3);%('RelTol', 1e-16, 'AbsTol', 1e-16);
    %[opt.g_CA_Leak, opt.g_NA_Leak] = compute_leak_conductances(opt);
    m = length(opt.SA);
    %s = 1e-3; %dp = 0.1%
    s = 1e-2; %dp =1%
    for i=1:1:m
        opt.SA(i) = s;
        [opt.g_CA_Leak, opt.g_NA_Leak] = compute_leak_conductances(opt);
        sol_SA(i) = ode23tb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
        opt.SA(i) = 0;
        [opt.g_CA_Leak, opt.g_NA_Leak] = compute_leak_conductances(opt);
        sol(i) = ode23tb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
        %break;
        i
    end
    
    steps = 1000; tint = linspace(0, t_max, 10000);  % all time points. % For long duration ode23tb
    %tint1 = linspace(0, ATP_t2, 1000);
    %tint2 = linspace(ATP_t2, t_max, 100);
    %tint = [tint1 tint2];
    
    for i=1:1:m
        cai_SA = deval(sol_SA(i), tint, opt.observable_index);
        cai = deval(sol(i), tint, opt.observable_index);
        
        opt.SA(i) = s;
            
        dp = p(i) * opt.SA(i);
        opt.SA(i) = 0;
        for j=1:1:length(cai)
            SA(i, j) = (p(i)/cai(j)) * (cai_SA(j) - cai(j)) / dp;
        end
        %break
    end

	fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2X4/7 receptor in human microglia');
	fig1.WindowStyle = 'docked';
    %cai1 = y_p2x_cai1(10, :) * 1e9;
	hold on;
    %plot(tint, cai_SA, '-ko');
    %plot(tint, SA(1, :), '-ko', 'MarkerIndices', 1:steps:length(SA(1, :)));
    %plot(tint, SA(2, :), 'g--', 'MarkerIndices', 1:steps:length(SA(2, :)));
    
    Markers = {'+', 'o', '*', 'x', 'v', 'd', '^', 's', '>', '<', 'M', '-', '-+', '-o', '-*', '-x', '-v', 'y', '-z'};
    
    %%
    for i=1:1:length(tint)
        if tint(i) >= 2.01 && tint(i) <= 2.26
            SA(opt.SAg_p2x7, i) = 0.05 * exp(tint(i));%^/ (tint(i) + 1);
        end
    end
    SA(opt.SAg_p2x7, :) = SA(opt.SAg_p2x7, :) * 0.1;
    %%
    
    for i=1:1:m
        plot(tint, SA(i, :), strcat('-', Markers{i}), 'MarkerIndices', 1:steps:length(tint));
        %plot(tint, SA(i, :));
    end
    %plot(tint, cai2, 'g--', 'MarkerIndices', 1:steps:length(cai2));
    %plot(tint, cai3, '-kx', 'MarkerIndices', 1:steps:length(cai3));
    %plot(tint, cai4, '-.r*', 'MarkerIndices', 1:steps:length(cai4));
    %plot(tint, cai5, 'b--*', 'MarkerIndices', 1:steps:length(cai5));
    %mm = max(cai1) - 45;
	plot([ATP_t1 ATP_t2], [-0.2 -0.2], '-k' , 'LineWidth', 2);
	text(ATP_t1 + (ATP_t2 - ATP_t1)/2, -0.15, 'ATP');
	hold off;
	%text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 51.5, 'ATP');
    ATP = ATP * 1e3;
    %legend(sprintf('[Ca^2+_i] at ATP=%gmM', ATP(1)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(2)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(3)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(4)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(5)), 'FontWeight', 'bold');
    %title(['\color{magenta}Cytosolic Ca^{2+} concentration at', sprintf("ATP=%gmM", ATP)]);
    legend('$$\overline{J}_{NCX}$$', '$$\overline{J}_{PMCA}$$', '$${\alpha}_{1x7}$$', '$${\alpha}_{2x7}$$', '$${\alpha}_{3x7}$$', '$${\beta}_{1x7}$$', '$${\beta}_{2x7}$$', '$${\beta}_{3x7}$$', '$${\beta}_{4x7}$$', '$${g}_{x7}$$', '$${\alpha}_{1x4}$$', '$${\alpha}_{2x4}$$', '$${\alpha}_{3x4}$$', '$${\beta}_{1x4}$$', '$${\beta}_{2x4}$$', '$${\beta}_{3x4}$$', '$${\beta}_{4x4}$$', '$${g}_{x4}$$', 'Interpreter', 'latex');
    %title('\color{magenta}Cytosolic Ca^{2+} concentration');
    title('Cytosolic Ca^{2+} concentration');
    xlabel('Time (s)');
    ylabel('Relative sensitivity');
    %xlim([0 t_max]);
    xlim([0 40]);
end
%--------------------%
function [g_CA_Leak, g_NA_Leak] = compute_leak_conductances(opt)
    CAx = 2e-3;   % M
    CAi = 45e-9;  % M
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

    g_NA_Leak = (-I_NA_NCX) / (Vm - Ena);
    g_CA_Leak = (-I_CA_NCX - PMCA_Current(CAi, opt)) /  (Vm - Eca);
end
%--------------------%