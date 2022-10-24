%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function hp2x7_cai_prediction()
    close all; clc; clear; clearvars;
    format long g;
    ATP = 10e-3;%[0.01e-3 0.1e-3 1e-3 3e-3 10e-3];
	ATP_t1 = 2; ATP_t2 = 4.6; t_max = 10;% unit in seconds - Original
    %ATP_t1 = 30; ATP_t2 = 330; t_max = 400;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 32; t_max = 35;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 2.5; t_max = 8.5;% unit in seconds
    %ATP_t1 = 2*60; ATP_t2 = 8 * 60; t_max = 9 * 60;% unit in seconds
    %ATP_t1 = 1*60; ATP_t2 = 11 * 60; t_max = 12 * 60;% unit in seconds
    
    global nT_P2X T;
    T = 0; % Period in seconds
    nT_P2X = 0;

    es_hp2x7_cai = load('variablescmaes-hp2x7-ATP-5.000000e-03-popsize-72-total-current.mat');
	K_hp2x7_cai = es_hp2x7_cai.out.solutions.bestever.x;
    %K_hp2x7_cai = [1.523384e+05 4.679995e+03 4.149279e+02 6.874058e-01 1.942445e+02 5.770408e+01 7.289234e-01 3.116321e+03 1.598147e+00 2.494154e+07 ];
	k_hp2x7_cai = K_hp2x7_cai(1:6);
	opt_hp2x7_cai.x0 = zeros(4, 1);
    opt_hp2x7_cai.x0(1, 1) = 1; % [C]_0
    opt_hp2x7_cai.x0(4, 1) = 45 * 1e-9; % Initial condtion for intracellular/cytosolic Ca+2 concentration 70 nM
	opt_hp2x7_cai.initial_condition_index = length(k_hp2x7_cai) + 1;
	opt_hp2x7_cai.x0 = load_initial_conditions_from_k_hp2x7_total_current(K_hp2x7_cai, opt_hp2x7_cai);
	%opt_hp2x7_cai.x0
	%K_hp2x7_cai
	opt_hp2x7_cai.K = k_hp2x7_cai;
	opt_hp2x7_cai.ATP = ATP;
	opt_hp2x7_cai.ATP_t1 = ATP_t1;
	opt_hp2x7_cai.ATP_t2 = ATP_t2;
	opt_hp2x7_cai.tspan = [0.0 t_max];
    opt_hp2x7_cai.observable_index = 4;

    odeopt = odeset;%('RelTol', 1e-12, 'AbsTol', 1e-12);
	%odeopt = odeset;%('Jacobian', @compute_jacobian_matrix_hp2x7_total_current);
    sol_hp2x7_cai = ode15s(@(t, x)(reaction_network_hp2x7_cai_prediction(t, x, opt_hp2x7_cai)), opt_hp2x7_cai.tspan, opt_hp2x7_cai.x0, odeopt);

	tint = linspace(0, t_max, 1000);  % all time points.
	y_hp2x7_cai = deval(sol_hp2x7_cai, tint);

	fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2X7 receptor in human microglia');
	fig1.WindowStyle = 'docked';
	cai = y_hp2x7_cai(opt_hp2x7_cai.observable_index, :) * 1e9;
    i_tot = y_hp2x7_cai(3, :) * (2.5 / 2.416) * 145;
	hold on;
	plot(tint, cai*1e0, '-ko', 'MarkerIndices', 1:50:length(cai));
    %plot(tint, i_tot, '--b', 'MarkerIndices', 1:40:length(i_tot));
    %plot(tint, y_hp2x7_cai(1, :), '-ko', 'MarkerIndices', 1:40:length(O));
	plot([ATP_t1/1 ATP_t2/1], [44 44], '-k' , 'LineWidth', 2);
	hold off;
	text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 45, 'ATP');
    ATP = ATP * 1e3;
    legend(sprintf('\\Delta[Ca^2+_i] at ATP=%gmM', ATP));
    
    xlabel('Time (s)');
    %ylabel('\Delta[Ca^2+_i] (\muM)');
    ylabel('\Delta[Ca^2+_i] (nM)');
    %ylabel('\Delta[Ca^2+_i] (M)');
    %ylabel('\Delta[Q] (M)');
    xlim([0 t_max]);
    ylim([43 Inf]);%
end
%--------------------%