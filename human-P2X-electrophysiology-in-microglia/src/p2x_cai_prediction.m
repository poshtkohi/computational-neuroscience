%------Functions -------%
function p2x_cai_prediction()
    close all; clc; clear; clearvars;
    format long g;
    ATP = [0.05e-3 0.1e-3 0.5e-3 1.0e-3 3.0e-3];
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
    ATP_t1 = 2; ATP_t2 = 17; t_max = 30;% unit in seconds
    %ATP_t1 = 5; ATP_t2 = 20; t_max = 25;% unit in seconds
    %ATP_t1 = 5; ATP_t2 = 30; t_max = 35;% unit in seconds
    
    baseline = 100;
    opt.baseline = baseline;
    
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
    
    global nT_P2X T;
    T = 0; % Period in seconds
    nT_P2X = 0;

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
    odeopt = odeset;%('InitialStep', 1e-3, 'MaxStep', 1e-3);%('RelTol', 1e-16, 'AbsTol', 1e-16);
	%odeopt = odeset('Iacobian', @compute_Iacobian_matrix_p2x);%, 'RelTol', 1e-12, 'AbsTol', 1e-12);
    for i=1:1:length(ATP)
        opt.ATP = ATP(i);
        my_ATP = ATP(i);
        sol_p2x_cai(i) = ode23tb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
        %sol_p2x_cai(i) = ode23t(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
    end

	%tint = linspace(0, t_max, 1000);  % all time points.
	%y_p2x_cai = deval(sol_p2x_cai, tint);
    

    %steps = 100; tint = linspace(0, t_max, 1000);  % all time points. % For short duration ode23t
    steps = 1000; tint = linspace(0, t_max, 10000);  % all time points. % For long duration ode23tb
    %tint1 = linspace(0, ATP_t2, 2000);  % all time points.
    %tint2 = linspace(ATP_t2, t_max, 50);  % all time points.
    %tint = [tint1 tint2]';
    
	y_p2x_cai1 = deval(sol_p2x_cai(1), tint);
    y_p2x_cai2 = deval(sol_p2x_cai(2), tint);
    y_p2x_cai3 = deval(sol_p2x_cai(3), tint);
    y_p2x_cai4 = deval(sol_p2x_cai(4), tint);
    y_p2x_cai5 = deval(sol_p2x_cai(5), tint);

	fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2X4/7 receptor in human microglia');
	fig1.WindowStyle = 'docked';
    cai1 = y_p2x_cai1(10, :) * 1e9;
    cai2 = y_p2x_cai2(10, :) * 1e9;
    cai3 = y_p2x_cai3(10, :) * 1e9;
    cai4 = y_p2x_cai4(10, :) * 1e9;
    cai5 = y_p2x_cai5(10, :) * 1e9;
    %%caii = 0.7 * y_p2x_cai5(3, :) + 1 * y_p2x_cai5(6, :);
    %i_tot_hp2x7 = y_p2x_cai(3, :) * (2.5 / 2.416) * 145;
    %i_tot_hp2x4 = y_p2x_cai(6, :) * (7.0 / 6.2) * 372;
    %i_tot = i_tot_hp2x7 + i_tot_hp2x4;
	hold on;

    %tint = tint / 60.0;
    gg = 1;
    alpha = 0;%- gg * 45 + 45;
    cai1x = cai1 * gg + alpha;
    cai2x = cai2 * gg + alpha;
    cai3x = cai3 * gg + alpha;
    cai4x = cai4 * gg + alpha;
    cai5x = cai5 * gg + alpha;
    plot(tint, cai1x, '-ko', 'MarkerIndices', 1:steps:length(cai1));
    plot(tint, cai2x, 'g--', 'MarkerIndices', 1:steps:length(cai2));
    plot(tint, cai3x, '-kx', 'MarkerIndices', 1:steps:length(cai3));
    plot(tint, cai4x, '-.r*', 'MarkerIndices', 1:steps:length(cai4));
    plot(tint, cai5x, 'b--*', 'MarkerIndices', 1:steps:length(cai5));
    mm = max(cai5x) - baseline;
	plot([ATP_t1 ATP_t2], [baseline + mm/2 baseline + mm/2], '-k' , 'LineWidth', 2);
	text(ATP_t1 + (ATP_t2 - ATP_t1)/2, baseline + mm/2 + 0.05 * mm/2, 'ATP');
	hold off;
	%text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 51.5, 'ATP');
    ATP = ATP * 1e3;
    legend(sprintf('[Ca^2+_i] at ATP=%gmM', ATP(1)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(2)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(3)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(4)), sprintf('[Ca^2+_i] at ATP=%gmM', ATP(5)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    %ylabel('\Delta[Ca^2+_i] (\muM)');
    ylabel('[Ca^2+_i] (nM)');
    %ylabel('[Ca^2+_i] (M)');
    %ylabel('[Q] (M)');
    xlim([0 t_max]);
    %ylim([44 90]);
    %ylim([44 max( [max(cai1) max(cai2) max(cai3) max(cai4) max(cai5)]) ]);%
    
    %return
    %figure
    %plot(tint, caii, 'g--', 'MarkerIndices', 1:40:length(caii));
    
    %f = figure;
    %f.WindowStyle = 'docked';
    %plot(tint, cai3, '-.r*', 'MarkerIndices', 1:40:length(cai2));
    
    %return;
    
    fig2 = figure('Name', 'Simulation of intracellular/cytosolic [NAi] concentration for P2X4/7 receptor in human microglia');
	fig2.WindowStyle = 'docked';
    nai1 = y_p2x_cai1(9, :) * 1e3;
    nai2 = y_p2x_cai2(9, :) * 1e3;
    nai3 = y_p2x_cai3(9, :) * 1e3;
    nai4 = y_p2x_cai4(9, :) * 1e3;
    nai5 = y_p2x_cai5(9, :) * 1e3;
	hold on;
    plot(tint, nai1, '-ko', 'MarkerIndices', 1:steps:length(nai1));
    plot(tint, nai2, '-.r*', 'MarkerIndices', 1:steps:length(nai2));
    plot(tint, nai3, '-kx', 'MarkerIndices', 1:steps:length(nai3));
    plot(tint, nai4, 'g--', 'MarkerIndices', 1:steps:length(nai4));
    plot(tint, nai5, 'b--*', 'MarkerIndices', 1:steps:length(nai5));
	plot([ATP_t1/1 ATP_t2/1], [8.001 8.001], '-k' , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 8.0018, 'ATP');
	hold off;
    legend(sprintf('[Na^+_i] at ATP=%gmM', ATP(1)), sprintf('[Na^+_i] at ATP=%gmM', ATP(2)), sprintf('[Na^+_i] at ATP=%gmM', ATP(3)), sprintf('[Na^+_i] at ATP=%gmM', ATP(4)), sprintf('[Na^+_i] at ATP=%dmM', ATP(5)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('[Na^+_i] (mM)');
    xlim([0 t_max]);
    ylim([8 Inf]);
    %ylim([7.9994 max( [max(nai1) max(nai2) max(nai3) max(nai4) max(nai5)]) ]);%
    
    fig3 = figure('Name', 'Simulation of the membrane potential induced by P2X4/7 receptor in human microglia');
	fig3.WindowStyle = 'docked';
    Vm1 = y_p2x_cai1(11, :) * 1e3;
    Vm2 = y_p2x_cai2(11, :) * 1e3;
    Vm3 = y_p2x_cai3(11, :) * 1e3;
    Vm4 = y_p2x_cai4(11, :) * 1e3;
    Vm5 = y_p2x_cai5(11, :) * 1e3;
	hold on;
    plot(tint, Vm1, '-ko', 'MarkerIndices', 1:steps:length(Vm1));
    plot(tint, Vm2, '-.r*', 'MarkerIndices', 1:steps:length(Vm2));
    plot(tint, Vm3, '-kx', 'MarkerIndices', 1:steps:length(Vm3));
    plot(tint, Vm4, 'g--', 'MarkerIndices', 1:steps:length(Vm4));
    plot(tint, Vm5, 'b--*', 'MarkerIndices', 1:steps:length(Vm5));
	plot([ATP_t1/1 ATP_t2/1], [-56 -56], '-k' , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), -55, 'ATP');
	hold off;
    legend(sprintf('Vm at ATP=%gmM', ATP(1)), sprintf('Vm at ATP=%gmM', ATP(2)), sprintf('Vm at ATP=%gmM', ATP(3)), sprintf('Vm at ATP=%gmM', ATP(4)), sprintf('Vm at ATP=%gmM', ATP(5)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('Vm (mV)');
    xlim([0 t_max]);
    ylim([-61 max( [max(Vm1) max(Vm2) max(Vm3) max(Vm4) max(Vm5)]) ]);%
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    cai1 = cai1 * 1e-9;
    cai2 = cai2 * 1e-9;
    cai3 = cai3 * 1e-9;
    cai4 = cai4 * 1e-9;
    cai5 = cai5 * 1e-9;
    
    nai1 = nai1 * 1e-3;
    nai2 = nai2 * 1e-3;
    nai3 = nai3 * 1e-3;
    nai4 = nai4 * 1e-3;
    nai5 = nai5 * 1e-3;
    
    Vm1 = Vm1 * 1e-3;
    Vm2 = Vm2 * 1e-3;
    Vm3 = Vm3 * 1e-3;
    Vm4 = Vm4 * 1e-3;
    Vm5 = Vm5 * 1e-3;
    
    O11 = y_p2x_cai1(4, :);
    O12 = y_p2x_cai2(4, :);
    O13 = y_p2x_cai3(4, :);
    O14 = y_p2x_cai4(4, :);
    O15 = y_p2x_cai5(4, :);
    
    O21 = y_p2x_cai1(8, :);
    O22 = y_p2x_cai2(8, :);
    O23 = y_p2x_cai3(8, :);
    O24 = y_p2x_cai4(8, :);
    O25 = y_p2x_cai5(8, :);
    
    R = 8.314; % Ideal gas constant; unit: I/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    NAx = 130e-3;   % M
    CAx = 2e-3;   % M
    Z_Na = 1;
    Z_Ca = 2;
    d = 25;
    s_p = 243.53 * 1e-12; % the area of processes
    s_b = 74.3 * 1e-12; % the area of cell body
    s_tot = s_p + s_b;
    
    for i=1:1:length(cai1)
        I_ncx1(i) = NA_NCX_Current(nai1(i), cai1(i), Vm1(i), opt);
        I_ncx2(i) = NA_NCX_Current(nai2(i), cai2(i), Vm2(i), opt);
        I_ncx3(i) = NA_NCX_Current(nai3(i), cai3(i), Vm3(i), opt);
        I_ncx4(i) = NA_NCX_Current(nai4(i), cai4(i), Vm4(i), opt);
        I_ncx5(i) = NA_NCX_Current(nai5(i), cai5(i), Vm5(i), opt);
        
        I_pmca1(i) = PMCA_Current(cai1(i), opt);
        I_pmca2(i) = PMCA_Current(cai2(i), opt);
        I_pmca3(i) = PMCA_Current(cai3(i), opt);
        I_pmca4(i) = PMCA_Current(cai4(i), opt);
        I_pmca5(i) = PMCA_Current(cai5(i), opt);
        
        I_leak_na1(i) = NA_Leak_Current(nai1(i), Vm1(i), opt);
        I_leak_na2(i) = NA_Leak_Current(nai2(i), Vm2(i), opt);
        I_leak_na3(i) = NA_Leak_Current(nai3(i), Vm3(i), opt);
        I_leak_na4(i) = NA_Leak_Current(nai4(i), Vm4(i), opt);
        I_leak_na5(i) = NA_Leak_Current(nai5(i), Vm5(i), opt);
        
        I_leak_ca1(i) = CA_Leak_Current(cai1(i), Vm1(i), opt);
        I_leak_ca2(i) = CA_Leak_Current(cai2(i), Vm2(i), opt);
        I_leak_ca3(i) = CA_Leak_Current(cai3(i), Vm3(i), opt);
        I_leak_ca4(i) = CA_Leak_Current(cai4(i), Vm4(i), opt);
        I_leak_ca5(i) = CA_Leak_Current(cai5(i), Vm5(i), opt);
        
        E_Na1(i) =  R * T * Z_Na^-1 * F^-1 * log(NAx/nai1(i));
        E_Na2(i) =  R * T * Z_Na^-1 * F^-1 * log(NAx/nai2(i));
        E_Na3(i) =  R * T * Z_Na^-1 * F^-1 * log(NAx/nai3(i));
        E_Na4(i) =  R * T * Z_Na^-1 * F^-1 * log(NAx/nai4(i));
        E_Na5(i) =  R * T * Z_Na^-1 * F^-1 * log(NAx/nai5(i));

        E_Ca1(i) =  R * T * Z_Ca^-1 * F^-1 * log(CAx/cai1(i));
        E_Ca2(i) =  R * T * Z_Ca^-1 * F^-1 * log(CAx/cai2(i));
        E_Ca3(i) =  R * T * Z_Ca^-1 * F^-1 * log(CAx/cai3(i));
        E_Ca4(i) =  R * T * Z_Ca^-1 * F^-1 * log(CAx/cai4(i));
        E_Ca5(i) =  R * T * Z_Ca^-1 * F^-1 * log(CAx/cai5(i));
        
        I_hP2X7_Na1(i) = 0.12 * 2.417 * O11(i) * (Vm1(i) - E_Na1(i)) * 1e3 * 1e-12 / s_tot/d;    % 0.37
        I_hP2X7_Na2(i) = 0.12 * 2.417 * O12(i) * (Vm2(i) - E_Na2(i)) * 1e3 * 1e-12 / s_tot/d;    % 0.37
        I_hP2X7_Na3(i) = 0.12 * 2.417 * O13(i) * (Vm3(i) - E_Na3(i)) * 1e3 * 1e-12 / s_tot/d;    % 0.37
        I_hP2X7_Na4(i) = 0.12 * 2.417 * O14(i) * (Vm4(i) - E_Na4(i)) * 1e3 * 1e-12 / s_tot/d;    % 0.37
        I_hP2X7_Na5(i) = 0.12 * 2.417 * O15(i) * (Vm5(i) - E_Na5(i)) * 1e3 * 1e-12 / s_tot/d;    % 0.37
        
        I_hP2X7_Ca1(i) = 0.1 * 2.417 * O11(i) * (Vm1(i) - E_Ca1(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hP2X7_Ca2(i) = 0.1 * 2.417 * O12(i) * (Vm2(i) - E_Ca2(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hP2X7_Ca3(i) = 0.1 * 2.417 * O13(i) * (Vm3(i) - E_Ca3(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hP2X7_Ca4(i) = 0.1 * 2.417 * O14(i) * (Vm4(i) - E_Ca4(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hP2X7_Ca5(i) = 0.1 * 2.417 * O15(i) * (Vm5(i) - E_Ca5(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
    
        I_hp2x4_Na1(i) = 0.12 * 6.2 * O21(i) * (Vm1(i) - E_Na1(i)) * 1e3 * 1e-12/ s_tot/d;    % 0.38
        I_hp2x4_Na2(i) = 0.12 * 6.2 * O22(i) * (Vm2(i) - E_Na2(i)) * 1e3 * 1e-12/ s_tot/d;    % 0.38
        I_hp2x4_Na3(i) = 0.12 * 6.2 * O23(i) * (Vm3(i) - E_Na3(i)) * 1e3 * 1e-12/ s_tot/d;    % 0.38
        I_hp2x4_Na4(i) = 0.12 * 6.2 * O24(i) * (Vm4(i) - E_Na4(i)) * 1e3 * 1e-12/ s_tot/d;    % 0.38
        I_hp2x4_Na5(i) = 0.12 * 6.2 * O25(i) * (Vm5(i) - E_Na5(i)) * 1e3 * 1e-12/ s_tot/d;    % 0.38
        
        I_hp2x4_Ca1(i) = 0.1 * 6.2 * O21(i) * (Vm1(i) - E_Ca1(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hp2x4_Ca2(i) = 0.1 * 6.2 * O22(i) * (Vm2(i) - E_Ca2(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hp2x4_Ca3(i) = 0.1 * 6.2 * O23(i) * (Vm3(i) - E_Ca3(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hp2x4_Ca4(i) = 0.1 * 6.2 * O24(i) * (Vm4(i) - E_Ca4(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
        I_hp2x4_Ca5(i) = 0.1 * 6.2 * O25(i) * (Vm5(i) - E_Ca5(i)) * 1e3 * 1e-12/ s_tot/d;      % 0.1
    end
    
    %return
    fig4 = figure('Name', 'Simulation of NCX current density induced by P2X4/7 receptor in human microglia');
	fig4.WindowStyle = 'docked';
	hold on;
    plot(tint, I_ncx1, '-ko', 'MarkerIndices', 1:steps:length(I_ncx1));
    plot(tint, I_ncx2, '-.r*', 'MarkerIndices', 1:steps:length(I_ncx2));
    plot(tint, I_ncx3, '-kx', 'MarkerIndices', 1:steps:length(I_ncx3));
    plot(tint, I_ncx4, 'g--', 'MarkerIndices', 1:steps:length(I_ncx4));
    plot(tint, I_ncx5, 'b--*', 'MarkerIndices', 1:steps:length(I_ncx5));
	plot([ATP_t1/1 ATP_t2/1], [0.009 0.009], '-k' , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.0092, 'ATP');
	hold off;
    legend(sprintf('I_{NCX} at ATP=%gmM', ATP(1)), sprintf('I_{NCX} at ATP=%gmM', ATP(2)), sprintf('I_{NCX} at ATP=%gmM', ATP(3)), sprintf('I_{NCX} at ATP=%gmM', ATP(4)), sprintf('I_{NCX} at ATP=%gmM', ATP(5)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('I_{NCX} (A/m^2)');
    xlim([0 t_max]);
    %ylim([0.007 max( [max(I_ncx1) max(I_ncx2) max(I_ncx3) max(I_ncx4) max(I_ncx5)]) ]);
   
    
    % Currents
    fig5 = figure('Name', 'Simulation of sodium/calcium currents induced by P2X4/7 receptors in human microglia');
	fig5.WindowStyle = 'docked';
    subplot(2, 3, 1);
	hold on;
    plot(tint, I_ncx2 * s_p * 1e12, '-ko', 'MarkerIndices', 1:steps:length(I_ncx2));
    plot(tint, I_leak_na2 * s_p * 1e12, '-kx', 'MarkerIndices', 1:steps:length(I_leak_na2));
    plot(tint, I_hP2X7_Na2 * s_tot * 1e12, 'b--*', 'MarkerIndices', 1:steps:length(I_hP2X7_Na2));
    plot(tint, I_hp2x4_Na2 * s_tot * 1e12, '-.r*', 'MarkerIndices', 1:steps:length(I_hp2x4_Na2));
	%%plot([ATP_t1/1 ATP_t2/1], [1e-12 * 1e12 1e-12 * 1e12], '-k');% , 'LineWidth', 2);
    %%text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 1.2e-12 * 1e12, sprintf('ATP=%gmM', ATP(2)), 'FontSize', 8);
    plot([ATP_t1/1 ATP_t2/1], [0.5 0.5], '-k');% , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.7, sprintf('ATP=%gmM', ATP(2)), 'FontSize', 8);
	hold off;
    legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X7}_{Na}', 'I_{hP2X4}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X7}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X4}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('I (pA)');
    xlim([0 t_max]);
    %ylim([-4 Inf]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    %fig6 = figure('Name', sprintf('Simulation of sodium currents induced by P2X4/7 receptors in human microglia at ATP=%gmM',  ATP(3)));
	%fig6.WindowStyle = 'docked';
    subplot(2, 3, 2);
	hold on;
    plot(tint, I_ncx4 * s_p * 1e12, '-ko', 'MarkerIndices', 1:steps:length(I_ncx4));
    plot(tint, I_leak_na4 * s_p * 1e12, '-kx', 'MarkerIndices', 1:steps:length(I_leak_na4));
    plot(tint, I_hP2X7_Na4 * s_tot * 1e12, 'b--*', 'MarkerIndices', 1:steps:length(I_hP2X7_Na4));
    plot(tint, I_hp2x4_Na4 * s_tot * 1e12, '-.r*', 'MarkerIndices', 1:steps:length(I_hp2x4_Na4));
    plot([ATP_t1/1 ATP_t2/1], [0.5 0.5], '-k');% , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.7, sprintf('ATP=%gmM', ATP(4)), 'FontSize', 8);
	hold off;
    legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X7}_{Na}', 'I_{hP2X4}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X7}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X4}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('I (pA)');
    xlim([0 t_max]);
    %ylim([-4 Inf]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    %fig7 = figure('Name', sprintf('Simulation of sodium currents induced by P2X4/7 receptors in human microglia at ATP=%gmM',  ATP(4)));
	%fig7.WindowStyle = 'docked';
    subplot(2, 3, 3);
	hold on;
    plot(tint, I_ncx5 * s_p * 1e12, '-ko', 'MarkerIndices', 1:steps:length(I_ncx5));
    plot(tint, I_leak_na5 * s_p * 1e12, '-kx', 'MarkerIndices', 1:steps:length(I_leak_na5));
    plot(tint, I_hP2X7_Na5 * s_tot * 1e12, 'b--*', 'MarkerIndices', 1:steps:length(I_hP2X7_Na5));
    plot(tint, I_hp2x4_Na5 * s_tot * 1e12, '-.r*', 'MarkerIndices', 1:steps:length(I_hp2x4_Na5));
    plot([ATP_t1/1 ATP_t2/1], [0.5 0.5], '-k');% , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.7, sprintf('ATP=%gmM', ATP(5)), 'FontSize', 8);
	hold off;
    legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X7}_{Na}', 'I_{hP2X4}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X7}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{NaNCX}', 'I_{NaLeak}', 'I_{hP2X4}_{Na}', 'FontSize', 8, 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('I (pA)');
    xlim([0 t_max]);
    %ylim([-4 Inf]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    %return
    %fig8 = figure('Name', sprintf('Simulation of calcium currents induced by P2X4/7 receptors in human microglia at ATP=%gmM',  ATP(2)));
	%fig8.WindowStyle = 'docked';
    subplot(2, 3, 4);
	hold on;
    plot(tint, (-2/3) * I_ncx2 * s_p * 1e12, '-ko', 'MarkerIndices', 1:steps:length(I_ncx2));
    plot(tint, I_pmca2 * s_p * 1e12, 'g--', 'MarkerIndices', 1:steps:length(I_pmca2));
    plot(tint, I_leak_ca2 * s_p * 1e12, '-kx', 'MarkerIndices', 1:steps:length(I_leak_ca2));
    plot(tint, I_hP2X7_Ca2 * s_tot * 1e12, 'b--*', 'MarkerIndices', 1:steps:length(I_hP2X7_Ca2));
    plot(tint, I_hp2x4_Ca2 * s_tot * 1e12, '-.r*', 'MarkerIndices', 1:steps:length(I_hp2x4_Ca2));
    plot([ATP_t1/1 ATP_t2/1], [0.5 0.5], '-k');% , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.8, sprintf('ATP=%gmM', ATP(2)), 'FontSize', 8);
	hold off;
    legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X7}_{Ca}', 'I_{hP2X4}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X7}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X4}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('I (pA)');
    xlim([0 t_max]);
    %ylim([-4 Inf]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    %fig9 = figure('Name', sprintf('Simulation of calcium currents induced by P2X4/7 receptors in human microglia at ATP=%gmM',  ATP(3)));
	%fig9.WindowStyle = 'docked';
    subplot(2, 3, 5);
	hold on;
    plot(tint, (-2/3) * I_ncx4 * s_p * 1e12, '-ko', 'MarkerIndices', 1:steps:length(I_ncx4));
    plot(tint, I_pmca4 * s_p * 1e12, 'g--', 'MarkerIndices', 1:steps:length(I_pmca4));
    plot(tint, I_leak_ca4 * s_p * 1e12, '-kx', 'MarkerIndices', 1:steps:length(I_leak_ca4));
    plot(tint, I_hP2X7_Ca4 * s_tot * 1e12, 'b--*', 'MarkerIndices', 1:steps:length(I_hP2X7_Ca4));
    plot(tint, I_hp2x4_Ca4 * s_tot * 1e12, '-.r*', 'MarkerIndices', 1:steps:length(I_hp2x4_Ca4));
    plot([ATP_t1/1 ATP_t2/1], [0.5 0.5], '-k');% , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.8, sprintf('ATP=%gmM', ATP(4)), 'FontSize', 8);
	hold off;
    legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X7}_{Ca}', 'I_{hP2X4}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X7}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X4}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('I (pA)');
    xlim([0 t_max]);
    %ylim([-4 Inf]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    %fig10 = figure('Name', sprintf('Simulation of calcium currents induced by P2X4/7 receptors in human microglia at ATP=%gmM',  ATP(4)));
	%fig10.WindowStyle = 'docked';
    subplot(2, 3, 6);
	hold on;
    plot(tint, (-2/3) * I_ncx5 * s_p * 1e12, '-ko', 'MarkerIndices', 1:steps:length(I_ncx5));
    plot(tint, I_pmca5 * s_p * 1e12, 'g--', 'MarkerIndices', 1:steps:length(I_pmca5));
    plot(tint, I_leak_ca5 * s_p * 1e12, '-kx', 'MarkerIndices', 1:steps:length(I_leak_ca5));
    plot(tint, I_hP2X7_Ca5 * s_tot * 1e12, 'b--*', 'MarkerIndices', 1:steps:length(I_hP2X7_Ca5));
    plot(tint, I_hp2x4_Ca5 * s_tot * 1e12, '-.r*', 'MarkerIndices', 1:steps:length(I_hp2x4_Ca5));
    plot([ATP_t1/1 ATP_t2/1], [0.5 0.5], '-k');% , 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.8, sprintf('ATP=%gmM', ATP(5)), 'FontSize', 8);
	hold off;
    legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X7}_{Ca}', 'I_{hP2X4}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X7}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    %legend('I_{CaNCX}', 'I_{PMCA}', 'I_{CaLeak}', 'I_{hP2X4}_{Ca}', 'FontSize', 8, 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('I (pA)');
    xlim([0 t_max]);
    %ylim([-4 Inf]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    % Reversal potentials    
    fig6 = figure('Name', sprintf('Simulation of sodium reversal potentials induced by P2X4/7 receptors in human microglia at ATP=%gmM',  ATP(3)));
	fig6.WindowStyle = 'docked';
    %subplot(1, 3, 2);
    hold on;
    plot(tint, E_Na2, '-ko', 'MarkerIndices', 1:steps:length(E_Na2));
    plot(tint, E_Na4, '-.r*', 'MarkerIndices', 1:steps:length(E_Na4));
    plot(tint, E_Na5, 'b--*', 'MarkerIndices', 1:steps:length(E_Na5));
	plot([ATP_t1/1 ATP_t2/1], [0.07446 0.07446], '-k', 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.074463, 'ATP');
	hold off;
    legend(sprintf('E_{Na} at ATP=%gmM', ATP(2)), sprintf('E_{Na} at ATP=%gmM', ATP(4)), sprintf('E_{Na} at ATP=%gmM', ATP(5)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('E_{Na} (mV)');
    xlim([0 t_max]);
    %ylim([-Inf 6]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    fig7 = figure('Name', 'Simulation of calcium reversal potentials induced by P2X4/7 receptors in human microglia');
	fig7.WindowStyle = 'docked';
    %subplot(1, 3, 1);
	hold on;
    plot(tint, E_Ca2, '-ko', 'MarkerIndices', 1:steps:length(E_Ca2));
    plot(tint, E_Ca4, '-.r*', 'MarkerIndices', 1:steps:length(E_Ca4));
    plot(tint, E_Ca5, 'b--*', 'MarkerIndices', 1:steps:length(E_Ca5));
	plot([ATP_t1/1 ATP_t2/1], [0.141 0.141], '-k', 'LineWidth', 2);
    text((ATP_t1 + (ATP_t2 - ATP_t1)/2), 0.1413, 'ATP');
	hold off;
    legend(sprintf('E_{Ca} at ATP=%gmM', ATP(2)), sprintf('E_{Ca} at ATP=%gmM', ATP(4)), sprintf('E_{Ca} at ATP=%gmM', ATP(5)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('E_{Ca} (mV)');
    xlim([0 t_max]);
    %ylim([-Inf 6]);
    %ylim([-0.023 max( [max(I_leak_na1) max(I_leak_na2) max(I_leak_na3) max(I_leak_na4) max(I_leak_na5)]) ]);
    
    opt.ATP = 1e-3;
    opt.test = false;
    sol_p2x_cai(1) = ode23tb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
    opt.test = true;
    sol_p2x_cai(2) = ode23tb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt)), opt.tspan, opt.x0, odeopt);
    opt.test = false;
    y_p2x_cai1 = deval(sol_p2x_cai(1), tint);
    y_p2x_cai2 = deval(sol_p2x_cai(2), tint);
    
    fig8 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2X4/7 receptor in human microglia by varying the volume');
	fig8.WindowStyle = 'docked';
    cai1 = y_p2x_cai1(10, :) * 1e9;
    cai2 = y_p2x_cai2(10, :) * 1e9;
	hold on;

    %tint = normalize(tint, 'range');
    %cai1 = normalize(cai1, 'range');
    %cai2 = normalize(cai2, 'range');
    plot(tint, cai1, '-ko', 'MarkerIndices', 1:steps:length(cai1));
    plot(tint, cai2, '-.r*', 'MarkerIndices', 1:steps:length(cai4));
    mm = max(cai5x) - 45;
	plot([ATP_t1 ATP_t2], [45 + mm/2 45 + mm/2], '-k' , 'LineWidth', 2);
	text(ATP_t1 + (ATP_t2 - ATP_t1)/2, 45 + mm/2 + 0.05 * mm/2, 'ATP');
	hold off;
    %legend('Cubic volume', 'Spherical volume', 'FontWeight', 'bold');
    legend('Original shape', '10% perturbation of the original shape', 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('[Ca^2+_i] (nM)');
    %xlim([0 t_max]);
    abs(max(cai1) - max(cai2))
end
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

    g_NA_Leak = (-I_NA_NCX) / (Vm - Ena);
    g_CA_Leak = (-I_CA_NCX - PMCA_Current(CAi, opt)) /  (Vm - Eca);
end
%--------------------%