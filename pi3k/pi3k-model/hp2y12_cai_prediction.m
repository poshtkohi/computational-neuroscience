%------Functions -------%
%function hp2y12_cai_prediction()
    close all; clc; clear; clearvars;
    format long g;
    %fclose('all');
	c = Constants; % Gets constants for the model
    c.should_use_matlab_solver = 1;
   
    ADP = 50; %uM
    %ADP_t1 = 2.0; ADP_t2 = 4.6; t_max = 40.0;% unit in seconds - Original
    %offset = 0; ADP_t1 = offset + 13.88; ADP_t2 = offset + 13.88 + 29.16; t_max = offset + 66.0;
    
    offset = 0; ADP_t1 = offset + 13.88; ADP_t2 = offset + 13.88 + 1200; t_max = offset + 2000;
    
    offset = 20;
    ADP_t1 = 0 + offset; ADP_t2 = 30 * 60 + offset; t_max = 9000;%t_max = 30 * 60 + offset;% unit in seconds
    
	offset = 20;
    %%%ADP_t1 = 0 + offset; ADP_t2 = 30 * 60 + offset; t_max = 2000;
    %ADP_t1 = 2.0; ADP_t2 = 41; t_max = 60.0;% unit in seconds
    %ADP_t1 = 10; ADP_t2 = 130.0; t_max = 150;% unit in seconds
    %ADP_t1 = 80.0; ADP_t2 = 480; t_max = 540.0;% unit in seconds
    %ADP_t1 = 30.0; ADP_t2 = (116.0 + 30.0); t_max = 180;% unit in seconds
    %ADP_t1 = 80.0; ADP_t2 = 480; t_max = 540.0;% unit in seconds
    
    %ATP_t1 = 2; ATP_t2 = 17; t_max = 200;% unit in seconds
    %----------------------
    global t1 t2 g enabled;
    g = 0.5;
    enabled = 1;
    t1 = ADP_t1;
    t2 = ADP_t2;
    t_max = t_max;
    opt_p2y_cai.delay = c.delay * 60;
    opt_p2y_cai.width = t2 - t1;
    opt_p2y_cai.amplitude = ADP;
    steps = 500;
    %----------------------
   
    es_hp2y12_cai = load(sprintf('variablescmaes-hp2y12-ADP-50-popsize-%d.mat', c.PopSize), '-mat');
	%K_hp2y12_cai = es_hp2y12_cai.out.solutions.bestever.x;
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
	opt_p2y_cai.ADP_t1 = ADP_t1;
	opt_p2y_cai.ADP_t2 = ADP_t2;
	opt_p2y_cai.tspan = [0.0 t_max];
    opt_p2y_cai.observable_index = 8;
    
    global my_k my_ADP;
    %my_k = k_hp2y12_i_tot;
    odeopt = odeset;%('InitialStep', 1e-2, 'MaxStep', 1e-2);%('RelTol', 1e-9, 'AbsTol', 1e-9);
    %odeopt = odeset;%('Jacobian', @compute_jacobian_matrix_hp2y12_total_current);
    opt_p2y_cai.ADP = ADP;
    sol_p2y_cai = ode15s(@(t, x)(reaction_network_hp2y12(t, x, opt_p2y_cai)), opt_p2y_cai.tspan, opt_p2y_cai.x0, odeopt);

    tint = linspace(0, t_max, 10000);  % all time points.
    y_p2y_cai1 = deval(sol_p2y_cai, tint);

	fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2Y12 receptor in human microglia');
	fig1.WindowStyle = 'docked'; %
    cai1 = y_p2y_cai1(opt_p2y_cai.observable_index, :);
    
	hold on;
    if c.periodic == 1
        g = 0.5;
        enabled = 1;
        t1 = ADP_t1;
        t2 = ADP_t2;
        A = zeros(length(tint), 1);
        G = zeros(length(tint), 1);
        for i=1:1:length(tint)
            [A(i), G(i)] = periodic_function(tint(i), ADP_t1, ADP_t2, opt_p2y_cai.delay, 0.5);
        end
        plot(tint, normalize(cai1, 'range'), '-ko', 'MarkerIndices', 1:steps:length(cai1));
        %plot(tint, cai1, '-ko', 'MarkerIndices', 1:steps:length(cai1));
        plot(tint, A, 'b--*', 'MarkerIndices', 1:steps:length(cai1));
        legend(sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP), 'ADP', 'FontWeight', 'bold');
    else
        plot(tint, cai1, '-ko', 'MarkerIndices', 1:steps:length(cai1));
        plot([ADP_t1/1 ADP_t2/1], [0.2 0.2], '-k' , 'LineWidth', 2);
        text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.22, 'ADP');
        legend(sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP), 'FontWeight', 'bold');
    end
    hold off;
    %ADP = ADP * 1e3;
    %legend(sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(1)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(2)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(3)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(4)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(5)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('[Ca^2+_i] (\muM)');
    %ylabel('\Delta[Ca^2+_i] (\muM)');
    %xlim([0 t_max]);
    %ylim([44 max( [max(cai1) max(cai2) max(cai3) max(cai4) max(cai5)]) ]);%
    %return
%end
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