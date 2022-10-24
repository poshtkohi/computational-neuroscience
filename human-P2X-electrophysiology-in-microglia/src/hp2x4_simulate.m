%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function hp2x4_simulate()
    close all; clc; clear; clearvars;
    format long g;
    ATP = [0.01e-3 0.1e-3 0.5e-3 1e-3 2e-3 3e-3];
    %ATP = [0.5e-3 0.5e-3 0.5e-3 0.5e-3 0.5e-3 0.5e-3];
    %ATP = [0.01e-3 0.01e-3 0.01e-3 0.01e-3 0.01e-3 0.01e-3];
    %ATP = 0.1e-3 * ones(5, 1);
    %ATP = [0.01e-3 0.3e-3 0.7e-3 3e-3 10e-3];
    %ATP = [0 0 0 0 0]; 
    ATP_t1 = 2; ATP_t2 = 4.6; t_max = 8.5;% unit in seconds - Original
    %ATP_t1 = 1; ATP_t2 = 1.6; t_max = 12;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 32; t_max = 40;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 122; t_max = 130;% unit in seconds
    %ATP_t1 = 80.0; ATP_t2 = 480; t_max = 540.0;% unit in seconds
    %ATP_t1 = 20; ATP_t2 = 350; t_max = 360;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 7; t_max = 20;% unit in seconds
    %ATP_t1 = 20; ATP_t2 = 80.0; t_max = 100;% unit in seconds
    %ATP_t1 = 19.56; ATP_t2 = 19.56 + 30.43; t_max = 70.0;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 30; t_max = 50;% unit in seconds
    %ATP_t1 = 13.46; ATP_t2 = 13.46 + 12.5; t_max = 51;% unit in seconds
    % ATP_t1 = 15.68; ATP_t2 = 25.255; t_max = 53;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 17; t_max = 50;% unit in seconds
    ATP_t1 = 2; ATP_t2 = 170; t_max = 215;% unit in seconds
    
    global nT_P2X T;
    T = 0; % Period in seconds
    nT_P2X = 0;

    es_hp2x4_i_tot = load('variablescmaes-hp2x4-ATP-1.000000e-04-popsize-18-total-current.mat');
	K_hp2x4_i_tot = es_hp2x4_i_tot.out.solutions.bestever.x;
    %K_hp2x4_i_tot = [1.523384e+05 4.679995e+03 4.149279e+02 6.874058e-01 1.942445e+02 5.770408e+01 7.289234e-01 3.116321e+03 1.598147e+00 2.494154e+07 ];
    K_hp2x4_i_tot = read_last_line();
	k_hp2x4_i_tot = K_hp2x4_i_tot;%(1:7);
    k = k_hp2x4_i_tot ;
	opt_hp2x4_i_tot.x0 = zeros(4, 1);
    opt_hp2x4_i_tot.x0(1, 1) = 1;
	opt_hp2x4_i_tot.initial_condition_index = length(k_hp2x4_i_tot) + 1;
	opt_hp2x4_i_tot.x0 = load_initial_conditions_from_k_hp2x4_total_current(K_hp2x4_i_tot, opt_hp2x4_i_tot);
	%opt_hp2x4_i_tot.x0
    %opt_hp2x4_i_tot.x0(4, 1) = 0;
	%K_hp2x4_i_tot
	opt_hp2x4_i_tot.K = k_hp2x4_i_tot;
	opt_hp2x4_i_tot.ATP = ATP;
	opt_hp2x4_i_tot.ATP_t1 = ATP_t1;
	opt_hp2x4_i_tot.ATP_t2 = ATP_t2;
	opt_hp2x4_i_tot.tspan = [0.0 t_max];
    observable_index = 4;
    %g = 2.416;
    a = - 1;%6.8 / 6.2) * 372;
    g = k(22);
    k(4:5)

	global my_k my_ATP;
	my_k = k_hp2x4_i_tot;
	odeopt = odeset;%('RelTol', 1e-12, 'AbsTol', 1e-12);
	%odeopt = odeset;%('Jacobian', @compute_jacobian_matrix_hp2x4_total_current);
    for i=1:1:length(ATP)
        opt_hp2x4_i_tot.ATP = ATP(i);
        my_ATP = ATP(i);
        sol_hp2x4_i_tot(i) = ode15s(@(t, x)(reaction_network_hp2x4_total_current(t, x, opt_hp2x4_i_tot)), opt_hp2x4_i_tot.tspan, opt_hp2x4_i_tot.x0, odeopt);
    end

	tint = linspace(0, t_max, 1000);  % all time points.
	y_hp2x4_i_tot1 = deval(sol_hp2x4_i_tot(1), tint);
    y_hp2x4_i_tot2 = deval(sol_hp2x4_i_tot(2), tint);
    y_hp2x4_i_tot3 = deval(sol_hp2x4_i_tot(3), tint);
    y_hp2x4_i_tot4 = deval(sol_hp2x4_i_tot(4), tint);
    y_hp2x4_i_tot5 = deval(sol_hp2x4_i_tot(5), tint);
    y_hp2x4_i_tot6 = deval(sol_hp2x4_i_tot(6), tint);
    
    steps = 100;

	fig1 = figure('Name', 'Simulation of total current for P2X4 receptor in human microglia');
	fig1.WindowStyle = 'docked';
    
    for i=1:1:length(tint)
        h = 1 - y_hp2x4_i_tot1(3, i); 
        O = y_hp2x4_i_tot1(4, i);
        O1(i) = g * h^0 * O^1 * (60e-3 - 0);
    end
    
    for i=1:1:length(tint)
        h = 1 - y_hp2x4_i_tot2(3, i); 
        O = y_hp2x4_i_tot2(4, i);
        O2(i) = g * h^0 * O^1 * (60e-3 - 0);
    end
    
    for i=1:1:length(tint)
        h = 1 - y_hp2x4_i_tot3(3, i); 
        O = y_hp2x4_i_tot3(4, i);
        O3(i) = g * h^0 * O^1 * (60e-3 - 0);
    end
    
    for i=1:1:length(tint)
        h = 1 - y_hp2x4_i_tot4(3, i); 
        O = y_hp2x4_i_tot4(4, i);
        O4(i) = g * h^0 * O^1 * (60e-3 - 0);
    end
    
    for i=1:1:length(tint)
        h = 1 - y_hp2x4_i_tot5(3, i); 
        O = y_hp2x4_i_tot5(4, i);
        O5(i) = g * h^0 * O^1 * (60e-3 - 0);
    end
    
    for i=1:1:length(tint)
        h = 1 - y_hp2x4_i_tot6(3, i); 
        O = y_hp2x4_i_tot6(4, i);
        O6(i) = g * h^0 * O^1 * (60e-3 - 0);
    end
	%O1 = y_hp2x4_i_tot1(observable_index, :);
    %O2 = y_hp2x4_i_tot2(observable_index, :);
    %O3 = y_hp2x4_i_tot3(observable_index, :);
    %O4 = y_hp2x4_i_tot4(observable_index, :);
    %O5 = y_hp2x4_i_tot5(observable_index, :);
    %O6 = y_hp2x4_i_tot5(observable_index, :);
	hold on;
	%plot(tint, O1 * a, 'g--', 'MarkerIndices', 1:steps:length(O1));
    plot(tint, O2 * a, '-ko', 'MarkerIndices', 1:steps:length(O2));
    %plot(tint, O3 * a, '-kx', 'MarkerIndices', 1:steps:length(O3));
    plot(tint, O4 * a, '-.r*', 'MarkerIndices', 1:steps:length(O4));
    %plot(tint, O5 * a, 'b--*', 'MarkerIndices', 1:steps:length(O5));
    plot(tint, O6 * a, '--ko', 'MarkerIndices', 1:steps:length(O6));
    gg = max(-O6 * a);
	plot([ATP_t1 ATP_t2], [-gg/2 -gg/2], '-k' , 'LineWidth', 2);
	text(ATP_t1 + (ATP_t2 - ATP_t1)/2, -gg/2 + 0.05 * gg/2, 'ATP');
	hold off;
    ATP = ATP * 1e3;
    %legend(sprintf('I_{tot} at ATP=%gmM', ATP(1)), sprintf('I_{tot} at ATP=%gmM', ATP(2)), sprintf('I_{tot} at ATP=%gmM', ATP(3)), sprintf('I_{tot} at ATP=%gmM', ATP(4)), sprintf('I_{tot} at ATP=%gmM', ATP(5)), sprintf('I_{tot} at ATP=%gmM', ATP(6)), 'FontWeight', 'bold');
    legend(sprintf('I_{tot} at ATP=%gmM', ATP(2)), sprintf('I_{tot} at ATP=%gmM', ATP(4)), sprintf('I_{tot} at ATP=%gmM', ATP(6)), 'FontWeight', 'bold');
	xlabel('Time (s)');
	ylabel('I_{tot} (pA)');
    xlim([0 t_max]);
    ylim([-Inf 0]);
    %ylim([-406 0]);
    %k(25:27)
    return;
    fig = figure('Name', 'State space of the P2X4 receptor in rat microglia');
	fig.WindowStyle = 'docked';
    str = {'[C]', '[S]', '[D]', '[O]'};
    for i=1:1:length(opt_hp2x4_i_tot.x0)
       subplot(2, 2, i);
       %str = { sprintf('s%d', i) };
       plot(tint, y_hp2x4_i_tot5(i, :)); xlabel('Time (s)'); ylabel(str(i));
       %xlim([0.0 8.5]);
    end
    return
    
    %{
    % Write 5mM response of the hp2x4 to file
    file = fopen("hp2x4-ATP-5mM-total-current.dat", "w");
    for i=1:1:length(tint)
        fprintf(file, "%g\t%g\n", tint(i), O5(i) * a);
    end
    fclose(file);
    %}
    %return
    
    % Repeated application of ADP/ATP
    T = 30; % Period in seconds
    nT_P2X = 0;
    
    ATP = 0.01e-3; ATP_t1 = 0; ATP_t2 = 3; t_max = 500;
    
    opt_hp2x4_i_tot.ATP = ATP;
    opt_hp2x4_i_tot.ATP_t1 = ATP_t1;
    opt_hp2x4_i_tot.ATP_t2 = ATP_t2;
    opt_hp2x4_i_tot.tspan = [0.0 t_max];
    
    %return;
    my_k = k_hp2x4_i_tot;
    my_ATP = ATP;
    odeopt = odeset;%('Jacobian', @compute_jacobian_matrix_hp2x4_total_current);
    sol = ode15s(@(t, x)(reaction_network_hp2x4_total_current(t, x, opt_hp2x4_i_tot)), opt_hp2x4_i_tot.tspan, opt_hp2x4_i_tot.x0, odeopt);
    
    fig2 = figure('Name', 'Simulation for P2X4 receptor in rat microglia for repeated application of ATP');
	fig2.WindowStyle = 'docked';
    
    tint = linspace(0, t_max, 2000);  % all time points.
	y_hp2x4_i_tot = deval(sol, tint);
    
	O = y_hp2x4_i_tot(observable_index, :);
	hold on;
	plot(tint, O * a, '-ko', 'MarkerIndices', 1:50:length(O));
	%plot([ATP_t1 ATP_t2], [-50 -50], '-k' , 'LineWidth', 2);
	hold off;
	%text(ATP_t1 + (ATP_t2 - ATP_t1)/2, -45, 'ATP');
    ATP = ATP * 1e3;
    legend(sprintf('I_{tot} at at a pulse train of ATP=%gmM', ATP));
	xlabel('Time (s)');
	ylabel('I_{tot} (pA)');
    
    fig3 = figure('Name', 'State space of the P2X receptor in human microglia');
	fig3.WindowStyle = 'docked';
    for i=1:1:length(opt_hp2x4_i_tot.x0)
       subplot(4, 2, i);
       %str = { sprintf('s%d', i) };
       plot(tint, y_hp2x4_i_tot(i, :)); xlabel('Time (s)'); ylabel(str(i));
       %xlim([0.0 8.5]);
    end
end
%--------------------%
function [K] = read_last_line()
    filename = "outcmaes-hp2x4-ATP-1.000000e-04-popsize-18-total-currentxmean.dat";
    fid = fopen(filename, 'r');     %# Open the file as a binary
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