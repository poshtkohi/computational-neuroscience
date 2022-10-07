%------Functions -------%
function hp2x7_simulate()
    close all; clc; clear; clearvars;
    format long g;
    %ATP = [0.01e-3 0.1e-3 0.3e-3 1.0e-3 3e-3 5e-3];
    %ATP = [0.01e-3 0.1e-3 0.3e-3 0.4e-3 0.5e-3 0.7e-3];
    ATP = [1e-3 2e-3 3e-3 5e-3 7e-3 10e-3];
    %ATP = [0 0 0 0 0];
	%ATP_t1 = 2; ATP_t2 = 4.6; t_max = 8.5;% unit in seconds - Original
    ATP_t1 = 1; ATP_t2 = 1.5; t_max = 4;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 32; t_max = 35;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 122; t_max = 130;% unit in seconds
    %ATP_t1 = 80.0; ATP_t2 = 480; t_max = 540.0;% unit in seconds
    %ATP_t1 = 20; ATP_t2 = 60; t_max = 80;% unit in seconds00
    %ATP_t1 = 2; ATP_t2 = 80.0; t_max = 100;% unit in seconds
    %ATP_t1 = 2; ATP_t2 = 17; t_max = 20;% unit in seconds
    
    global nT_P2X T;
    T = 0; % Period in seconds
    nT_P2X = 0;

    %es_hp2x7_i_tot = load('variablescmaes-hp2x7-ATP-5.000000e-03-popsize-18-total-current.mat');
	%K_hp2x7_i_tot = es_hp2x7_i_tot.out.solutions.bestever.x;
    %K_hp2x7_i_tot = [2.001021e+03 9.992772e+02 6.506595e+00 2.299729e-03 3.357356e+00 3.019729e+00 1.998572e+00 9.929779e-02 9.577838e-01 1.490466e+00 1.986822e-03 1.598367e+00 1.684461e+00 2.652513e+00 1.123690e+00 1.771947e+00 1.311239e+00 1.528492e-01 1.635376e+00 1.249865e-01 1.167005e+00 1.773232e+00 8.447642e+00 2.046450e-01 ];
    K_hp2x7_i_tot = read_last_line();
	k_hp2x7_i_tot = K_hp2x7_i_tot;%K_hp2x7_i_tot(1:24);
	opt_hp2x7_i_tot.x0 = zeros(4, 1);
    opt_hp2x7_i_tot.x0(1, 1) = 1;
	opt_hp2x7_i_tot.initial_condition_index = length(k_hp2x7_i_tot) + 1;
	opt_hp2x7_i_tot.x0 = load_initial_conditions_from_k_hp2x7_total_current(K_hp2x7_i_tot, opt_hp2x7_i_tot);
	%opt_hp2x7_i_tot.x0
    %opt_hp2x7_i_tot.x0(4, 1) = 0;
	%K_hp2x7_i_tot
    k = k_hp2x7_i_tot;
	opt_hp2x7_i_tot.K = k_hp2x7_i_tot;
	opt_hp2x7_i_tot.ATP = ATP;
	opt_hp2x7_i_tot.ATP_t1 = ATP_t1;
	opt_hp2x7_i_tot.ATP_t2 = ATP_t2;
	opt_hp2x7_i_tot.tspan = [0.0 t_max];
    observable_index = 4;
    %g = 2.416;
    a = - 145;%(2.61 / 2.416) * 145;
    %a = -145;

	global my_k my_ATP;
	my_k = k_hp2x7_i_tot;
	odeopt = odeset;%('RelTol', 1e-12, 'AbsTol', 1e-12);
	%odeopt = odeset;%('Jacobian', @compute_jacobian_matrix_hp2x7_total_current);
    for i=1:1:length(ATP)
        opt_hp2x7_i_tot.ATP = ATP(i);
        my_ATP = ATP(i);
        sol_hp2x7_i_tot(i) = ode15s(@(t, x)(reaction_network_hp2x7_total_current(t, x, opt_hp2x7_i_tot)), opt_hp2x7_i_tot.tspan, opt_hp2x7_i_tot.x0, odeopt);
    end

	tint = linspace(0, t_max, 1000);  % all time points.
	y_hp2x7_i_tot1 = deval(sol_hp2x7_i_tot(1), tint);
    y_hp2x7_i_tot2 = deval(sol_hp2x7_i_tot(2), tint);
    y_hp2x7_i_tot3 = deval(sol_hp2x7_i_tot(3), tint);
    y_hp2x7_i_tot4 = deval(sol_hp2x7_i_tot(4), tint);
    y_hp2x7_i_tot5 = deval(sol_hp2x7_i_tot(5), tint);
    y_hp2x7_i_tot6 = deval(sol_hp2x7_i_tot(6), tint);

	fig1 = figure('Name', 'Simulation of total current for P2X7 receptor in human microglia');
	fig1.WindowStyle = 'docked';
	O1 = y_hp2x7_i_tot1(observable_index, :);
    O2 = y_hp2x7_i_tot2(observable_index, :);
    O3 = y_hp2x7_i_tot3(observable_index, :);
    O4 = y_hp2x7_i_tot4(observable_index, :);
    O5 = y_hp2x7_i_tot5(observable_index, :);
    O6 = y_hp2x7_i_tot6(observable_index, :);
	hold on;
	plot(tint, O1 * a, '-ko', 'MarkerIndices', 1:40:length(O1));
    plot(tint, O2 * a, 'g--', 'MarkerIndices', 1:40:length(O2));
    plot(tint, O3 * a, '-kx', 'MarkerIndices', 1:40:length(O3));
    plot(tint, O4 * a, '-.r*', 'MarkerIndices', 1:40:length(O4));
    plot(tint, O5 * a, 'b--*', 'MarkerIndices', 1:40:length(O5));
    plot(tint, O6 * a, '--ko', 'MarkerIndices', 1:40:length(O6));
    %plot(tint, y_hp2x7_i_tot(1, :), '-ko', 'MarkerIndices', 1:40:length(O));
	%plot([ATP_t1 ATP_t2], [-10 -10], '-k' , 'LineWidth', 2);
    gg = max(-O6 * a);
	plot([ATP_t1 ATP_t2], [-gg/2 -gg/2], '-k' , 'LineWidth', 2);
	text(ATP_t1 + (ATP_t2 - ATP_t1)/2, -gg/2 + 0.05 * gg/2, 'ATP');
	hold off;
	%text(ATP_t1 + (ATP_t2 - ATP_t1)/2, -7, 'ATP')
    ATP = ATP * 1e3;
    legend(sprintf('I_{tot} at ATP=%gmM', ATP(1)), sprintf('I_{tot} at ATP=%gmM', ATP(2)), sprintf('I_{tot} at ATP=%gmM', ATP(3)), sprintf('I_{tot} at ATP=%gmM', ATP(4)), sprintf('I_{tot} at ATP=%gmM', ATP(5)), sprintf('I_{tot} at ATP=%gmM', ATP(6)), 'FontWeight', 'bold');
	xlabel('Time (s)');
	ylabel('I_{tot} (pA)');
    xlim([0 t_max]);
    ylim([-Inf 0]);
    %k(25:27)
    return
    
    fig = figure('Name', 'State space of the P2X receptor in human microglia');
	fig.WindowStyle = 'docked';
    str = {'[C]', '[S]', '[D]', '[O]'};
    for i=1:1:length(opt_hp2x7_i_tot.x0)
       subplot(4, 2, i);
       %str = { sprintf('s%d', i) };
       plot(tint, y_hp2x7_i_tot6(i, :)); xlabel('Time (s)'); ylabel(str(i));
       %xlim([0.0 8.5]);
    end
    
    return
    
    % Repeated application of ADP/ATP
    T = 15; % Period in seconds
    nT_P2X = 0;
    
    ATP = 10e-3; ATP_t1 = 0; ATP_t2 = 3; t_max = 2000;
    
    opt_hp2x7_i_tot.ATP = ATP;
    opt_hp2x7_i_tot.ATP_t1 = ATP_t1;
    opt_hp2x7_i_tot.ATP_t2 = ATP_t2;
    opt_hp2x7_i_tot.tspan = [0.0 t_max];
    
    %return;
    my_k = k_hp2x7_i_tot;
    my_ATP = ATP;
    odeopt = odeset;%('Jacobian', @compute_jacobian_matrix_hp2x7_total_current);
    sol = ode23s(@(t, x)(reaction_network_hp2x7_total_current(t, x, opt_hp2x7_i_tot)), opt_hp2x7_i_tot.tspan, opt_hp2x7_i_tot.x0, odeopt);
    
    fig2 = figure('Name', 'Simulation for P2X7 receptors in microglia for repeated application of ADP/ATP');
	fig2.WindowStyle = 'docked';
    
    tint = linspace(0, t_max, 2000);  % all time points.
	y_hp2x7_i_tot = deval(sol, tint);
    
	O = y_hp2x7_i_tot(observable_index, :);
	hold on;
	plot(tint, O * a, '-ko', 'MarkerIndices', 1:40:length(O));
    gg = max(-O * a)
	plot([ATP_t1 ATP_t2], [-gg -gg], '-k' , 'LineWidth', 2);
	hold off;
	text(ATP_t1 + (ATP_t2 - ATP_t1)/2, -gg + 0.01 * gg/2, 'ATP');
    ATP = ATP * 1e3;
    legend(sprintf('I_{tot} at at a pulse train of ATP=%gmM', ATP));
	xlabel('Time (s)');
	ylabel('I_{tot} (pA)');
    
    fig3 = figure('Name', 'State space of the P2X receptor in human microglia');
	fig3.WindowStyle = 'docked';
    for i=1:1:length(opt_hp2x7_i_tot.x0)
       subplot(4, 2, i);
       %str = { sprintf('s%d', i) };
       plot(tint, y_hp2x7_i_tot(i, :)); xlabel('Time (s)'); ylabel(str(i));
       %xlim([0.0 8.5]);
    end
end
%--------------------%
function [K] = read_last_line()
    filename = "outcmaes-hp2x7-ATP-5.000000e-03-popsize-18-total-currentxmean.dat";
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
%--------------------%