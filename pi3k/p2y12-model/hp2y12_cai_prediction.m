%------Functions -------%
%function hp2y12_cai_prediction()
    close all; clc; clear; clearvars;
    format long g;
    %fclose('all');
	c = Constants; % Gets constants for the model
    c.should_use_matlab_solver = 1;
   
    ADP = [0 10 20 30 50 100]; % uM
    %ADP = [0 10 30 40 50 100]; % uM
    %ADP = [0 5 10 100 500 1000]; % uM
    %ADP = [50 50 50 50 50 50];
    %ADP = [50 50 50 50 50 50]; % uM
    %ADP_t1 = 2.0; ADP_t2 = 4.6; t_max = 40.0;% unit in seconds - Original
    %offset = 0; ADP_t1 = offset + 13.88; ADP_t2 = offset + 13.88 + 29.16; t_max = offset + 66.0;
    %ADP_t1 = 2.0; ADP_t2 = 41; t_max = 60.0;% unit in seconds
    %ADP_t1 = 10; ADP_t2 = 130.0; t_max = 150;% unit in seconds
    ADP_t1 = 80.0; ADP_t2 = 480; t_max = 540.0;% unit in seconds
    %ADP_t1 = 30.0; ADP_t2 = (116.0 + 30.0); t_max = 180;% unit in seconds
    %ADP_t1 = 80.0; ADP_t2 = 480; t_max = 540.0;% unit in seconds
    steps = 1000;
    if c.should_use_matlab_solver == 0
        steps = 200;
    end
   
    es_hp2y12_cai = load(sprintf('variablescmaes-hp2y12-ADP-50-popsize-%d.mat', c.PopSize), '-mat');
	%K_hp2y12_cai = es_hp2y12_cai.out.solutions.bestever.x
    %return
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
    
	if c.should_use_matlab_solver == 1
		global my_k my_ADP;
		%my_k = k_hp2y12_i_tot;
		odeopt = odeset;%('InitialStep', 1e-2, 'MaxStep', 1e-2);%('RelTol', 1e-9, 'AbsTol', 1e-9);
		%odeopt = odeset;%('Jacobian', @compute_jacobian_matrix_hp2y12_total_current);
		for i=1:1:length(ADP)
			opt_p2y_cai.ADP = ADP(i);
			my_ADP = ADP(i);
			sol_p2y_cai(i) = ode23s(@(t, x)(reaction_network_hp2y12(t, x, opt_p2y_cai)), opt_p2y_cai.tspan, opt_p2y_cai.x0, odeopt);
		end
		
		tint = linspace(0, t_max, 10000);  % all time points.
		y_p2y_cai1 = deval(sol_p2y_cai(1), tint);
		y_p2y_cai2 = deval(sol_p2y_cai(2), tint);
		y_p2y_cai3 = deval(sol_p2y_cai(3), tint);
		y_p2y_cai4 = deval(sol_p2y_cai(4), tint);

        y_p2y_cai5 = deval(sol_p2y_cai(5), tint);
		y_p2y_cai6 = deval(sol_p2y_cai(6), tint);

    else
		ts = opt_p2y_cai.tspan(1);
		tf = opt_p2y_cai.tspan(length(opt_p2y_cai.tspan));
		dt = 1e-2;
		for i=1:1:length(ADP)
			opt_p2y_cai.ADP = ADP(i);
			my_ADP = ADP(i);
			[tint, y] = rk38(@(t, x)(reaction_network_hp2y12(t, x, opt_p2y_cai)), opt_p2y_cai.x0, ts, tf, dt);
			sol_p2y_cai(:, :, i) = y;
		end
		
		y_p2y_cai1 = sol_p2y_cai(:, :, 1);
		y_p2y_cai2 = sol_p2y_cai(:, :, 2);
		y_p2y_cai3 = sol_p2y_cai(:, :, 3);
		y_p2y_cai4 = sol_p2y_cai(:, :, 4);
		y_p2y_cai5 = sol_p2y_cai(:, :, 5);
		y_p2y_cai6 = sol_p2y_cai(:, :, 6);
	end

	fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2Y12 receptor in human microglia');
	fig1.WindowStyle = 'docked'; %
    cai1 = y_p2y_cai1(opt_p2y_cai.observable_index, :);
    cai2 = y_p2y_cai2(opt_p2y_cai.observable_index, :);
    cai3 = y_p2y_cai3(opt_p2y_cai.observable_index, :);
    cai4 = y_p2y_cai4(opt_p2y_cai.observable_index, :);
    cai5 = y_p2y_cai5(opt_p2y_cai.observable_index, :);
    cai6 = y_p2y_cai6(opt_p2y_cai.observable_index, :);
	hold on;

    baseline = 0;
    plot(tint, cai1 + baseline, '-ko', 'MarkerIndices', 1:steps:length(cai1));
    plot(tint, cai2 + baseline, '-.r*', 'MarkerIndices', 1:steps:length(cai2));
    plot(tint, cai3 + baseline, '-kx', 'MarkerIndices', 1:steps:length(cai3));
    plot(tint, cai4 + baseline, 'g--', 'MarkerIndices', 1:steps:length(cai4));
    plot(tint, cai5 + baseline, 'b--*', 'MarkerIndices', 1:steps:length(cai5));
    plot(tint, cai6 + baseline, '--ko', 'MarkerIndices', 1:steps:length(cai6));
	%plot(tint, cai, '-ko', 'MarkerIndices', 1:100:length(cai));
    %plot(tint, i_tot/1, '--b', 'MarkerIndices', 1:100:length(i_tot));
    %plot(tint, cumtrapz(tint, i_tot) , '--b', 'MarkerIndices', 1:100:length(i_tot));
    %plot(tint, y_hp2y12_cai(1, :), '-ko', 'MarkerIndices', 1:40:length(O));
	plot([ADP_t1/1 ADP_t2/1], [0.2 0.2], '-k' , 'LineWidth', 2);
	%text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.22, 'ADP');
    text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.22, 'ADP');
    hold off;
    %ADP = ADP * 1e3;
    %legend(sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(1)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(2)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(3)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(4)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(5)), sprintf('\\Delta[Ca^2+_i] at ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    legend(sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(1)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(2)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(3)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(4)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(5)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    
    xlabel('Time (s)');
    ylabel('[Ca^2+_i] (\muM)');
    %ylabel('\Delta[Ca^2+_i] (\muM)');
    xlim([0 t_max]);
    ylim([0.09 Inf]);
    %ylim([44 max( [max(cai1) max(cai2) max(cai3) max(cai4) max(cai5)]) ]);%
    %return
    fig3 = figure('Name', 'Transient of the state space for P2Y12 receptor at ADP 100uM in human microglia');
	fig3.WindowStyle = 'docked';
    str = {'[IP3^{ATP}]', '[IP3]', '[C]', '[S]', '[D]', '[O]', '[Ca^2+_E_R]', '[Ca^2+_i]'};
    for i=1:1:length(opt_p2y_cai.x0)
        subplot(3, 3, i);
        %str = { sprintf('s%d', i) };
        plot(tint, y_p2y_cai6(i, :)); xlabel('Time (s)'); ylabel(str(i));
    end
    %return;
    
    if c.just_show_calcium == 1
        return;
    end
    
    k = K_hp2y12_cai;
    for i=1:1:length(ADP)
        if c.should_use_matlab_solver == 1
            y_p2y_cai = deval(sol_p2y_cai(i), tint);
        else
			y_p2y_cai = sol_p2y_cai(:, :, i);
        end
        
        D = y_p2y_cai(5, :);
        O = y_p2y_cai(6, :);
        CAs = y_p2y_cai(7, :);
        CAi = y_p2y_cai(8, :);
        alpha_leak = k(26);
        alpha_ip3r = k(27);
        V_SECRA = k(28);
        K_SECRA = k(29);
        for j=1:1:length(tint)
            J_LEAK(i, j) = leak_current(CAs(j), CAi(j), alpha_leak);
            J_IP3R(i, j) = ip3r_current(O(j), CAs(j), CAi(j), D(j), alpha_ip3r);
            J_SECRA(i, j) = secra_current(CAi(j), V_SECRA, K_SECRA);
        end
    end
    
    if c.should_use_matlab_solver == 0
        steps = 500;
    end
    fig2 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2Y12 receptor in human microglia');
	fig2.WindowStyle = 'docked';

    %subplot(2, 2, 1);
    %hold on;
    %plot(tint, cai1, '-ko', 'MarkerIndices', 1:steps:length(tint));
    %plot(tint, cai2, '-.r*', 'MarkerIndices', 1:steps:length(tint));
    %plot(tint, cai3, '-kx', 'MarkerIndices', 1:steps:length(tint));
    %plot(tint, cai4, 'g--', 'MarkerIndices', 1:steps:length(tint));
    %plot(tint, cai5, 'b--*', 'MarkerIndices', 1:steps:length(tint));
    %plot(tint, cai6, '--ko', 'MarkerIndices', 1:steps:length(tint));
    %xlabel('Time (s)');
    %ylabel('[Ca^2+_i] (\muM)');
    %xlim([0 t_max]);
    %hold off;
    %legend(sprintf('ADP=%g\\muM', ADP(1)), sprintf('ADP=%g\\muM', ADP(2)), sprintf('ADP=%g\\muM', ADP(3)), sprintf('ADP=%g\\muM', ADP(4)), sprintf('ADP=%g\\muM', ADP(5)), sprintf('ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    
    subplot(1, 3, 1);
    hold on;
    plot(tint, J_IP3R(1, :), '-ko', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_IP3R(2, :), '-.r*', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_IP3R(3, :), '-kx', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_IP3R(4, :), 'g--', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_IP3R(5, :), 'b--*', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_IP3R(6, :), '--ko', 'MarkerIndices', 1:steps:length(tint));
	plot([ADP_t1/1 ADP_t2/1], [0.4 0.4], '-k' , 'LineWidth', 2);
    text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.5, 'ADP');
    xlabel('Time (s)');
    ylabel('J_{IP3R} (\muMs^{−1})');
    xlim([0 t_max]);
    hold off;
    legend(sprintf('ADP=%g\\muM', ADP(1)), sprintf('ADP=%g\\muM', ADP(2)), sprintf('ADP=%g\\muM', ADP(3)), sprintf('ADP=%g\\muM', ADP(4)), sprintf('ADP=%g\\muM', ADP(5)), sprintf('ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    
    subplot(1, 3, 2);
    hold on;
    plot(tint, J_SECRA(1, :), '-ko', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_SECRA(2, :), '-.r*', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_SECRA(3, :), '-kx', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_SECRA(4, :), 'g--', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_SECRA(5, :), 'b--*', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_SECRA(6, :), '--ko', 'MarkerIndices', 1:steps:length(tint));
	plot([ADP_t1/1 ADP_t2/1], [0.4 0.4], '-k' , 'LineWidth', 2);
    text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.5, 'ADP');
    xlabel('Time (s)');
    ylabel('J_{SECRA} (\muMs^{−1})');
    xlim([0 t_max]);
    hold off;
    legend(sprintf('ADP=%g\\muM', ADP(1)), sprintf('ADP=%g\\muM', ADP(2)), sprintf('ADP=%g\\muM', ADP(3)), sprintf('ADP=%g\\muM', ADP(4)), sprintf('ADP=%g\\muM', ADP(5)), sprintf('ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    
    subplot(1, 3, 3);
    hold on;
    plot(tint, J_LEAK(1, :), '-ko', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_LEAK(2, :), '-.r*', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_LEAK(3, :), '-kx', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_LEAK(4, :), 'g--', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_LEAK(5, :), 'b--*', 'MarkerIndices', 1:steps:length(tint));
    plot(tint, J_LEAK(6, :), '--ko', 'MarkerIndices', 1:steps:length(tint));
	plot([ADP_t1/1 ADP_t2/1], [0.045 0.045], '-k' , 'LineWidth', 2);
    text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.0457, 'ADP');
    xlabel('Time (s)');
    ylabel('J_{Leak} (\muMs^{−1})');
    xlim([0 t_max]);
    hold off;
    legend(sprintf('ADP=%g\\muM', ADP(1)), sprintf('ADP=%g\\muM', ADP(2)), sprintf('ADP=%g\\muM', ADP(3)), sprintf('ADP=%g\\muM', ADP(4)), sprintf('ADP=%g\\muM', ADP(5)), sprintf('ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
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