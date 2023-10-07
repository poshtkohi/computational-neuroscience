%------Functions -------%
%function hp2y12_cai_prediction()
    %--------------- clean up simulation enviroment ---------------%
    close all; clc; clear; clearvars;
    format long g;
    %fclose('all');
	%--------------- global settings ---------------%
	%ADP = 1000; % unit in micromollars
	%ADP = [0 10 30 50 80 50]; % uM
    ADP = [0 20 30 50 60 100]; % uM
    %ADP = [50 50 50 50 50 50]; % uM
	offset = 20 / 60;
    ADP_t1 = 0 + offset; ADP_t2 = 30 + offset; t_max = (2000 + 100)/60;%t_max = 30 * 60 + offset;% unit in seconds
    %ADP_t1 = 2*60; ADP_t2 = 30 * 60 ; t_max = 35 * 60;% unit in seconds
    %ADP_t1 = 100/60; ADP_t2 = 200/60 ; t_max = 1000/60;% unit in seconds
	c = Constants; % Gets constants for the model
    steps = c.steps;
    baseline = 100; %45
    opt_p2x_cai.baseline = baseline;
    %--------------- simulate the P2Y model ---------------%
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
	opt_p2y_cai.ADP_t1 = ADP_t1 * 60; % Convert mintues to seconds
	opt_p2y_cai.ADP_t2 = ADP_t2 * 60; % Convert mintues to seconds
	opt_p2y_cai.tspan = [0.0 t_max * 60]; % Convert mintues to seconds
    opt_p2y_cai.observable_index = 8;
	odeopt = odeset;%('InitialStep', 1e-2, 'MaxStep', 1e-2);%('RelTol', 1e-9, 'AbsTol', 1e-9);
	
	for i=1:1:length(ADP)
		opt_p2y_cai.ADP = ADP(i);
		sol_p2y_cai(i) = ode23s(@(t, x)(reaction_network_hp2y12(t, x, opt_p2y_cai)), opt_p2y_cai.tspan, opt_p2y_cai.x0, odeopt);
	end
		
	tint = linspace(opt_p2y_cai.tspan(1), opt_p2y_cai.tspan(2), 100000);  % all time points.
	y_p2y_cai1 = deval(sol_p2y_cai(1), tint);
	y_p2y_cai2 = deval(sol_p2y_cai(2), tint);
	y_p2y_cai3 = deval(sol_p2y_cai(3), tint);
	y_p2y_cai4 = deval(sol_p2y_cai(4), tint);
	y_p2y_cai5 = deval(sol_p2y_cai(5), tint);
	y_p2y_cai6 = deval(sol_p2y_cai(6), tint);
	
	caiy(1, :) = y_p2y_cai1(opt_p2y_cai.observable_index, :);
    caiy(2, :) = y_p2y_cai2(opt_p2y_cai.observable_index, :);
    caiy(3, :) = y_p2y_cai3(opt_p2y_cai.observable_index, :);
    caiy(4, :) = y_p2y_cai4(opt_p2y_cai.observable_index, :);
    caiy(5, :) = y_p2y_cai5(opt_p2y_cai.observable_index, :);
    caiy(6, :) = y_p2y_cai6(opt_p2y_cai.observable_index, :);

	fig1 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2Y12 receptor in human microglia');
	fig1.WindowStyle = 'docked';
	hold on;
	plot(tint/60, caiy(1, :), '-ko', 'MarkerIndices', 1:steps:length(caiy(1, :)));
    plot(tint/60, caiy(2, :), '-.r*', 'MarkerIndices', 1:steps:length(caiy(2, :)));
    plot(tint/60, caiy(3, :), '-kx', 'MarkerIndices', 1:steps:length(caiy(3, :)));
    plot(tint/60, caiy(4, :), 'g--', 'MarkerIndices', 1:steps:length(caiy(4, :)));
    plot(tint/60, caiy(5, :), 'b--*', 'MarkerIndices', 1:steps:length(caiy(5, :)));
    plot(tint/60, caiy(6, :), '--ko', 'MarkerIndices', 1:steps:length(caiy(6, :)));
	plot([ADP_t1/1 ADP_t2/1], [0.2 0.2], '-k' , 'LineWidth', 2);
	text((ADP_t1 + (ADP_t2 - ADP_t1)/2), 0.22, 'ADP');
    hold off;
    legend(sprintf('[Ca^2+_i] at ATP=%g\\muM', ADP(1)), sprintf('[Ca^2+_i] at ATP=%g\\muM', ADP(2)), sprintf('[Ca^2+_i] at ATP=%g\\muM', ADP(3)), sprintf('[Ca^2+_i] at ATP=%g\\muM', ADP(4)), sprintf('[Ca^2+_i] at ATP=%g\\muM', ADP(5)), sprintf('[Ca^2+_i] at ATP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    xlabel('t (m)');
    ylabel('[Ca^2+_i] (\muM)');
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
    %opt_p2x_cai.ATP = ADP * 1e-6; % uM to M
    opt_p2x_cai.ATP_t1 = ADP_t1 * 60; % Convert mintues to seconds
    opt_p2x_cai.ATP_t2 = ADP_t2 * 60; % Convert mintues to seconds
    opt_p2x_cai.tspan = [0.0 t_max * 60]; % Convert mintues to seconds
    opt_p2x_cai.observable_index = 10;
    opt_p2x_cai.test = false;
    opt_p2x_cai.CAiB = opt_p2x_cai.x0(10, 1);
    
    odeopt = odeset;%('InitialStep', 1e-2, 'MaxStep', 1e-2);%('RelTol', 1e-9, 'AbsTol', 1e-9);
    
    for i=1:1:length(ADP)
		opt_p2x_cai.ATP = ADP(i) * 1e-6;
		sol_p2x_cai(i) = ode23tb(@(t, x)(reaction_network_p2x_cai_prediction(t, x, opt_p2x_cai)), opt_p2x_cai.tspan, opt_p2x_cai.x0, odeopt);
    end
    
    tintx = linspace(opt_p2x_cai.tspan(1), opt_p2x_cai.tspan(2), 100000);  % all time points.
    
    y_p2y_cai1 = deval(sol_p2x_cai(1), tintx);
    y_p2y_cai2 = deval(sol_p2x_cai(2), tintx);
    y_p2y_cai3 = deval(sol_p2x_cai(3), tintx);
    y_p2y_cai4 = deval(sol_p2x_cai(4), tintx);
    y_p2y_cai5 = deval(sol_p2x_cai(5), tintx);
    y_p2y_cai6 = deval(sol_p2x_cai(6), tintx);
    
    caix(1, :) = y_p2y_cai1(opt_p2x_cai.observable_index, :);
    caix(2, :) = y_p2y_cai2(opt_p2x_cai.observable_index, :);
    caix(3, :) = y_p2y_cai3(opt_p2x_cai.observable_index, :);
    caix(4, :) = y_p2y_cai4(opt_p2x_cai.observable_index, :);
    caix(5, :) = y_p2y_cai5(opt_p2x_cai.observable_index, :);
    caix(6, :) = y_p2y_cai6(opt_p2x_cai.observable_index, :);
    %opt_p2x_cai.cai = cai;
    %opt_p2x_cai.tint_cai = tintx;
    
    fig2 = figure('Name', 'Simulation of intracellular/cytosolic [CAi] concentration for P2X receptors in human microglia');
	fig2.WindowStyle = 'docked';
	hold on;
	plot(tintx/60, caix(1, :) * 1e6, '-ko', 'MarkerIndices', 1:steps:length(caix(1, :)));
    plot(tintx/60, caix(2, :) * 1e6, '-.r*', 'MarkerIndices', 1:steps:length(caix(2, :)));
    plot(tintx/60, caix(3, :) * 1e6, '-kx', 'MarkerIndices', 1:steps:length(caix(3, :)));
    plot(tintx/60, caix(4, :) * 1e6, 'g--', 'MarkerIndices', 1:steps:length(caix(4, :)));
    plot(tintx/60, caix(5, :) * 1e6, 'b--*', 'MarkerIndices', 1:steps:length(caix(5, :)));
    plot(tintx/60, caix(6, :) * 1e6, '--ko', 'MarkerIndices', 1:steps:length(caix(6, :)));
    %mm = max([caix1' * 1e6; caix2' * 1e6; caix3' * 1e6; caix4' * 1e6; caix5' * 1e6; caix6' * 1e6]);
    %%plot([ADP_t1 ADP_t2], [0.065 0.065], '-k' , 'LineWidth', 2);
    %%text((ADP_t1 + (ADP_t2 - ADP_t1) /2), 0.064, 'ADP');
    hold off;
    legend(sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(1)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(2)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(3)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(4)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(5)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    xlabel('t (m)');
    ylabel('[Ca^2+_i] (\muM)');
    %------ Simulate the PI3K model activated by P2X receptor -------%
    es_pi3k = load(sprintf('variablescmaes-pi3k-ADP-50-popsize-%d.mat', c.PopSize), '-mat');
	K_pi3k = es_pi3k.out.solutions.bestever.x;
    %K_pi3k = [1.286440e+02 3.218806e-01 1.636501e+00 8.463720e+00 9.117736e+00 2.392495e-01 2.267273e+00 1.484790e-02 3.637311e-01 1.221563e+01 6.249985e-01 2.999939e+01 1.828263e+01 4.238167e+01 4.826946e+01 7.079491e+00 2.857487e+00 1.468633e+05 4.815618e+03 4.147397e+00 1.681926e-02 3.637156e-01 1.627701e+01 1.203207e+01 1.419787e+02 1.123689e+00 -1.488650e-03 7.390988e+01 1.074026e+00 8.031701e+00 5.951535e+00 8.397641e+00 1.000689e+01 6.241643e+01 6.916468e+02 1.564680e+01 2.680433e+02 6.940629e+02 1.269101e+02 7.868626e+02 7.624969e+00 2.887579e+02 1.852644e+02];
    K_pi3k = [2.223930e+03 1.017092e+00 2.622925e+01 2.255357e+01 5.876924e+01 2.482221e+00 5.350945e+00 3.216568e+01 8.663748e-02 1.344500e+00 4.146379e+00 7.480792e-01 4.524940e-03 1.046063e+00 1.408488e+01 1.280787e+01 8.048355e+00 2.985973e+00 1.074286e+02 6.621960e-01 7.337991e+00 1.070576e+02 1.486438e+02 6.937988e+01 3.227989e+00 2.781007e+01 1.709805e+01 1.024117e+03 1.078698e+02 2.589290e+02 6.933637e+02 1.572939e+02 9.981844e+02];
    K_init = load('data/initial-guess-k-pi3k.dat');
    %K_pi3k = K_init;%%%%
    opt_pi3k.nn = length(K_init);
    
    k_pi3k = K_pi3k(1:opt_pi3k.nn);
	opt_pi3k.x0 = zeros(15, 1);
    opt_pi3k.initial_condition_index = length(k_pi3k) + 1;
	opt_pi3k.x0 = load_initial_conditions_from_k_pi3k(K_pi3k, opt_pi3k);
    %%opt_pi3k.sol_p2y_cai = sol_p2y_cai;
    opt_pi3k.K = k_pi3k;
	opt_pi3k.ADP = ADP;
	opt_pi3k.ADP_t1 = ADP_t1; % unit in seconds
	opt_pi3k.ADP_t2 = ADP_t2; % unit in seconds
	opt_pi3k.tspan = [0 t_max]; % unit in seconds
    opt_pi3k.observable_index = 10;
    %opt_pi3k.mode = c.pi3k_mode;
    opt_pi3k.CAiB_p2y = opt_p2y_cai.x0(opt_p2y_cai.observable_index); % Steady-state baseline of CAi: units in uM
    opt_pi3k.CAiB_p2x = opt_p2x_cai.x0(opt_p2x_cai.observable_index); % Steady-state baseline of CAi: units in uM
    t_d = k_pi3k(25) * 60; % minutes to seconds
    k_pi3k(25)
    
    odeopt = odeset;%('InitialStep', 1e-2, 'MaxStep', 1e-2);%('RelTol', 1e-9, 'AbsTol', 1e-9);
    for i=1:1:length(ADP)
		%opt_pi3k_p2x.sol_p2y_cai = sol_p2y_cai(i);
        %opt_p2x_cai.cai = caix(i, :);
        cai = caix(i, :);
        tint = tintx;
        tint_new = zeros(1, 100002);
        cai_new = baseline * 1e-9 * ones(1, 100002);
        cai_new(2:100001) = cai;
        tint_new(2:100001) = tint + t_d;
        tint_new(100002) = tint_new(100001) + 1 * 60;
        mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
        opt_pi3k.f = fit(tint_new', cai_new', mode);
        if i == 1
            f1 = opt_pi3k.f;
        elseif i == 2
            f2 = opt_pi3k.f;
        elseif i == 3
            f3 = opt_pi3k.f;
        elseif i == 4
            f4 = opt_pi3k.f;
        elseif i == 5
            f5 = opt_pi3k.f;
        elseif i == 6
            f6 = opt_pi3k.f;
        end
        opt_pi3k.sol_p2y_cai = sol_p2y_cai(i);
        sol_pi3k(i) = ode15s(@(t, x)(reaction_network_pi3k(t, x, opt_pi3k)), opt_pi3k.tspan, opt_pi3k.x0, odeopt);
	end
    
    tint = linspace(opt_pi3k.tspan(1), opt_pi3k.tspan(2), 100000);  % all time points.
    
    pi3k1 = deval(sol_pi3k(1), tint);
	pi3k2 = deval(sol_pi3k(2), tint);
	pi3k3 = deval(sol_pi3k(3), tint);
	pi3k4 = deval(sol_pi3k(4), tint);
	pi3k5 = deval(sol_pi3k(5), tint);
	pi3k6 = deval(sol_pi3k(6), tint);
	
	pAkt1 = pi3k1(opt_pi3k.observable_index, :);
    pAkt2 = pi3k2(opt_pi3k.observable_index, :);
    pAkt3 = pi3k3(opt_pi3k.observable_index, :);
    pAkt4 = pi3k4(opt_pi3k.observable_index, :);
    pAkt5 = pi3k5(opt_pi3k.observable_index, :);
    pAkt6 = pi3k6(opt_pi3k.observable_index, :);
    
    %CaMK_active = calculate_CaMK_active(cai, CaMK_trap, opt_pi3k.K, opt_pi3k);
    fig3 = figure('Name', 'Simulation of pAkt concentration in microglia');
	fig3.WindowStyle = 'docked';
    hold on;
    plot(tint, pAkt1, '-ko', 'MarkerIndices', 1:steps:length(pAkt1));
    plot(tint, pAkt2, '-.r*', 'MarkerIndices', 1:steps:length(pAkt2));
    plot(tint, pAkt3, '-kx', 'MarkerIndices', 1:steps:length(pAkt3));
    plot(tint, pAkt4, 'g--', 'MarkerIndices', 1:steps:length(pAkt4));
    plot(tint, pAkt5, 'b--*', 'MarkerIndices', 1:steps:length(pAkt5));
    plot(tint, pAkt6, '--ko', 'MarkerIndices', 1:steps:length(pAkt6));
    mm = max([pAkt1'; pAkt2'; pAkt3'; pAkt4'; pAkt5'; pAkt6']);
    %plot([ADP_t1 ADP_t2], [mm/2 mm/2], '-k' , 'LineWidth', 2);   
    %text((ADP_t1 + (ADP_t2 - ADP_t1) /2), mm/2 - 0.07*mm/2, 'ADP');
    hold off;
    legend(sprintf('[pAkt] at ADP=%g\\muM', ADP(1)), sprintf('[pAkt] at ADP=%g\\muM', ADP(2)), sprintf('[pAkt] at ADP=%g\\muM', ADP(3)), sprintf('[pAkt] at ADP=%g\\muM', ADP(4)), sprintf('[pAkt] at ADP=%g\\muM', ADP(5)), sprintf('[pAkt] at ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    %legend(sprintf('[pAkt] at ADP=%g\\muM', ADP(1)), 'FontWeight', 'bold');
    xlabel('t (m)');
    ylabel('[pAkt] (\muM)');
    ylim([0 Inf]);
    %------ Total Ca2+ -------%
    fig3 = figure('Name', 'Simulation of total intracellular/cytosolic [CAi] concentration for P2 receptors in microglia');
	fig3.WindowStyle = 'docked';
	hold on;
	plot(tint, baseline * 1e-3 + calculate_total_calcium(caiy(1, :) - baseline * 1e-3, f1(tint * 60)' * 1e6 - baseline * 1e-3, opt_pi3k.K), '-ko', 'MarkerIndices', 1:steps:length(caix(1, :)));
    plot(tint, baseline * 1e-3 + calculate_total_calcium(caiy(2, :) - baseline * 1e-3, f2(tint * 60)' * 1e6 - baseline * 1e-3, opt_pi3k.K), '-.r*', 'MarkerIndices', 1:steps:length(caix(2, :)));
    plot(tint, baseline * 1e-3 + calculate_total_calcium(caiy(3, :) - baseline * 1e-3, f3(tint * 60)' * 1e6 - baseline * 1e-3, opt_pi3k.K), '-kx', 'MarkerIndices', 1:steps:length(caix(3, :)));
    plot(tint, baseline * 1e-3 + calculate_total_calcium(caiy(4, :) - baseline * 1e-3, f4(tint * 60)' * 1e6 - baseline * 1e-3, opt_pi3k.K), 'g--', 'MarkerIndices', 1:steps:length(caix(4, :)));
    plot(tint, baseline * 1e-3 + calculate_total_calcium(caiy(5, :) - baseline * 1e-3, f5(tint * 60)' * 1e6 - baseline * 1e-3, opt_pi3k.K), 'b--*', 'MarkerIndices', 1:steps:length(caix(5, :)));
    plot(tint, baseline * 1e-3 + calculate_total_calcium(caiy(6, :) - baseline * 1e-3, f6(tint * 60)' * 1e6 - baseline * 1e-3, opt_pi3k.K), '--ko', 'MarkerIndices', 1:steps:length(caix(6, :)));
    %mm = max([caix1' * 1e6; caix2' * 1e6; caix3' * 1e6; caix4' * 1e6; caix5' * 1e6; caix6' * 1e6]);
    plot([ADP_t1 ADP_t2], [0.4 0.4], '-k' , 'LineWidth', 2);
    text((ADP_t1 + (ADP_t2 - ADP_t1) /2), 0.43, 'ADP');
    hold off;
    legend(sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(1)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(2)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(3)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(4)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(5)), sprintf('[Ca^2+_i] at ADP=%g\\muM', ADP(6)), 'FontWeight', 'bold');
    xlabel('t (m)');
    ylabel('[Ca^2+_i]_{tot} (\muM)');
    ylim([0.1 Inf]);
    %------ State variables -------%
    fig4 = figure('Name', 'State variables of the PI3K/Akt model in microglia');
	fig4.WindowStyle = 'docked';
    str = {'CaMK_{trap}', 'PIP_2 (\muM)', 'PI3K (\muM)', 'PIP_2-PI3K (\muM)', 'PIP_3 (\muM)', 'PTEN (\muM)', 'PIP_3-PTEN (\muM)', 'AKT (\muM)', 'PIP_3-AKT (\muM)', 'pAKT (\muM)', 'PP2A (\muM)', 'pAKT-PP2A (\muM)', 'PDK1 (\muM)', 'PDK1-PIP_3 (\muM)', 'PIP_3-PDK1-Akt (\muM)', 'CaMK_{active}', '[Ca^2+_i]_{tot} (\muM)'};
    for i=2:1:length(opt_pi3k.x0)
        subplot(5, 4, i-1);
        plot(tint, pi3k6(i, :)); xlabel('t (m)'); ylabel(str(i));
        %xlim([0.0 8.5]);
    end

    %%{
    y_p2y_cai = deval(sol_p2y_cai(6), tint * 60);  % minute to seconds
    cai_p2y = y_p2y_cai(8, :) - baseline * 1e-3;
    cai_p2x = opt_pi3k.f(tint * 60)' * 1e6 - baseline * 1e-3; % M to uM
    cai_tot = calculate_total_calcium(cai_p2y, cai_p2x, opt_pi3k.K);
    i = i + 1;
    CaMK_active  = calculate_CaMK_active(cai_tot, pi3k6(1, :), opt_pi3k.K, opt_pi3k);
    subplot(5, 4, i - 1);
    plot(tint, CaMK_active); xlabel('t (m)'); ylabel(str(i));
    i = i + 1;
    subplot(5, 4, i - 1);
    plot(tint, cai_tot + baseline * 1e-3); xlabel('t (m)'); ylabel(str(i));

%end
%--------------------%
function [g_CA_Leak, g_NA_Leak] = compute_leak_conductances(opt_p2x_cai)
    CAx = 2e-3;   % M
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