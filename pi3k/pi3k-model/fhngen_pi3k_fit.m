%------Functions -------%
function [err] = fhngen_pi3k_fit(x, opt)%, w_experimental)
    global rt_fig fig_states_rt err_prev counter total_elapsed_time;
    global f1 w_experimental first;
    
    c = Constants; % Gets constants for the model
    opt.K(1:opt.initial_condition_index-1) = x(1:opt.initial_condition_index-1);
    opt.x0 = load_initial_conditions_from_k_pi3k(x, opt);
    k = x;
    
    % Prepares delayed P2X-mediated Ca2+
    td = k(25) + 1e-6;% 2.7703;%
    %%td
    %td = 2.7703;
    td = td * 60;
    %td = 2.7703 * 60;
    tint_new = zeros(1, 1002);
    cai_new = opt.baseline * 1e-9 * ones(1, 1002);
    cai_new(2:1001) = opt.cai_p2x;
    tint_new(2:1001) = opt.tint_p2x + td;
    tint_new(1002) = tint_new(1001) + 1 * 60;
    mode = 'linearinterp'; % linearinterp, smoothingspline, cubicinterp
    f = fit(tint_new', cai_new', mode);
    opt.f = f;
    
    %disp('hello world');
    %pause
    tStart = tic;
    % https://uk.mathworks.com/help/matlab/math/choose-an-ode-solver.html
    
    err = 0.0;
	
	tspan = [opt.tint(1) opt.tint(length(opt.tint))];
    %{
    if c.solver == 15

        sol = ode15s(@(t, x)(reaction_network_pi3k(t, x, opt)), tspan, opt.x0, opt.odeopt);
    else
        sol = ode23s(@(t, x)(reaction_network_pi3k(t, x, opt)), tspan, opt.x0, opt.odeopt);
    end
	y = deval(sol, opt.tint);
	w = y(opt.observable_index, :)';
    %w = k(14)^0 * y(3, :)' + y(4, :)';
    %}
    %----------------------------------
    if c.should_use_matlab_solver == 1
        warning('error', 'MATLAB:ode15s:IntegrationTolNotMet');
        %x(40)
        %tic;
        try
           %disp('ode15s');
           if c.explicit_solver == 1
                sol = ode45(@(t, x)(reaction_network_pi3k(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
           else
               sol = ode15s(@(t, x)(reaction_network_pi3k(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
               %disp('ode15s');
           end
        catch
            try
                disp('ode23s');
                sol = ode23s(@(t, x)(reaction_network_pi3k(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
            catch
                %disp(lasterr);
                c.should_use_matlab_solver = 0;
            end
        end
        %if c.matlab_solver == 45
        %    sol = ode45(@(t, x)(reaction_network_hp2y12(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
        %elseif c.matlab_solver == 15
        %    sol = ode15s(@(t, x)(reaction_network_hp2y12(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
        %elseif c.matlab_solver == 23
        %    sol = ode23s(@(t, x)(reaction_network_hp2y12(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
        %end
        if c.should_use_matlab_solver == 1
            y = deval(sol, opt.tint);
            w = y(opt.observable_index, :)';
            %mm = max(y(2, :));
        end
        %toc;
	end
    
    if c.should_use_matlab_solver == 0
		ts = opt.tspan(1);
		tf = opt.tspan(length(opt.tspan));
        %%[ts tf]
		dt = 1e-2;%heun rk4 euler rk38
        %tic;
		[tint, y] = rk38(@(t, x)(reaction_network_pi3k(t, x, opt)), opt.x0, ts, tf, dt);
        %toc;
        w = y(opt.observable_index, :)';
        %mm = max(y(2, :));
        if first == true
            w_experimental = f1(tint);
            first = false;
        end
    end
    %----------------------------------
	err = loss_function(w, w_experimental, c.loss_function);
    
    tEnd = toc(tStart);
    
    if c.should_plot == 0
        return;
    end
    counter = counter + 1;
    if counter == opt.PopSize
        counter = 0;
        figure(rt_fig);
        clf(rt_fig, 'reset');
        %tt = linspace(0.0, opt.tint(length(opt.tint)), 500)';  % all time points
        %yy = deval(sol, tt);
        hold on;
        %plot(tt, yy(opt.observable_index, :), '-.r*', 'MarkerIndices', 1:30:length(yy(opt.observable_index, :)));
        %%%plot(opt.tint, w, '-.r*', 'MarkerIndices', 1:30:length(w));
        %%%plot(opt.tint, w_experimental, '-bo', 'MarkerIndices', 1:30:length(w_experimental));
		if c.should_use_matlab_solver == 1
			tt = linspace(opt.tspan(1), opt.tspan(length(opt.tspan)), 1000)';  % all time points
			yy = deval(sol, tt);
			plot(tt, yy(opt.observable_index, :), '-.r*', 'MarkerIndices', 1:30:length(yy));
            plot(opt.tint, w_experimental, '-bo', 'MarkerIndices', 1:30:length(w_experimental));
        end
		if c.should_use_matlab_solver == 0
			%plot(opt.tint, w, '-x', 'MarkerIndices', 1:5:length(w));
            plot(tint, w, '-.r*', 'MarkerIndices', 1:500:length(w));
            plot(tint, w_experimental, '-bo', 'MarkerIndices', 1:500:length(w_experimental));
            yy = y;
            tt = tint;
        end
		
        mm = max(w_experimental);
        plot([opt.ADP_t1 opt.ADP_t2], [mm/2 mm/2], '-k' , 'LineWidth', 2);
        text((opt.ADP_t1 + (opt.ADP_t2 - opt.ADP_t1) /2), mm/2 - 0.07*mm/2, 'ADP');
        xlabel('Time (m)');
        ylabel('pAkt (\muM)');
        %ylabel('pAkt (\muM)');
        %ylabel('Normalised [O]');
        %xlim([0.0 8.5]);
        %ylim([0.0 0.05]);
        legend(sprintf('Fitted ADP=%g\\muM', opt.ADP), sprintf('Experimental ADP=%g\\muM', opt.ADP), 'FontWeight', 'bold');
        %drawnow update;
        hold off;
        
        %%{
        if c.should_plot_state_space == 1
            str = {'CaMK_{trap}', 'PIP_2 (\muM)', 'PI3K (\muM)', 'PIP_2-PI3K (\muM)', 'PIP_3 (\muM)', 'PTEN (\muM)', 'PIP_3-PTEN (\muM)', 'AKT (\muM)', 'PIP_3-AKT (\muM)', 'pAKT (\muM)', 'PP2A (\muM)', 'pAKT-PP2A (\muM)', 'PDK1 (\muM)', 'PDK1-PIP_3 (\muM)', 'PIP_3-PDK1-Akt (\muM)', 'CaMK_{active}', '[Ca^2+_i]_{tot} (\muM)'};
            figure(fig_states_rt);
            clf(fig_states_rt, 'reset');
            for i=2:1:length(opt.x0)
               subplot(4, 4, i-1);
               plot(tt, yy(i, :)); xlabel('t (m)'); ylabel(str(i));
               %xlim([0.0 8.5]);
            end

            %%{
            y_p2y_cai = deval(opt.sol_p2y_cai, tt * 60);  % minute to seconds
            cai_p2y = y_p2y_cai(8, :) - opt.baseline * 0.001;
            cai_p2x = opt.f(tt * 60)' * 1e6 - opt.baseline * 0.001; % M to uM
            cai_tot = calculate_total_calcium(cai_p2y, cai_p2x, k);
            i = i + 1;
            CaMK_active  = calculate_CaMK_active(cai_tot, yy(1, :), opt.K, opt);
            subplot(4, 4, i - 1);
            plot(tt, CaMK_active); xlabel('t (m)'); ylabel(str(i));
            i = i + 1;
            subplot(4, 4, i - 1);
            plot(tt, cai_tot + opt.baseline * 0.001); xlabel('t (m)'); ylabel(str(i));  
            %%}

        end
        %%}
        drawnow update;
        
        err_prev = err;
        %pause %%
    end
end
%--------------------%