function [err] = fhngen_hp2y12_fit(x, opt)%, w_experimental)
    global rt_fig fig_states_rt err_prev counter total_elapsed_time;
    c = Constants; % Gets constants for the model
    global f1 w_experimental first;
    
    opt.K(1:opt.nn) = x(1:opt.nn);
    opt.x0 = load_initial_conditions_from_k_hp2y12(x, opt);
    
    %%tStart = tic;
    % https://uk.mathworks.com/help/matlab/math/choose-an-ode-solver.html
    % Set a couple of warnings to temporarily issue errors (exceptions)
    % ode23s, ode15s, ode23tb, ode23t
    %tic;ode23tb
    %%{
	if c.should_use_matlab_solver == 1
        warning('error', 'MATLAB:ode15s:IntegrationTolNotMet');
        %x(40)
        %tic;
        try
           %disp('ode15s');
           if c.explicit_solver == 1
                sol = ode45(@(t, x)(reaction_network_hp2y12(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
           else
               sol = ode15s(@(t, x)(reaction_network_hp2y12(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
           end
        catch
            try
                disp('ode23s');
                sol = ode23s(@(t, x)(reaction_network_hp2y12(t, x, opt)), opt.tspan, opt.x0, opt.odeopt);
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
    %%}
	
    %%{
	if c.should_use_matlab_solver == 0
		ts = opt.tspan(1);
		tf = opt.tspan(length(opt.tspan));
		dt = 1e-2;%heun rk4 euler rk38
        %tic;
		[tint, y] = rk38(@(t, x)(reaction_network_hp2y12(t, x, opt)), opt.x0, ts, tf, dt);
        %toc;
        w = y(opt.observable_index, :)';
        %mm = max(y(2, :));
        if first == true
            w_experimental = f1(tint);
            first = false;
        end
	end
    %%}
    
    err = loss_function(w, w_experimental, c.loss_function);
    
    %return;
    % https://uk.mathworks.com/matlabcentral/answers/384848-get-figures-and-use-them-to-build-a-video-avi
    %return;
    %if err < err_prev

    if c.should_plot == 0
        return
    end
    counter = counter + 1;
    %if counter == opt.PopSize
        counter = 0;
        figure(rt_fig);
        clf(rt_fig, 'reset');
		hold on
        baseline = 0;
		if c.should_use_matlab_solver == 1
			tt = linspace(opt.tspan(1), opt.tspan(length(opt.tspan)), 1000)';  % all time points
			yy = deval(sol, tt);
			plot(tt, yy(opt.observable_index, :) + baseline, '-x', 'MarkerIndices', 1:40:length(yy));
            plot(opt.tint, w_experimental + baseline, '-o', 'MarkerIndices', 1:100:length(w_experimental));
        end
		if c.should_use_matlab_solver == 0
			%plot(opt.tint, w, '-x', 'MarkerIndices', 1:5:length(w));
            plot(tint, w, '-x', 'MarkerIndices', 1:500:length(w));
            plot(tint, w_experimental, '-o', 'MarkerIndices', 1:500:length(w_experimental));
            yy = y;
            tt = tint;
        end
        plot([opt.ADP_t1/1 opt.ADP_t2/1], [0.2 0.2], '-k' , 'LineWidth', 2);
        %text((opt.ADP_t1 + (opt.ADP_t2 - opt.ADP_t1)/2), 0.22, 'ATP');
        text((opt.ADP_t1 + (opt.ADP_t2 - opt.ADP_t1)/2), 0.22, 'ADP');
        xlabel('Time (s)');
        %ylabel('\Delta[Ca^2+_i] (\muM)');
        ylabel('[Ca^2+_i] (\muM)');
        legend(sprintf('Fitted ADP=%d\\muM', opt.ADP), sprintf('Experimental ADP=%d\\muM', opt.ADP));
        %legend(sprintf('Fitted ATP=%d\\muM', opt.ADP), sprintf('Experimental ATP=%d\\muM', opt.ADP));
        xlim([0.0 opt.tspan(2)]);
        ylim([0.09 Inf]);
        hold off;
        
        if c.should_plot_state_space == 1
            str = {'[IP3^{ATP}]', '[IP3]', 'C', 'S', 'h', 'O', '[Ca^2+_E_R]', '[Ca^2+_i]'};
            figure(fig_states_rt);
            clf(fig_states_rt, 'reset');
            j = 0;
            for i=1:1:length(opt.x0)-1
                if i == 7
                    j = 1;
                    continue;
                end
                subplot(2, 3, i - j);
                %str = { sprintf('s%d', i) };
                if i == 5
                    plot(tt, 1 - yy(i, :)); xlabel('Time (s)'); ylabel(str(i));
                else
                    plot(tt, yy(i, :)); xlabel('Time (s)'); ylabel(str(i));
                end
                xlim([0.0 opt.tspan(2)]);
            end    
        end
        
        drawnow update;
        
        err_prev = err;
    %end
    
    %pause
    %counter = counter + 1;
    %if counter == opt.PopSize
    %    disp(sprintf("Total elapsed time for ODE solver is %g", total_elapsed_time));
    %    counter = 0;
    %    total_elapsed_time = 0.0;
    %else
    %    total_elapsed_time = total_elapsed_time + tEnd;
    %end
        
end


% Heun's method
function [t, y] = heun(odefun, y0, ts, tf, dt)
t = ts:dt:tf;
y = zeros(length(y0),length(t));
y(:,1) = y0;

for k = 2 : length(t)
    yn = y(:,k-1);
    tn = t(k-1);
    k1 = dt * odefun(tn, yn);
    k2 = dt * odefun(tn + dt, yn + k1);
    
    y(:,k) = yn + 1 / 2 * k1 + 1 / 2 * k2;
end
end


% Euler method
function [t, y] = euler(odefun, y0, ts, tf, dt)
t = ts:dt:tf;
y = zeros(length(y0),length(t));
y(:,1) = y0;
for k = 2 : length(t)
    ydot = odefun(t(k), y(:,k-1));
    y(:,k) = y(:,k-1)+ydot.*dt;
end
end

% Fourth oder Runge-Kutta
function [t, y] = rk4(odefun, y0, tstart, tfinal, dt)
t = tstart:dt:tfinal;
y = zeros(length(y0),length(t));
y(:,1) = y0;

for k = 2 : length(t)
    yn = y(:,k-1);
    tn = t(k-1);
    k1 = dt * odefun(tn, yn);
    k2 = dt * odefun(tn + dt / 2 , yn + k1 / 2);
    k3 = dt * odefun(tn + dt / 2, yn + k2 / 2);
    k4 = dt * odefun(tn + dt, yn + k3);
    
    y(:,k) = yn + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6;
end
end