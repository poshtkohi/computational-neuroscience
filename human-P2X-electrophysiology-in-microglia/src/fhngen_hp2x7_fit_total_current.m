%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [err] = fhngen_hp2x7_fit_total_current(x, opt, w_experimental)
    global rt_fig fig_states_rt err_prev counter total_elapsed_time;
    
    c = Constants; % Gets constants for the model
    opt.K(1:opt.initial_condition_index-1) = x(1:opt.initial_condition_index-1);
    opt.x0 = load_initial_conditions_from_k_hp2x7_total_current(x, opt);
    
    tStart = tic;
    % https://uk.mathworks.com/help/matlab/math/choose-an-ode-solver.html
    
    err = 0.0;
	
	tspan = [opt.tint(1) opt.tint(length(opt.tint))];
	sol = ode15s(@(t, x)(reaction_network_hp2x7_total_current(t, x, opt)), tspan, opt.x0, opt.odeopt);
	y = deval(sol, opt.tint);
	w = y(opt.observable_index, :)';
    %w = y(5, :)' + y(6, :)';
    
    %for i=1:1:length(w_experimental)
    %    h = 1 - y(3, i); 
    %    O = y(4, i);
    %    w(i) = h * O^3;
    %end
    %w
    %g = 1e-8;
    %V = 60 * 1e-3;
    %w = g * V * w;
    %w(1)
    %pause
	err = loss_function(w, w_experimental, c.loss_function);
    
    tEnd = toc(tStart);
    
    
    % https://uk.mathworks.com/matlabcentral/answers/384848-get-figures-and-use-them-to-build-a-video-avi
    %counter = counter + 1;
    %if counter == opt.PopSize
    %    %disp(sprintf("Total elapsed time for ODE solver is %g", total_elapsed_time));
    %    counter = 0;
    %    total_elapsed_time = 0.0;
    %else
    %    total_elapsed_time = total_elapsed_time + tEnd;
    %end    
    
    counter = counter + 1;
    if counter == opt.PopSize
        counter = 0;
    %if err < err_prev
        figure(rt_fig);
        clf(rt_fig, 'reset');
        %tt = linspace(0.0, opt.tint(length(opt.tint)), 500)';  % all time points
        %yy = deval(sol, tt);
        hold on;
        a = -145;
        %a = -  ;%(2.61 / 2.416);% * 145;
        %a =  -2.5 / 2.416;
        plot(opt.tint, a * w, '-.r*', 'MarkerIndices', 1:30:length(w));
        plot(opt.tint, a * w_experimental, '-bo', 'MarkerIndices', 1:30:length(w_experimental));
        plot([opt.ATP_t1 opt.ATP_t2], [a/2 a/2], '-k' , 'LineWidth', 2);
        text((opt.ATP_t1 + (opt.ATP_t2 - opt.ATP_t1) /2), a/2 - 0.07*a/2, 'ATP');
        xlabel('Time (s)');
        ylabel('I_{tot} (pA)'); 
        %ylabel('Normalised [O]');
        %xlim([0.0 8.5]);
        %ylim([0.0 0.05]);
        legend(sprintf('Fitted ATP=%dmM', opt.ATP * 1e3), sprintf('Experimental ATP=%dmM', opt.ATP * 1e3), 'FontWeight', 'bold');
        %drawnow update;
        hold off;
        %saveas(rt_fig, 'combined_max.png');
        
        if c.should_plot_state_space == 1
            str = {'C_7', 'S_7', 'D_7', 'O_7'};
            %str = {'C_1', 'S_1', 'D_1', 'O_1'};
            %str = {'[C1]', '[C2]', '[Q1]', '[Q2]', '[Q3]', '[Q4]', '[C3]', '[C4]'};
            figure(fig_states_rt);
            clf(fig_states_rt, 'reset');
            for i=1:1:length(opt.x0)
               subplot(2, 2, i);
               %str = { sprintf('s%d', i) };
               plot(opt.tint, y(i, :)); xlabel('Time (s)'); ylabel(str(i));
               %xlim([0.0 8.5]);
            end
        end
        %pause
        drawnow update;
        
        err_prev = err;
    end
    
    %pause
end
%--------------------%