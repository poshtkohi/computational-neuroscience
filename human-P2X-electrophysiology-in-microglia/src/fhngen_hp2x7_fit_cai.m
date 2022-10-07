%------Functions -------%
function [err] = fhngen_hp2x7_fit_cai(x, opt, w_experimental)
    global rt_fig fig_states_rt err_prev counter total_elapsed_time;
    
    c = Constants; % Gets constants for the model
    opt.K(1:opt.initial_condition_index-1) = x(1:opt.initial_condition_index-1);
    opt.x0 = load_initial_conditions_from_k_hp2x7_cai(x, opt);
    
    tStart = tic;
    % https://uk.mathworks.com/help/matlab/math/choose-an-ode-solver.html
    
    err = 0.0;
	
	tspan = [opt.tint(1) opt.tint(length(opt.tint))];
	sol = ode15s(@(t, x)(reaction_network_hp2x7_cai(t, x, opt)), tspan, opt.x0, opt.odeopt);
	y = deval(sol, opt.tint);
	w = y(opt.observable_index, :)';
    
    %for i=1:1:length(w)
    %    w(i) = w(i)^2;
    %end
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
        tt = linspace(0.0, opt.tint(length(opt.tint)), 500)';  % all time points
        yy = deval(sol, tt);
        hold on;
        a = 1;
        plot(tt, a * yy(opt.observable_index, :), '-.r*', 'MarkerIndices', 1:30:length(yy(opt.observable_index, :)));
        plot(opt.tint, a * w_experimental, '-bo', 'MarkerIndices', 1:30:length(w_experimental));
        %plot([opt.ATP_t1 opt.ATP_t2], [0.05 0.05], '-k' , 'LineWidth', 2);
        %text(3.0, 0.1, 'ATP');
        plot([opt.ATP_t1 opt.ATP_t2], [0.5 0.5], '-k' , 'LineWidth', 2);
        text(90/60, 0.65, 'ATP');
		xlabel('Time (m)');
        ylabel('\Delta[Ca^2+_i] (AU)');
        %ylabel('Normalised [O]');
        %xlim([0.0 8.5]);
        %ylim([0.0 0.05]);
        legend(sprintf('Fitted ATP=%dmM', opt.ATP * 1e3), sprintf('Experimental ATP=%dmM', opt.ATP * 1e3));
        %drawnow update;
        hold off;
        %saveas(rt_fig, 'combined_max.png');
        
        str = {'\Delta[C]', '\Delta[S]', '\Delta[D]', '\Delta[O]', '\Delta[CAi]', '\Delta[NCX]', '\Delta[NCX3-CA]', '\Delta[PMCA]', '\Delta[PMCA-CA]'};
        figure(fig_states_rt);
        clf(fig_states_rt, 'reset');
        for i=1:1:length(opt.x0)
           subplot(5, 2, i);
           %str = { sprintf('s%d', i) };
           plot(tt, yy(i, :)); xlabel('Time (m)'); ylabel(str(i));
           xlim([0.0 160/60]);
        end
        
        drawnow update;
        
        err_prev = err;
    end
end
%--------------------%