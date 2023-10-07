% Global Parameters
function [val] = Constants
% -- Simulation Parameters --%
val.PopSize = 18;
val.loss_function = 2;
val.loss_function_division = true;
val.resume = 'yes'; % yes or no.
val.should_use_matlab_solver = 1;
val.explicit_solver = 0;
val.should_plot = 1;
val.should_plot_state_space = 1;
% PI3K settings
val.pi3k_mode = 0; % both humps (0), first hump (1) and second hump (2)
val.steps = 10000;
val.sigma = 0.5;
val.periodic = 0;
val.delay = 20; % minute
end