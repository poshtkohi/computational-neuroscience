% Global Parameters
function [val] = Constants
% -- Simulation Parameters --%
val.presimulation_mode = 0;
val.curefitting_mode = 1;
val.simulation_mode = 2;

val.PopSize = 18;
val.loss_function = 1;
val.loss_function_division = true;
val.resume = 'yes'; % yes or no.
val.should_use_matlab_solver = 1;
val.explicit_solver = 0;
val.should_plot = 1;
val.should_plot_state_space = 1;
val.just_show_calcium = 0;
val.should_show_extra_k = 0;
%val.matlab_solver = 23;
end