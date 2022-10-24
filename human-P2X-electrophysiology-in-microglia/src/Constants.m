%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

% Global Parameters
function [val] = Constants
% -- Simulation Parameters --%
val.presimulation_mode = 0;
val.curefitting_mode = 1;
val.simulation_mode = 2;

val.PopSize = 18;
val.loss_function = 2;
val.loss_function_division = true;
val.resume = 'no'; % yes or no.
val.should_plot = 0;
val.should_plot_state_space = 0;
end