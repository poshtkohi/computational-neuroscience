%------Functions -------%
function [K_new, num_added] = load_initial_conditions_from_file_hp2y12(K)
    i = load('data/initial_guess_x0_hp2y12.dat'); % Loads non-zero initial conditions
    num_added = length(i);
    K_new = [K; i];
%--------------------%

%--------------------%