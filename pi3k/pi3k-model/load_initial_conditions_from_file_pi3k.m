%------Functions -------%
function [K_new, num_added] = load_initial_conditions_from_file_pi3k(K)
    i = load('data/initial-guess-x0-pi3k.dat'); % Loads non-zero initial conditions
    num_added = length(i);
    K_new = [K; i];
%--------------------%

%--------------------%