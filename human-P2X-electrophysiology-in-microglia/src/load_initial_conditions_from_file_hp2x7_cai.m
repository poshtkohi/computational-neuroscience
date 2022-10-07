%------Functions -------%
function [K_new, num_added] = load_initial_conditions_from_file_hp2x7_cai(K)
    x0 = load('data/initial-guess-x0-hp2x7-cai.dat'); % Loads non-zero initial conditions
    num_added = length(x0);
    K_new = [K; x0];
end
%--------------------%