%------Functions -------%
function [x0_new] = load_initial_conditions_from_k_rp2x4_total_current(K, opt)
    x0_new = opt.x0;
   
    %x0_new(1) = K(opt.initial_condition_index + 0);    % C
    %x0_new(2) = K(opt.initial_condition_index + 1);    % S
    %x0_new(3) = K(opt.initial_condition_index + 2);    % D
    %x0_new(3) = K(opt.initial_condition_index + 1);    % O
end
%--------------------%