%------Functions -------%
function [x0_new] = load_initial_conditions_from_k_hp2y12(K, opt)
    x0_new = opt.x0;
    x0_new(1) = 0.0;%K(opt.initial_condition_index + 0);       % IP3ATP
    x0_new(2) = 0.0;%K(opt.initial_condition_index + 1);       % IP3
    x0_new(3) = 1.0;%K(opt.initial_condition_index + 2);       % C
    x0_new(4) = 0;%K(opt.initial_condition_index + 3);       % S
    x0_new(5) = 0.0;%K(opt.initial_condition_index + 4);       % D
    x0_new(6) = 0;%K(opt.initial_condition_index + 5);       % O
    x0_new(7) = 50;%K(opt.initial_condition_index + 6);      % CAs
    x0_new(8) = 0.1;%K(opt.initial_condition_index + 7);   % CAi
end
%--------------------%

%--------------------%