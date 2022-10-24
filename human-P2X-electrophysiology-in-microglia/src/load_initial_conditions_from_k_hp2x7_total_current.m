%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [x0_new] = load_initial_conditions_from_k_hp2x7_total_current(K, opt)
    x0_new = opt.x0;
   
    x0_new(1) = 1;%K(opt.initial_condition_index + 0);    % C
    x0_new(2) = 0;%K(opt.initial_condition_index + 1);    % S
    x0_new(3) = 0;%K(opt.initial_condition_index + 2);    % D
    x0_new(4) = 0;%K(opt.initial_condition_index + 3);    % O
end
%--------------------%