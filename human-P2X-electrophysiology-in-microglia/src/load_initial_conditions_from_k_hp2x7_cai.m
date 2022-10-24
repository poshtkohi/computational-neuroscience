%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [x0_new] = load_initial_conditions_from_k_hp2x7_cai(K, opt)
    x0_new = opt.x0;
   
    %x0_new(1) = K(opt.initial_condition_index + 0);		% s1
    %x0_new(5) = K(opt.initial_condition_index + 1);		% CAi
	%x0_new(6) = K(opt.initial_condition_index + 0);		% NCX
	%x0_new(8) = K(opt.initial_condition_index + 1);	% PMCA
end
%--------------------%