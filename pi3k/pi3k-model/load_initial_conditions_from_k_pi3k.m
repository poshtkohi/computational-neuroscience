%------Functions -------%
function [x0_new] = load_initial_conditions_from_k_pi3k(K, opt)
    x0_new = opt.x0;
	
    x0_new(2) = 700;    % PIP2          nM
    x0_new(3) = 100;    % PI3K         nM
    x0_new(6) = 270;    % PTEN         nM
    x0_new(8) = 700;    % AKT           nM
    x0_new(11) = 150;   % PP2A         nM
	x0_new(13) = 1000;  % PDK         	nM
    
    x0_new = x0_new * 1e-3; % nM to uM
end
%--------------------%

%--------------------%