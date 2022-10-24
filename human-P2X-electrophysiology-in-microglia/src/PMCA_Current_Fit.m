%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [I_PMCA] = PMCA_Current_Fit(CAi, I_PMCA_Bar)
    % PMCA model parameters 
    %%I_PMCA_Bar = 0.4e-6;  % A/m^2
    n_PMCA = 2;
    K_PMCA_CAi = 0.1e-6;   % M
    I_PMCA = I_PMCA_Bar * (CAi^n_PMCA/(CAi^n_PMCA + K_PMCA_CAi^n_PMCA));
    
    %I_PMCA = 0;
end
%--------------------%