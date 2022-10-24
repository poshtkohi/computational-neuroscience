%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [I_PMCA] = PMCA_Current(CAi, opt)
    % PMCA model parameters
    %%I_PMCA_Bar = 10e-6;  % A
    %%n_PMCA = 2;
    %%K_PMCA_CAi = 2.5e-6;   % M
    
    %%I_PMCA = I_PMCA_Bar * (CAi^n_PMCA/(CAi^n_PMCA + K_PMCA_CAi^n_PMCA)) * S_mem;
    %I_PMCA = 0.0;
    
    
    %%I_PMCA_Bar = 0.4e-6;  % A/m^2
    %I_PMCA_Bar = 0.740341573296545;
    %%I_PMCA_Bar = 0.305418960755253;
    %%I_PMCA_Bar = 0.05;%0.05;%1 * 0.03;%0.2e-6;%0.05;%0.2e-6;
    I_PMCA_Bar = 0.06;
    I_PMCA_Bar = I_PMCA_Bar + I_PMCA_Bar * opt.SA(opt.SApmca);
    n_PMCA = 1;
    K_PMCA_CAi = 0.1e-6;   % M
    I_PMCA = I_PMCA_Bar * (CAi^n_PMCA/(CAi^n_PMCA + K_PMCA_CAi^n_PMCA));
    %I_PMCA = 0;
end

%%%function [I_PMCA] = PMCA_Current(CAi)
%%%    % PMCA model parameters
%%%    %I_PMCA_Bar = 2.67e-12;  % A
%%%    I_PMCA_Bar = 0.04e-2;  % A
%%%    n_PMCA = 1.4;
%%%    K_PMCA_CAi = 0.260e-6;   % M
   
%%%    %CAi = abs(CAi);
%%%    I_PMCA = I_PMCA_Bar * (CAi^n_PMCA/(CAi^n_PMCA + K_PMCA_CAi^n_PMCA));
%%%    %I_PMCA = 0.0;
%%%end
%--------------------%