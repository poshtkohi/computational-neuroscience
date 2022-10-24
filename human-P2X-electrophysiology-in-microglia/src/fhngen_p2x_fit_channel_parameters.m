%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [err] = fhngen_p2x_fit_channel_parameters(x, opt)

    %g_NA_leak = x(1);
    g_CA_leak = x(1);
    I_PMCA_Bar = x(2);
    I_NCX_Bar = x(3);
    NAi = opt.NAi;
    CAi = opt.CAi;
    
    I_NA_NCX = NA_NCX_Current_Fit(NAi, CAi, I_NCX_Bar);
    %I_NA_NCX = NA_NCX_Current(NAi, CAi);
    I_CA_NCX = -2 * 3^-1 * I_NA_NCX;
    
    %i_nai_tot = -(I_NA_NCX + NA_Leak_Current_Fit(NAi, g_NA_leak));
    i_cai_tot = -(I_CA_NCX + PMCA_Current_Fit(CAi, I_PMCA_Bar) + CA_Leak_Current_Fit(CAi, g_CA_leak));
    
    %err = abs(i_nai_tot) + abs(i_cai_tot);
    err = abs(i_cai_tot);
end
%--------------------%