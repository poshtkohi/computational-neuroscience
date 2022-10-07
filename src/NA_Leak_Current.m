%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info
%------Functions -------%
function [I_NA_Leak] = NA_Leak_Current(NAi, Vm, opt)
    NAx = 130e-3;   % M
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    %Vm = -60e-3;
    Z_Na = 1;
    %INA_LEAK_p_soma =  R * T * Z_Na^-1 * F^-1 * log(NAx/NAi) * g*A_p;
    
    %INA_LEAK_x_p =  R * T * Z_Na^-1 * F^-1 * log(NAx/NAi) * g*A_s;
    
    E = R * T * Z_Na^-1 * F^-1 * log(NAx/NAi);
    
    %g =  0.027670538220688;
    g = opt.g_NA_Leak;
    I_NA_Leak = g * (Vm - E);
    %I_NA_Leak = 0;
end
%--------------------%