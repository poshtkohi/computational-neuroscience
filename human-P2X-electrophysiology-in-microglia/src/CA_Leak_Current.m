%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [I_CA_Leak] = CA_Leak_Current(CAi, Vm, opt)
    CAx = 2e-3;     % M
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    %Vm = -60e-3;
    Z_Ca = 2;
    E =  R * T * Z_Ca^-1 * F^-1 * log(CAx/CAi);
    
    %%g =  0.0606259308788846;
    g = opt.g_CA_Leak;
    I_CA_Leak = g * (Vm - E);
    %I_CA_Leak = 0;
end
%--------------------%