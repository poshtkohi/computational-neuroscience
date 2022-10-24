%------Functions -------%
function [I_Na_Leak] = NA_Leak_Current_Fit(NAi, g)
    NAx = 145e-3;   % M
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    Vm = -60e-3;
    Z_Na = 1;
    E =  R * T * Z_Na^-1 * F^-1 * log(NAx/NAi);
    %g = 0.05;
    %g = 0.013954; % s−1
    I_Na_Leak = g * (Vm - E);
    %I_Ca_Leak = 0;
end
%--------------------%