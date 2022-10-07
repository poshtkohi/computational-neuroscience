%------Functions -------%
function [I_Ca_Leak] = CA_Leak_Current_Fit(CAi, g)
    CAx = 2e-3;     % M
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    Vm = -60e-3;
    Z_Ca = 2;
    E =  R * T * Z_Ca^-1 * F^-1 * log(CAx/CAi);
    %g = 0.05;
    %g = 0.013954; % s−1
    I_Ca_Leak = g * (Vm - E);
    %I_Ca_Leak = 0;
end
%--------------------%