%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
%%{
function [I_Na_NCX] = NA_NCX_Current(NAi, CAi, Vm, opt)
    CAx = 2e-3;     % M
    NAx = 130e-3;   % M
    %%NAi = 8e-3;   % M
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    %Vm = -60e-3;
    I_NCX_BAR1 = 50;%25;%100%30, 1; % A/m2
    I_NCX_BAR1 = I_NCX_BAR1 + I_NCX_BAR1 * opt.SA(opt.SAncx);
    I_NCX_BAR2 = 1;
    %%I_NCX_BAR = 100;%25;%100%30, 1; % A/m2
    %%I_NCX_BAR = 1.97477417445101;
    %I_NCX_BAR = 1.02531742582684;
    gamma = 0.5; % NCX partition parameter
    
    beta1 = I_NCX_BAR1 * (NAi / NAx)^3 * exp(gamma * F * Vm * R^-1 * T^-1);
    beta2 = I_NCX_BAR2 * (CAi / CAx) * exp((gamma - 1) * F * Vm * R^-1 * T^-1);
    I_Na_NCX = (beta1 - beta2);
    %%I_Ca_NCX = -2 * I_Na_NCX * 3^-1;
    
%%%%    %I_Ca_NCX = 0;
    
%%%%    %I_Ca_NCX = abs(1000 * I_Ca_NCX);
%%%%    %I_Ca_NCX = abs(I_Ca_NCX);
end
%%}
%{
function [I_NCX] = NA_NCX_Current(NAi, CAi, Vm)
    % NCX model parameters
    CAx = 2e-3;     % M
    NAx = 130e-3;   % M
    %NAi = 8.3e-3;   % M
    K_NCX_CAi = 0.502e-3; % M
    n_NCX = 3;
    n_NCX_h = 1.5;
    d_NCX = 3.04e-4;
    z_NCX = 1;
    %Vm = -60e-3;
    %gamma = 0.483;
    gamma = 0.483;
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    %g_NCX = 1.99e-9; % NCX conductance; unit: s
    g_NCX = 1; % NCX conductance; unit: s
    
   alpha = Vm * (n_NCX - 2) * z_NCX * F * R^-1 * T^-1;
   Phi_F = exp(gamma * alpha);
   Phi_R = exp((gamma - 1) * alpha);
    
   beta1 = 1 / (1 + (K_NCX_CAi / CAi)^n_NCX_h);
   beta2 = NAi^n_NCX * CAx * Phi_F - NAx^n_NCX * CAi * Phi_R;
   beta3 = 1 + d_NCX * (NAx^n_NCX * CAi + NAi^n_NCX * CAx);
    
   I_NCX = beta1 * g_NCX * beta2 * beta3^-1;
end
%}
%--------------------%