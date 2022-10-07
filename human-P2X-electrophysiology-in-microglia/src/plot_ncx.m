%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

maxNumCompThreads(1);

CAx = 2e-3;   % M
CAi = 45e-9;  % M
NAi = 8e-3;   %M
NAx = 130e-3; % M

Vm = -60e-3;
gamma = 0.5;
my_na_ncx_current(NAi, CAi, Vm, gamma)
return
gamma = linspace(0, 1, 1000);  % all time points.
I = zeros(length(gamma), 1);

for i=1:1:length(gamma)
    I(i) = my_na_ncx_current(NAi, CAi, Vm, gamma(i));
end
figure;
plot(gamma, I);
legend(sprintf('Vm=%gmV', Vm * 1e3), 'FontWeight', 'bold');
xlabel('\gamma');
ylabel('I_{NA}\__{NCX}');
%--------------------%
function [I_Na_NCX] = my_na_ncx_current(NAi, CAi, Vm, gamma)
    CAx = 2e-3;     % M
    NAx = 130e-3;   % M
    %%NAi = 8e-3;   % M
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    %Vm = -60e-3;
    I_NCX_BAR = 1;%1; % A/m2
    %gamma = 0.5; % NCX partition parameter
    
    beta1 = (NAi / NAx)^3 * exp(gamma * F * Vm * R^-1 * T^-1);
    beta2 = (CAi / CAx) * exp((gamma - 1) * F * Vm * R^-1 * T^-1);
    I_Na_NCX = I_NCX_BAR * (beta1 - beta2);
end
%--------------------%