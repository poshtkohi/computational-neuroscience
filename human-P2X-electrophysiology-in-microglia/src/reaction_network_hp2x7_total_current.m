%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------
function [dxdt] = reaction_network_hp2x7_total_current(t, x, opt)
    %-- Simulation Settings --%
    global nT_P2X T;
    k = opt.K;
    %if t > opt.ATP_t1 && t < opt.ATP_t2
    %    ATP = opt.ATP;
    %else
    %    ATP = 0.0;
    %end
    
    if t > (opt.ATP_t1 + (nT_P2X + 1) * T) && t < (opt.ATP_t2 + (nT_P2X + 1) * T) % T is the period in minutes
        nT_P2X = nT_P2X + 1;
    end
    if t > (opt.ATP_t1 + nT_P2X * T) && t < (opt.ATP_t2 + nT_P2X * T)
        ATP = opt.ATP;
        %disp('xxx');
    %elseif ismembertol(t, c.ATPPulseFallTime)
    else
        ATP = 0.0;
        %disp(t);
    end
    
    %str = {sprinT_P2Xf('@time %g  ATP %g', t, ATP) }; disp(str);
    % hP2X receptor
	C = x(1);
	S = x(2);
    D = x(3);
	O = x(4);
    
    %k(2) = k(1);
    %k(4) = k(3);
    
    k1f = k(1);
    k2f = k(2);
    k3f = k(3);
    k1b = k(4);
    k2b = k(5);
    k3b = k(6);
    kc = k(7:9);
    ks = k(10:12);
    kd = k(19:21);
    ko = k(22:24);
    k4b = k(28);
    
    %k1b = k1b * (ATP * 1e3 + 1e-6)^-1;
    %k2b = k2b * (ATP + 1e-6)^-1;
    %k3b = k3b * (ATP + 1e-6)^-1;
    
    %k1b = k1b * exp(-k(19) * O);
    %k2b = k2b * exp(-k(20) * O);
    %k3b = k3b * exp(-k(21) * O);
    
    a = k(19);
    b = k(20);
    %a = 10;
    %b = 10;
    %i = (O + a) / (exp((O + a) / b)- 1);
    %i = (100 - O) * exp(-20 * O);
    %i = (ATP * 1e3 + 1e-6)^-1;
    ki1 = k(25);
    ki2 = k(26);
    %k(25:26)
    %ki3 = k(27);
    
    A = ATP * 1e3;
    ki1 = 1.0;
    ki2 = 1.0;
    %i1 = exp(-ki1 * A * O * t);
    %i2 = exp(-ki2 * A * O * t);
    i1 = exp(-ki1 * A);
    i2 = exp(-ki2 * A);
    %i3 = exp(-ki3 * A);
    
    k1b = k1b * i1;
    k2b = k2b * i2;
    %k3b = k3b * i3;
    k3f = k3f * exp(A);
    k3b = k3b * exp(A);
    
    k4b = k4b * exp(-1000 * O);
    
    %ki
    
    dCdt = k1b * F(S, ks) - k1f * ATP * F(C, kc) + k4b * D;
	dSdt = k3b * F(D, kd) + k1f * ATP * F(C, kc) + k2b * F(O, ko) - k1b * F(S, ks) - k2f * ATP * F(S, ks) - k3f * F(S, ks);
    dDdt = k3f * F(S, ks) - k3b * F(D, kd) - k4b * D;
	dOdt = k2f * ATP * F(S, ks) - k2b * F(O, ko);
	
    dxdt(1) = dCdt;
	dxdt(2) = dSdt;
    dxdt(3) = dDdt;
	dxdt(4) = dOdt;
    
    dxdt = dxdt';
end

%--------------------%