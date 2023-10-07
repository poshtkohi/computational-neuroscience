%------Functions -------%
function [dxdt] = reaction_network_hp2y12(t, x, opt)
    %-- Simulation Settings --%
    %c = Constants; % Gets constants for the model
    ADP = 0.0;
    k = opt.K;
	
	if t > opt.ADP_t1 && t < opt.ADP_t2
        ADP = opt.ADP;
    else
        ADP = 0.0;
        %disp(t);
    end
    
    %ADP = ATP_function(t/60, opt);
    
    IP3ATP = x(1);
    IP3 = x(2);
    C = x(3);
    S = x(4);
	D = x(5);
    %h = x(5);
    O = x(6);
    CAs = x(7);
    CAi = x(8);
       
    %dCAidt
    
    % G-protein cascade and IP3 generation
    %IP3Star = k(1); 	% 0.16 microMolar
	tau_IP3 = k(2);% decay rate of IP3, 7(s) 0.016
	r_IP3 = k(3);	% rate of IP3 production 0.5uM Sec^-1 0.1
    dIP3ATPdt = - IP3ATP/tau_IP3 + r_IP3 * ADP;
    
	PLC_DELTA_PRIME_BAR = k(4);
    K_DELTA = k(5);
    K_PLC_DELTA = k(6);
    r_5P_BAR = k(7);
    v_3K_BAR = k(8);
    K_D = k(9);
    K_3 = k(10);
	
    PLC_DELTA_PRIME = PLC_DELTA_PRIME_BAR / ( 1 + IP3/K_DELTA);
    PLC_DELTA = PLC_DELTA_PRIME * Hill(CAi, K_PLC_DELTA, 2);
    IP3_5P = r_5P_BAR * IP3;
    IP3_3K = v_3K_BAR * Hill(CAi, K_D, 4) * Hill(IP3, K_3, 1);
    dIP3dt = IP3ATP + PLC_DELTA - IP3_5P - IP3_3K;
    
    % IP3R  
    alpha1 = k(11);
    alpha2 = k(12);
    alpha2_bar = k(13);
    alpha3 = k(14);
    beta1 = k(15);
    beta2 = k(16);
    beta3 = k(17);
    K1 = k(18);
    K2 = k(19);
    K3 = k(20);
    K4 = k(21);
    K5 = k(22);
    K6 = k(23);
    K7 = k(24);
    K8 = k(25);
    
    CAii = CAi;
    
    h = (CAi + 0 * K1) / (CAi + 1.0);
    alpha3 = alpha3 * h;
    beta3 = beta3 * h;
    beta4 = 10 * exp(-1000 * (CAi - 0.1));

    dCdt = beta1 * S - alpha1 * IP3 * C + beta4 * D;
    dSdt = alpha1 * IP3 * C + beta3 * D + beta2 * O  - beta1 * S - alpha3 * S - alpha2 * S * CAi - 0 * alpha2_bar * S;
    dDdt = alpha3 * S - beta3 * D - beta4 * D;
    dOdt = alpha2 * S * CAi + 0 * alpha2_bar * S - beta2 * O;
	
	% Calcium dynamics
    alpha_leak = k(26);
    alpha_ip3r = k(27);
    V_SECRA = k(28);
    K_SECRA = k(29);
    J_LEAK = leak_current(CAs, CAi, alpha_leak);
    J_IP3R = ip3r_current(O, CAs, CAi, D, alpha_ip3r);
    J_SECRA = secra_current(CAi, V_SECRA, K_SECRA);

    dCAsdt = J_SECRA - J_IP3R - J_LEAK;
    dCAidt = J_IP3R + J_LEAK - J_SECRA;
    
    dxdt(1) = dIP3ATPdt;
    dxdt(2) = dIP3dt;
    dxdt(3) = dCdt;
    dxdt(4) = dSdt;
    dxdt(5) = dDdt;
    %dxdt(5) = dhdt;
    dxdt(6) = dOdt;
    dxdt(7) = dCAsdt;
    dxdt(8) = dCAidt;
    
    dxdt = dxdt';
end
%%--------------------%