%------Functions -------%
function [dxdt] = reaction_network_hp2x7_cai_prediction(t, x, opt)
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
    
    % Human microglia morphology
    MtoL = 1e3; % this is volume scale convert from cubic meter to liter
    %%%%%%% important
    cytosol_coeff = 1; % [L] The volume of cytoplasm in microglia is 0.2 times of the whole cell volume
    Vol_p = 123.4 * 1e-18; % the volume of processes in M^3
    Vol_b = 43.5769 * 1e-18; % the volume of cell body in M^3
    Vol_tot = MtoL * (Vol_p + Vol_b) * cytosol_coeff;
    s_p = 243.53 * 1e-12; % the area of processes
    s_b = 74.3 * 1e-12; % the area of cell body
    s_tot = s_p + s_b;
    volume = Vol_tot;
    area = s_tot;
    
    % Cell Dimensions
    %rMiG = 3.68e-6; % activated: 3.68 uM -> Davis et al (2017)
    %area = 4 * pi* rMiG^2; % surface area of microglia in [sq. meter]
    %volume = MtoL *(4/3) * pi * rMiG^3 * cytosol_coeff; % volume of microglia in [L]
    %s_tot = area;
    %pause
    
    %str = {sprinT_P2Xf('@time %g  ATP %g', t, ATP) }; disp(str);
    % hP2X receptor
	C = x(1);
	S = x(2);
	O = x(3);
    CAi = x(4);
    
    % Current calculations
    alpha = (2.5 / 2.416) * 145 * 1e-12;
    i_hp2x7_tot = -alpha * O;
    C_mem = 12e-12; % F
    
    Z_ca = 2;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    
    %fprintf('t=%g\tCAi=%gnM\tI_P2X7_CAi=%g\tI_NCX=%g\tI_PMCA=%g\tI_LEAK=%g\n', t, CAi * 1e9, 0.1*i_hp2x7_tot, NCX_Current(CAi)*s_tot, PMCA_Current(CAi)*s_tot, Leak_Current(CAi)*s_tot);
    
    %i_cai_tot = -(0.1 * (i_hp2x7_tot/s_tot)*s_p + NCX_Current(CAi, s_p) + PMCA_Current(CAi, s_p));
   
    i_cai_tot = -(0.1 * i_hp2x7_tot + NCX_Current(CAi) * s_tot + PMCA_Current(CAi) * s_tot + Leak_Current(CAi) * s_tot);
    
    
    %alpha = 2.44;
    %CAx = 2e-3;     % M
    %J_leak = 0;%-alpha  * (CAx - CAi);
    
    %i_eflux = 0.0004 * CAi;
    %i_cai_tot = 0.1 * i_hp2x7_tot - i_eflux;
    
    dCdt = k(3) * S - k(1) * ATP * C + k(5) * O;
	dSdt = k(1) * ATP * C + k(4) * O - k(3) * S - k(2) * S;
	dOdt = k(2) * S - k(4) * O - k(5) * O^2;
    
    alpha = 1;
    dCAidt = alpha * i_cai_tot / (Z_ca * F * volume);
	
    dxdt(1) = dCdt;
	dxdt(2) = dSdt;
	dxdt(3) = dOdt;
    dxdt(4) = dCAidt;
    
    dxdt = dxdt';
end

%--------------------%