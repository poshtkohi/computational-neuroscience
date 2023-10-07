%------Functions -------%
function [dxdt] = reaction_network_p2x_cai_prediction(t, x, opt)
    %-- Simulation Settings --%
    k_p2x7 = opt.K_p2x7;
    k_p2x4 = opt.K_p2x4;
    
    c = Constants; % Gets constants for the model
    
    if c.periodic == 0
        if t > opt.ATP_t1 && t < opt.ATP_t2
            ATP = opt.ATP;
        else
            ATP = 0.0;
        end
        G = 1;
    else
        %opt.amplitude = 50;
        %[ATP, G] = periodic_function(t, opt.delay, opt.width, opt.amplitude);
        [ATP, G] = periodic_function(t, opt.ATP_t1, opt.ATP_t2, opt.delay, opt.amplitude);
    end
    
    % Human microglia morphology
    MtoL = 1e3; % this is volume scale convert from cubic meter to liter
    %%%%%%% important
    cytosol_coeff = 1.0; % [L] The volume of cytoplasm in microglia is 0.2 times of the whole cell volume
    Vol_p = MtoL * 123.4 * 1e-18; % the volume of processes in M^3
    Vol_b = MtoL * 43.5769 * 1e-18; % the volume of cell body in M^3
    Vol_tot = (Vol_p + Vol_b) * cytosol_coeff;
    s_p = 243.53 * 1e-12; % the area of processes
    s_b = 74.3 * 1e-12; % the area of cell body
    s_tot = s_p + s_b;
    %volume = Vol_tot;
    %area = s_tot;
    r_p = 0.5278e-6;
    ss = pi * r_p^2;
    % Cell Dimensions
    rMiG = 3.68e-6; % activated: 3.68 uM -> Davis et al (2017)
    area = 4 * pi* rMiG^2; % surface area of microglia in [sq. meter]
    volume = MtoL *(4/3) * pi * rMiG^3 * cytosol_coeff; % volume of microglia in [L]
    %s_tot = area;
    %pause
     
    %str = {sprinT_P2Xf('@time %g  ATP %g', t, ATP) }; disp(str);
    % hP2X receptor
	C1 = x(1);
	S1 = x(2);
	D1 = x(3);
    O1 = x(4);
    C2 = x(5);
	S2 = x(6);
    D2 = x(7);
	O2 = x(8);
    NAi = x(9);
    CAi = x(10);
    Vm = x(11);
    
    % Current calculations
    R = 8.314; % Ideal gas constant; unit: J/K.mol;
    T = 310; % Absolute temperature; unit: K;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
    NAx = 130e-3;   % M
    CAx = 2e-3;   % M
    Z_Na = 1;
    Z_Ca = 2;
    %E =  R * T * Z_Na^-1 * F^-1 * log((NAx + CAx)/(2 * NAi + CAi));
    E_Na =  R * T * Z_Na^-1 * F^-1 * log(NAx/NAi);
    E_Ca =  R * T * Z_Ca^-1 * F^-1 * log(CAx/CAi);
    %E = E_Na + E_Ca;
    %Vm = -60e-3;
    %E = 0;
    %alpha1 = 2.5 * (Vm - E) * 1e3 * 1e-12;
    %i_hp2x7_tot = alpha1 * O1;
    %alpha2 = 7.0 * (Vm - E) * 1e3 * 1e-12;
    %i_rp2x4_tot = alpha2 * O2;
    %[Vm E_Na E_Ca]
    
    %pause
    %E_Na = 0;
    %E_Ca = 0;
    %E_Na = E_Na * 0.1;
    %E_Ca = 0;
    %E_Ca = E_Ca * 0.1;
    %O1 = 0;
    %O2 = 0;
    g_p2x7 = 2.417 * G;
    I_hP2X7_Na = (0.12 * g_p2x7 + 0.12 * g_p2x7 * opt.SA(opt.SAg_p2x7)) * O1 * (Vm - E_Na) * 1e3 * 1e-12;
    I_hP2X7_Ca = (0.08 * g_p2x7 + 0.08 * g_p2x7 * opt.SA(opt.SAg_p2x7)) * O1 * (Vm - E_Ca) * 1e3 * 1e-12;
    %%I_hP2X7_Na = 0.37 * (145/134) * O1 * (Vm - E_Na) * 1e3 * 1e-12;
    %%I_hP2X7_Ca = 0.1 * (145/202) * O1 * (Vm - E_Ca) * 1e3 * 1e-12;
    
    %%I_rP2X4_Na = 0.38 * (372/134) * O2 * (Vm - E_Na) * 1e3 * 1e-12;
    %%I_rP2X4_Ca = 0.08 * (372/202) * O2 * (Vm - E_Ca) * 1e3 * 1e-12;
    g_p2x4 = k_p2x4(22) * G;%6.2;
    %h2 = 1 - D2;
    %0.12 * g_p2x4 / 25
    %0.08 * g_p2x4 / 25
    %g_p2x4
    I_rP2X4_Na = (0.12 * g_p2x4 + 0.12 * g_p2x4 * opt.SA(opt.SAg_p2x7)) * O2 * (Vm - E_Na) * 1e0 * 1e-12;
    I_rP2X4_Ca = (0.08 * g_p2x4 + 0.08 * g_p2x4 * opt.SA(opt.SAg_p2x7)) * O2 * (Vm - E_Ca) * 1e0 * 1e-12;
    
    i_hp2x7_tot = I_hP2X7_Na + I_hP2X7_Ca;
    i_rp2x4_tot = I_rP2X4_Na + I_rP2X4_Ca;
    
    i_na_p2x_tot = I_hP2X7_Na + I_rP2X4_Na;
    i_ca_p2x_tot = I_hP2X7_Ca + I_rP2X4_Ca;
    
    Z_na = 1;
    Z_ca = 2;
    F = 96485.33212; % Faraday's constant; unit: C/mol;
   
    %NAi = 8e-3;%%
    %Vm = -60e-3;
    %%I_NA_NCX = 1000*NA_NCX_Current(NAi, CAi);
    %%I_CA_NCX = -2 * 3^-1 * I_NA_NCX;
    %I_NA_NCX = 3 * NCX_Current_Fit(opt.K_rncx, Vm, NAi, CAi);
    %I_CA_NCX = -200000 * NCX_Current_Fit(opt.K_rncx, Vm, NAi, CAi) * 1e-12;
    I_NA_NCX = NA_NCX_Current(NAi, CAi, Vm, opt);
    I_CA_NCX = -2 * 3^-1 * I_NA_NCX;
    
    %s_tot = s_tot^-1;
    %s_tot = 1;
    %volume = Vol_p;
    %s_tot = 1;
    
    if opt.test == false
        Sq = s_p; % s_tot
        Vq = Vol_p; % Vol_tot
    else
        %Sq = area;
        %Vq = volume;
        Sq = s_p + s_p * 0.1;
        Vq = Vol_p + Vol_p * 0.1;
    end
    %%s_tot = area;
    %%Vol_tot = volume;
    %%Sq = s_tot;
    %%Vq = Vol_tot;
    
    %d1 = d2 = 200;
    %d = 50;
    d = 40;%25
    d_na = d;
    d_ca = d;
    %d = 200;
    %d1 = d;%20; %d = 10
    %d2 = d;
    
    %%fprintf('t=%g\tCAi=%gnM\tI_P2X_CAi=%gpA\tI_CA_NCX=%gpA\tI_NA_NCX=%gpA\tI_PMCA=%gpA\tI_NA_LEAK=%gpA\tI_CA_LEAK=%gpA\n', t, CAi * 1e9, (i_hp2x7_tot*s_tot^-1*Sq/d1+i_rp2x4_tot*s_tot^-1*Sq/d2)*1e12, I_CA_NCX*Sq*1e12, I_NA_NCX*Sq*1e12, PMCA_Current(CAi)*Sq*1e12, NA_Leak_Current(NAi, Vm)*Sq*1e12, CA_Leak_Current(CAi, Vm)*Sq*1e12);
    
    %i_na_eflux = 0;%300*NAi;
    %g_soma_na = 0.1;%1e-5;
    %g_soma_ca = 0.1;
    
    i_na_eflux = 0;%-g_soma_na * R * T * Z_Na^-1 * F^-1 * log(8e-3/NAi);
    i_ca_eflux = 0;%-g_soma_ca * R * T * Z_Na^-1 * F^-1 * log(45e-9/CAi);
    
    i_nai_tot = (i_na_p2x_tot/s_tot/d_na + I_NA_NCX + NA_Leak_Current(NAi, Vm, opt));
    i_cai_tot = (i_ca_p2x_tot/s_tot/d_ca + I_CA_NCX + PMCA_Current(CAi, opt) + CA_Leak_Current(CAi, Vm, opt));
    
    %i_eflux = 0.0006 * (-CAi^3 + CAi);
    %i_cai_tot = -(0.1 * i_hp2x7_tot + 0.08 * i_rp2x4_tot + i_eflux);
    
    
    % hP2X7 reaction network
    k1f_p2x7 = k_p2x7(1);
    k2f_p2x7 = k_p2x7(2);
    k3f_p2x7 = k_p2x7(3);
    k1b_p2x7 = k_p2x7(4);
    k2b_p2x7 = k_p2x7(5);
    k3b_p2x7 = k_p2x7(6);
    k4b_p2x7 = k_p2x7(28);
    
    k1f_p2x7 = k1f_p2x7 + k1f_p2x7 * opt.SA(opt.SAk1f_p2x7);
    k2f_p2x7 = k2f_p2x7 + k2f_p2x7 * opt.SA(opt.SAk2f_p2x7);
    k3f_p2x7 = k3f_p2x7 + k3f_p2x7 * opt.SA(opt.SAk3f_p2x7);
    k1b_p2x7 = k1b_p2x7 + k1b_p2x7 * opt.SA(opt.SAk1b_p2x7);
    k2b_p2x7 = k2b_p2x7 + k2b_p2x7 * opt.SA(opt.SAk2b_p2x7);
    k3b_p2x7 = k3b_p2x7 + k3b_p2x7 * opt.SA(opt.SAk3b_p2x7);
    k4b_p2x7 = k4b_p2x7 + k4b_p2x7 * opt.SA(opt.SAk4b_p2x7);
    %ki1_p2x7 = k_p2x7(25);
    %ki2_p2x7 = k_p2x7(26);
    
    A = ATP * 1e3;
    
    ki1_p2x7 = 1;
    ki2_p2x7 = 1;
    
    i1_p2x7 = exp(-ki1_p2x7 * A);
    i2_p2x7 = exp(-ki2_p2x7 * A);
    %i3 = exp(-ki3 * A);
    
    k1b_p2x7 = k1b_p2x7 * i1_p2x7;
    k2b_p2x7 = k2b_p2x7 * i2_p2x7;
    %k3b_p2x7 = k3b_p2x7 * i3_p2x7;
    k3f_p2x7 = k3f_p2x7 * exp(A);
    k3b_p2x7 = k3b_p2x7 * exp(A);
    
    k4b_p2x7 = k4b_p2x7 * exp(-1000 * O1);
    
    dC1dt = k1b_p2x7 * S1 - k1f_p2x7 * ATP * C1 + k4b_p2x7 * D1;
	dS1dt = k3b_p2x7 * D1 + k1f_p2x7 * ATP * C1 + k2b_p2x7 * O1 - k1b_p2x7 * S1 - k2f_p2x7 * ATP * S1 - k3f_p2x7 * S1;
    dD1dt = k3f_p2x7 * S1 - k3b_p2x7 * D1 - k4b_p2x7 * D1;
	dO1dt = k2f_p2x7 * ATP * S1 - k2b_p2x7 * O1;
    
    % rP2X4
    k1f_p2x4 = k_p2x4(1);
    k2f_p2x4 = k_p2x4(2);
    k3f_p2x4 = k_p2x4(3);
    k1b_p2x4 = k_p2x4(4);
    k2b_p2x4 = k_p2x4(5);
    k3b_p2x4 = k_p2x4(6);
    k4b_p2x4 = k_p2x4(28);
    
    k1f_p2x4 = k1f_p2x4 + k1f_p2x4 * opt.SA(opt.SAk1f_p2x4);
    k2f_p2x4 = k2f_p2x4 + k2f_p2x4 * opt.SA(opt.SAk2f_p2x4);
    k3f_p2x4 = k3f_p2x4 + k3f_p2x4 * opt.SA(opt.SAk3f_p2x4);
    k1b_p2x4 = k1b_p2x4 + k1b_p2x4 * opt.SA(opt.SAk1b_p2x4);
    k2b_p2x4 = k2b_p2x4 + k2b_p2x4 * opt.SA(opt.SAk2b_p2x4);
    k3b_p2x4 = k3b_p2x4 + k3b_p2x4 * opt.SA(opt.SAk3b_p2x4);
    k4b_p2x4 = k4b_p2x4 + k4b_p2x4 * opt.SA(opt.SAk4b_p2x4);
    
    k4b_p2x4 = k4b_p2x4 * exp(-200 * O2);    
    A = ATP * 0.2e3;% * 1e3;
    k2f_p2x4 = k2f_p2x4 * exp(A);
    k1b_p2x4 = k1b_p2x4 * exp(-A);
    k2b_p2x4 = k2b_p2x4 * exp(-A);
    k3f_p2x4 = k3f_p2x4 * exp(A);
    k3b_p2x4 = k3b_p2x4 * exp(A);
    
    dC2dt = k1b_p2x4 * S2 - k1f_p2x4 * ATP * C2 + k4b_p2x4 * D2;
	dS2dt = k3b_p2x4 * D2 + k1f_p2x4 * ATP * C2 + k2b_p2x4 * O2 - k1b_p2x4 * S2 - k2f_p2x4 * ATP^0 * S2 - k3f_p2x4 * S2;
    dD2dt = k3f_p2x4 * S2 - k3b_p2x4 * D2 - k4b_p2x4 * D2;
	dO2dt = k2f_p2x4 * ATP^0 * S2 - k2b_p2x4 * O2;
    
    dNAidt = -(i_nai_tot * Sq + ss * i_na_eflux) / (Z_na * F * Vq);
    dCAidt = -(i_cai_tot * Sq + ss * i_ca_eflux) / (Z_ca * F * Vq);
    
    Cm = 12e-12;
    dVmdt = -Cm^-1 * (i_na_p2x_tot/d_na + i_ca_p2x_tot/d_ca + (I_NA_NCX + I_CA_NCX) * s_tot + PMCA_Current(CAi, opt) * s_tot + NA_Leak_Current(NAi, Vm, opt) * s_tot + CA_Leak_Current(CAi, Vm, opt) * s_tot + ss * i_na_eflux + ss * i_ca_eflux);
    %dVmdt = -Cm^-1 * (i_hp2x7_tot/d1 + i_rp2x4_tot/d2 + (I_NA_NCX + I_CA_NCX) * s_tot + PMCA_Current(CAi) * s_tot + NA_Leak_Current(NAi, Vm) * s_tot + CA_Leak_Current(CAi, Vm) * s_tot + ss * i_na_eflux + ss * i_ca_eflux);
	
    dxdt(1) = dC1dt;
	dxdt(2) = dS1dt;
    dxdt(3) = dD1dt;
	dxdt(4) = dO1dt;
    dxdt(5) = dC2dt;
	dxdt(6) = dS2dt;
    dxdt(7) = dD2dt;
	dxdt(8) = dO2dt;
    dxdt(9) = dNAidt;
    dxdt(10) = dCAidt;
    dxdt(11) = dVmdt;
    
    dxdt = dxdt';
end

%--------------------%