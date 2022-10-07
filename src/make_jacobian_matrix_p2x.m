%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info
%------Functions -------%
function [string, sparsity_pattern] = make_jacobian_matrix_p2x()
    string = '';
    syms k1_p2x7 k2_p2x7 k3_p2x7 k4_p2x7 k5_p2x7 k6_p2x7 k1_p2x4 k2_p2x4 k3_p2x4 k4_p2x4 k5_p2x4 k6_p2x4 k7_p2x4 real
	syms C1 S1 O1 C2 S2 O2 NAi CAi Vm real
    syms ATP real
	syms beta1 beta2 I_NA_NCX I_CA_NCX I_PMCA I_CA_Leak I_NA_Leak alpha1 alpha2 i_hp2x7_tot i_rp2x4_tot i_nai_tot i_cai_tot E real
    syms I_NCX_BAR I_PMCA_Bar g_ca_leak g_na_leak d1 d2
	
	str = {'C1', 'S1', 'O1', 'C2', 'S2' , 'O2', 'NAi', 'CAi', 'Vm'};
	
	MtoL = 1e3;
	cytosol_coeff = 1.0; 
	Vol_p = MtoL * 123.4 * 1e-18;
	Vol_b = MtoL * 43.5769 * 1e-18;
	Vol_tot = (Vol_p + Vol_b) * cytosol_coeff;
	s_p = 243.53 * 1e-12;
	s_b = 74.3 * 1e-12;
	s_tot = s_p + s_b;
	Cm = 12e-12;
	%E = 0;
	Sq = s_p;
	Vq = Vol_p;
	R = 8.314;
	T = 310;
	F = 96485.33212;
	%%I_NCX_BAR = 3000;
	gamma = 0.5;
	%%I_PMCA_Bar = 0.05;
	n_PMCA = 2;
	K_PMCA_CAi = 0.1e-6;
	CAx = 2e-3;
	Z_CA = 2;
	%%g_ca_leak =-0.0238702843361464;
	NAx = 130e-3;
	Z_NA = 1;
	%%g_na_leak = 0.147953211839816;
	%%d = 20;


    %E = R * T * Z_NA^-1 * F^-1 * log(NAx/NAi) + R * T * Z_CA^-1 * F^-1 * log(CAx/CAi);
    E = R * T * Z_NA^-1 * F^-1 * log(NAx/NAi) + R * T * Z_CA^-1 * F^-1 * log(CAx/CAi);
    
	beta1 = I_NCX_BAR * (NAi / NAx)^3 * exp(gamma * F * Vm * R^-1 * T^-1);
	beta2 = (CAi / CAx) * exp((gamma - 1) * F * Vm * R^-1 * T^-1);
	I_NA_NCX = 1 * (beta1 - beta2);
	I_CA_NCX = -2 * 3^-1 * I_NA_NCX;

	I_PMCA = I_PMCA_Bar * (CAi^n_PMCA/(CAi^n_PMCA + K_PMCA_CAi^n_PMCA));

	I_CA_Leak = g_ca_leak * (Vm - R * T * Z_CA^-1 * F^-1 * log(CAx/CAi));

	I_NA_Leak = g_na_leak * (Vm - R * T * Z_NA^-1 * F^-1 * log(NAx/NAi));



	alpha1 = 2.5 * (Vm - E) * 1e3 * 1e-12;
	i_hp2x7_tot = alpha1 * O1;
	alpha2 = 7.0 * (Vm - E) * 1e3 * 1e-12;
	i_rp2x4_tot = alpha2 * O2;
	i_nai_tot = -(0.37 * i_hp2x7_tot/s_tot/d1 + 0.38 * i_rp2x4_tot/s_tot/d2 + I_NA_NCX + I_NA_Leak);
	i_cai_tot = -(0.1 * i_hp2x7_tot/s_tot/d1 + 0.08 * i_rp2x4_tot/s_tot/d2 + I_CA_NCX + I_PMCA + I_CA_Leak);

	f(1) = k3_p2x7 * S1 - k1_p2x7 * ATP * C1 + k5_p2x7 * O1^2;
	f(2) = k1_p2x7 * ATP * C1 + k4_p2x7 * O1 - k3_p2x7 * S1 - k2_p2x7 * ATP * S1 - k6_p2x7 * S1;
	f(3) = k2_p2x7 * ATP * S1 - k4_p2x7 * O1 - k5_p2x7 * O1^4;
	f(4) = k3_p2x4 * S2 - k1_p2x4 * ATP * C2 + k5_p2x4 * O2 - k7_p2x4 * C2 * O2;
	f(5) = k1_p2x4 * ATP * C2 + k4_p2x4 * O2 - k3_p2x4 * S2 - k2_p2x4 * S2 - k6_p2x4 * S2;
	f(6) = k2_p2x4 * S2 - k4_p2x4 * O2 - k5_p2x4 * (O2^2 + O2);
	f(7) = i_nai_tot * Sq / (Z_NA * F * Vq);
	f(8) = i_cai_tot * Sq / (Z_CA * F * Vq);
	f(9) = -Cm^-1 * (0.47 * i_hp2x7_tot/d1 + 0.46 * i_rp2x4_tot/d2 + (I_NA_NCX + I_CA_NCX) * s_tot + I_PMCA * s_tot + I_NA_Leak * s_tot + I_CA_Leak * s_tot);
    
    diff(f(9), CAi)
    
    %return
    
    sparsity_pattern = zeros(9, 9);


    for i=1:1:length(f)
        %s = sprintf('dx%d/dt=%s', i, f(i)); disp(s); continue;
        for j=1:1:9
            %str(j)
            x = sym(str(j));
            s = char(diff(f(i), x));
            s = sprintf('J(%d, %d) = %s;', i, j, s);
            string = sprintf('%s\n%s', string, s);
            disp(s);
            %s = strrep(s, '_', '(');
            %s = strrep(s, '', ')');
            %if s ~=	'0'
            %    s = sprintf('J(%d, %d) = %s;', i, j, s);
            %    string = sprintf('%s\n%s', string, s);
            %    sparsity_pattern(i, j) = 1;
            %    %disp(s);
            %end
        end
    end
    
    sparsity_pattern = sparse(sparsity_pattern);
end
%--------------------%