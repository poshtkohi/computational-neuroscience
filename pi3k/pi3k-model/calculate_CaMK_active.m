function CaMK_active = calculate_CaMK_active(CAi, CaMK_trap, k, opt)
    %alpha_CaMK = k(1);
	%beta_CaMK = k(2);
	CaMK0 = k(3);
	Km_CaM = k(4);
    
	%i = (1 + Km_CaM ./ CA_ss).^-1;
	%CaMK_bound = CaMK0 * (1 - CaMK_trap) * i;
	%CaMK_active = CaMK_bound + CaMK_trap;
    CaMK_bound = zeros(length(CaMK_trap), 1);
    CaMK_active = zeros(length(CaMK_trap), 1);
    for i=1:1:length(CaMK_trap)
	    %%CA_ss = CAi(i) - opt.CAiB;% * 1e-3;
        CA_ss = CAi(i);
	    %%ii = (1 + Km_CaM / CA_ss)^-1;
        ii = Hill(CA_ss, Km_CaM, 1)^2;%CA_ss / (CA_ss + Km_CaM);
        CaMK_bound(i) = CaMK0 * (1 - CaMK_trap(i)) * ii;
        CaMK_active(i) = CaMK_bound(i) + CaMK_trap(i);
    end
end