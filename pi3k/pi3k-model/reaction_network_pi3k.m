%------Functions -------%
function [dxdt] = reaction_network_pi3k(t, x, opt)
    %-- Simulation Settings --%
    k = opt.K;
    %str = {sprintf('@time %g  ATP %g', t, ADP) }; disp(str);
	%-------------- State variables -------------%
    % CaMKII
    CaMK_trap = x(1);
    % PI3K pathway
	x1 = x(2);      % PIP2
	x2 = x(3);      % PI3K
	y1 = x(4);      % PIP2-PI3K
	x3 = x(5);      % PIP3
	x4 = x(6);      % PTEN
	y2 = x(7);      % PIP3-PTEN
	x5 = x(8);      % AKT
    y3 = x(9);      % PIP3-AKT
	x6 = x(10);     % pAKT
	x7 = x(11);     % PP2A
	y4 = x(12);     % pAKT-PP2A
	x8 = x(13);     % PDK1
	y5 = x(14);     % PDK1-PIP3
	y6 = x(15);     % PDK1-PIP3-Akt
    %-------------- Parameters -------------%
    % CaMKII
    alpha_CaMK = k(1);
	beta_CaMK = k(2);
	CaMK0 = k(3);
	Km_CaM = k(4);
    %alpha_CaMK = 36;
    %beta_CaMK = 1.8;
    %CaMK0 = 0.05;
    %Km_CaM = 1;
    % PI3k pathway
	a1 = k(5);
	a2 = k(6);
	a3 = k(7);
	a4 = k(8);
    a5 = k(9);
	a6 = k(10);
	a7 = k(11);
    
	b1 = k(12);
	b2 = k(13);
	b3 = k(14);
	b4 = k(15);
    b5 = k(16);
	b6 = k(17);
	b7 = k(18);
	
	c1 = k(19);
	c2 = k(20);
	c3 = k(21);
	c4 = k(22);
    c5 = k(23);
    c6 = k(24);
	%-------------- Calcium dynamics -------------%
    y_p2y_cai = deval(opt.sol_p2y_cai, t * 60);
    cai_p2y = y_p2y_cai(8, 1);
	CA_ss_p2y = cai_p2y - opt.CAiB_p2y;
	
	cai_p2x = opt.f(t * 60);
    CA_ss_p2x = cai_p2x - opt.CAiB_p2x;
    CA_ss_p2x = CA_ss_p2x * 1e6; % M to uM
	
	%CAi_tot_ss = CA_ss_p2y + CA_ss_p2x;
    CAi_tot_ss = calculate_total_calcium(CA_ss_p2y, CA_ss_p2x, k);
    %str = {sprintf('@time %g  ATP %g', t, ADP) }; disp(str);
    %-------------- Kintetic reactions -------------%
    % CaMKII
    i = Hill(CAi_tot_ss, Km_CaM, 1)^2;
	CaMK_bound = CaMK0 * (1 - CaMK_trap) * i;
	CaMK_active = CaMK_bound + CaMK_trap;
    dCaMKtrap = alpha_CaMK * CaMK_bound * CaMK_active - beta_CaMK * CaMK_trap;
    % PI3K pathway
    x2 = x2 * CaMK_active;%Hill(CaMK_active, 0.5, 1);	% PI3K
    dx1dt = -a1 * x1 * x2 + b1 * y1 + c2 * y2;
	dx2dt = -a1 * x1 * x2 + (b1 + c1) * y1;
	dy1dt = a1 * x1 * x2 - (b1 + c1) * y1;
	dx3dt = -a2 * x3 * x4 + b2 * y2 + c1 * y1 - a3 * x3 * x5 + (b3 + c3) * y3 - a5 * x3 * x8 + b5 * y5 + c6 * y5; % PIP3
	dx4dt = -a2 * x3 * x4 + (b2 + c2) * y2; % PTEN
	dy2dt = a2 * x3 * x4 - (b2 + c2) * y2;
	dx5dt = -a3 * x3 * x5 + b3 * y3 + c4 * y4 - a7 * y5 * x5 + b7 * y6; % Akt
	dy3dt = a3 * x3 * x5 - (b3 + c3) * y3 - a6 * y3 * x8 + b6 * y6; % PIP3-AKT
	dx6dt = -a4 * x6 * x7 + b4 * y4 + c3 * y3 + c5 * y6; % pAkt
	dx7dt = -a4 * x6 * x7 + (c4 + b4) * y4;
	dy4dt = a4 * x6 * x7 - (b4 + c4) * y4;
    dx8dt = - a5 * x3 * x8 + b5 * y5 - a6 * y3 * x8 + b6 * y6 + c6 * y5; % PDK1
    dy5dt = a5 * x3 * x8 - b5 * y5 - a7 * y5 * x5 + b7 * y6 + c5 * y6 - c6 * y5; % PDK1-PIP3
	dy6dt = a6 * y3 * x8 - b6 * y6 + a7 * y5 * x5 - b7 * y6 - c5 * y6; %PDk1-PIP3-Akt
    
    dxdt(1) = dCaMKtrap;
    dxdt(2) = dx1dt;
	dxdt(3) = dx2dt;
	dxdt(4) = dy1dt;
	dxdt(5) = dx3dt;
	dxdt(6) = dx4dt;
	dxdt(7) = dy2dt;
	dxdt(8) = dx5dt;
	dxdt(9) = dy3dt;
    dxdt(10) = dx6dt;
	dxdt(11) = dx7dt;
	dxdt(12) = dy4dt;
	dxdt(13) = dx8dt;
    dxdt(14) = dy5dt;
	dxdt(15) = dy6dt;
    
    dxdt = dxdt';
end

%--------------------%