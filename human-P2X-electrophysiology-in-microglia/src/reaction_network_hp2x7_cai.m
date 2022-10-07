%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [dxdt] = reaction_network_hp2x7_cai(t, x, opt)
    %-- Simulation Settings --%
    k = opt.K;
    if t > opt.ATP_t1 && t < opt.ATP_t2
        ATP = opt.ATP;
    else
        ATP = 0.0;
    end
    
    
    %str = {sprinT_P2Xf('@time %g  ATP %g', t, ATP) }; disp(str);
    % hP2X receptor
	C = x(1);
	S = x(2);
	D = x(3);
	O = x(4);
	CAi = x(5);
    NCX = x(6);
    NCX3_CA = x(7);
	PMCA = x(8);
	PMCA_CA = x(9);
	
    dCdt = k(9) * O + k(8) * D + k(3) * S - k(1) * ATP * C + k(10) * O;
	dSdt = k(1) * ATP * C + k(4) * O + k(6) * D - k(5) * S - k(3) * S - k(2) * ATP * S;% + k(10) * S;
	dDdt = k(7) * O + k(5) * S - k(8) * D - k(6) * D;% + k(10) * D;
	dOdt =  k(2) * ATP * S - k(4) * O - k(7) * O - k(9) * O - k(10) * O^2;
    
    %y = O * S;
    y = O;
    
	%dCAidt =  k(14) * (k(11) * y - k(12) * CAi - k(13) * CAi^2);
    dCAidt =  k(11) * y - k(12)  * CAi - k(13) * CAi * NCX^3 + k(14) * NCX3_CA - k(15) * CAi * PMCA + k(16) * PMCA_CA ;
    
    % NCX
    dNCXdt =  k(14) * NCX3_CA - 3 * k(13) * CAi * NCX^3;
	dNCX3_CAdt = k(13) * CAi * NCX^3 - k(14) * NCX3_CA;
	% PMCA
	dPMCAdt =  k(16) * PMCA_CA - k(15) * CAi * PMCA;
	dPMCA_CAdt = k(15) * CAi * PMCA - k(16) * PMCA_CA;
	
    dxdt(1) = dCdt;
	dxdt(2) = dSdt;
	dxdt(3) = dDdt;
	dxdt(4) = dOdt;
	dxdt(5) = dCAidt;
    dxdt(6) = dNCXdt;
    dxdt(7) = dNCX3_CAdt;
	dxdt(8) = dPMCAdt;
    dxdt(9) = dPMCA_CAdt;
    
    dxdt = dxdt';
    
    %t
end

%--------------------%