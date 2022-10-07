%------Functions -------%
function [J] = compute_jacobian_matrix_p2x(t, x)
    global my_k_p2x7 my_k_p2x4 my_ATP J;
    k_p2x7 = my_k_p2x7;
	k_p2x4 = my_k_p2x4;
    ATP = my_ATP;
	
    %disp('compute_jacobian_matrix_p2x');
    
    %J = zeros(73, 73);
    
    C1 = x(1);
    S1 = x(2);
    O1 = x(3);
    C2 = x(4);
    S2 = x(5);
    O2 = x(6);
    NAi = x(7);
    CAi = x(8);
    Vm = x(9);
	
	MtoL = 1e3;
	cytosol_coeff = 1.0; 
	Vol_p = MtoL * 123.4 * 1e-18;
	Vol_b = MtoL * 43.5769 * 1e-18;
	Vol_tot = (Vol_p + Vol_b) * cytosol_coeff;
	s_p = 243.53 * 1e-12;
	s_b = 74.3 * 1e-12;
	s_tot = s_p + s_b;
	Cm = 12e-12;
	E = 0;
	Sq = s_p;
	Vq = Vol_p;
	R = 8.314;
	T = 310;
	F = 96485.33212;
	gamma = 0.5;
	n_PMCA = 2;
	K_PMCA_CAi = 0.1e-6;
	CAx = 2e-3;
	Z_CA = 2;
    NAx = 130e-3;
	Z_NA = 1;
    I_NCX_BAR = 300;
    I_PMCA_Bar = 0.04;
	g_na_leak = 0.16859513252379;
    g_ca_leak = -0.0412872846808662;
	d1 = 180;
    d2 = 180;
    
    k_p2x7(2) = k_p2x7(1);
    k_p2x7(4) = k_p2x7(3);
    
    k_p2x4(2) = k_p2x4(1);
    k_p2x4(4) = k_p2x4(3);
    
     J(1, 1) = -ATP*k_p2x7(1);
     J(1, 2) = k_p2x7(3);
     J(1, 3) = 2*O1*k_p2x7(5);
     J(2, 1) = ATP*k_p2x7(1);
     J(2, 2) = - k_p2x7(3) - k_p2x7(6) - ATP*k_p2x7(7);
     J(2, 3) = k_p2x7(4);
     J(3, 2) = ATP*k_p2x7(2);
     J(3, 3) = - k_p2x7(4) - 4*O1^3*k_p2x7(5);
     J(4, 4) = - ATP*k_p2x4(1) - O2*k_p2x4(7);
     J(4, 5) = k_p2x4(3);
     J(4, 6) = k_p2x4(5) - C2*k_p2x4(7);
     J(5, 4) = ATP*k_p2x4(1);
     J(5, 5) = - k_p2x4(2) - k_p2x4(3) - k_p2x4(6);
     J(5, 6) = k_p2x4(4);
     J(6, 5) = k_p2x4(2);
     J(6, 6) = - k_p2x4(4) - k_p2x4(5)*(2*O2 + 1);
     J(7, 3) = (5267608507340082536959730371314977865728*((7699280930566099*log(1/(500*CAi)))/230584300921369395200000000 - Vm/400000000 + (7699280930566099*log(13/(100*NAi)))/115292150460684697600000000))/(221223185303294868540385555130685*d1);
     J(7, 6) = (5409976304835760443364047408377544835072*((53894966513962693*log(1/(500*CAi)))/576460752303423488000000000 - (7*Vm)/1000000000 + (53894966513962693*log(13/(100*NAi)))/288230376151711744000000000))/(221223185303294868540385555130685*d2);
     J(7, 7) = - (18133944210647822582035090624495*g_na_leak)/(33189892398032426259080586417143808*NAi) - (137411775650604940397629102716173312*O1)/(86415306759099558023588107472923828125*NAi*d1) - (1975758503949238602474018449864978432*O2)/(432076533795497790117940537364619140625*NAi*d2) - (147204852425375312500*I_NCX_BAR*NAi^2*exp((2673557070961925*Vm)/142833432395776))/5270537594212508423;
     J(7, 8) = (294409704850750625*exp(-(2673557070961925*Vm)/142833432395776))/28787642754005508 - (68705887825302470198814551358086656*O1)/(86415306759099558023588107472923828125*CAi*d1) - (987879251974619301237009224932489216*O2)/(432076533795497790117940537364619140625*CAi*d2);
     J(7, 9) = - (2355277638806005*g_na_leak)/115150571016022032 - (5144148932949299352499736690737283072*O1)/(86415306759099558023588107472923828125*d1) - (73964519792676412311617835661411745792*O2)/(432076533795497790117940537364619140625*d2) - (787121148163537683614227044953125*CAi*exp(-(2673557070961925*Vm)/142833432395776))/4111837825137996553152739934208 - (98390143520442210451778380619140625*I_NCX_BAR*NAi^3*exp((2673557070961925*Vm)/142833432395776))/564606731364261151704785602215936;
     J(8, 3) = (142367797495677906404317037062566969344*((7699280930566099*log(1/(500*CAi)))/230584300921369395200000000 - Vm/400000000 + (7699280930566099*log(13/(100*NAi)))/115292150460684697600000000))/(44244637060658973708077111026137*d1);
     J(8, 6) = (569471189982711625617268148250267877376*((53894966513962693*log(1/(500*CAi)))/576460752303423488000000000 - (7*Vm)/1000000000 + (53894966513962693*log(13/(100*NAi)))/288230376151711744000000000))/(221223185303294868540385555130685*d2);
     J(8, 7) = (147204852425375312500*I_NCX_BAR*NAi^2*exp((2673557070961925*Vm)/142833432395776))/15811612782637525269 - (207974579363077747628844047354208256*O2)/(432076533795497790117940537364619140625*NAi*d2) - (3713831774340674064800786559896576*O1)/(17283061351819911604717621494584765625*NAi*d1);
     J(8, 8) = (2355277638806005*CAi^3*I_PMCA_Bar)/(115150571016022032*(CAi^2 + 1/100000000000000)^2) - (18133944210647822582035090624495*g_ca_leak)/(132759569592129705036322345668575232*CAi) - (1856915887170337032400393279948288*O1)/(17283061351819911604717621494584765625*CAi*d1) - (103987289681538873814422023677104128*O2)/(432076533795497790117940537364619140625*CAi*d2) - (2355277638806005*CAi*I_PMCA_Bar)/(115150571016022032*(CAi^2 + 1/100000000000000)) - (294409704850750625*exp(-(2673557070961925*Vm)/142833432395776))/86362928262016524;
     J(8, 9) = (787121148163537683614227044953125*CAi*exp(-(2673557070961925*Vm)/142833432395776))/12335513475413989659458219802624 - (139031052241872955472965856506413056*O1)/(17283061351819911604717621494584765625*d1) - (7785738925544885506486087964359131136*O2)/(432076533795497790117940537364619140625*d2) - (2355277638806005*g_ca_leak)/230301142032044064 + (98390143520442210451778380619140625*I_NCX_BAR*NAi^3*exp((2673557070961925*Vm)/142833432395776))/1693820194092783455114356806647808;
     J(9, 3) = (256682666666666651*((7699280930566099*log(1/(500*CAi)))/230584300921369395200000000 - Vm/400000000 + (7699280930566099*log(13/(100*NAi)))/115292150460684697600000000))/(6553600*d1);
     J(9, 6) = (125610666666666659*((53894966513962693*log(1/(500*CAi)))/576460752303423488000000000 - (7*Vm)/1000000000 + (53894966513962693*log(13/(100*NAi)))/288230376151711744000000000))/(3276800*d2);
     J(9, 7) = - (258501682823733443159904278862341824057604111927*g_na_leak)/(365375409332725729550921208179070754913983135744*NAi) - (1976271960673521013628598754464449*O1)/(755578637259143234191360000000000*NAi*d1) - (6769782673796529429663923392952687*O2)/(944473296573929042739200000000000*NAi*d2) - (524605976914763146647258895990203125*I_NCX_BAR*NAi^2*exp((2673557070961925*Vm)/142833432395776))/43516068260959687423254014722048;
     J(9, 8) = (4196847815318105173178071167921625*exp(-(2673557070961925*Vm)/142833432395776))/950737950171172051122527404032 - (258501682823733443159904278862341824057604111927*g_ca_leak)/(730750818665451459101842416358141509827966271488*CAi) - (1976271960673521013628598754464449*O1)/(1511157274518286468382720000000000*CAi*d1) - (6769782673796529429663923392952687*O2)/(1888946593147858085478400000000000*CAi*d2) - (33574782522544841385424569343373*CAi*I_PMCA_Bar)/(633825300114114700748351602688*(CAi^2 + 1/100000000000000)) + (33574782522544841385424569343373*CAi^3*I_PMCA_Bar)/(633825300114114700748351602688*(CAi^2 + 1/100000000000000)^2);
     J(9, 9) = - (33574782522544841385424569343373*g_ca_leak)/1267650600228229401496703205376 - (33574782522544841385424569343373*g_na_leak)/1267650600228229401496703205376 - (256682666666666651*O1)/(2621440000000000*d1) - (879274666666666613*O2)/(3276800000000000*d2) - (11220512152394827219503674858669333832841759128125*CAi*exp(-(2673557070961925*Vm)/142833432395776))/135797164731872754491255831337010603682168832 - (1402564019049353402437959357333666729105219891015625*I_NCX_BAR*NAi^3*exp((2673557070961925*Vm)/142833432395776))/18646648182245277601080566340463268518107807744;
end
%--------------------%