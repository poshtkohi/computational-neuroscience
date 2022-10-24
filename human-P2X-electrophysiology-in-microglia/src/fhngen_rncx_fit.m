%	All rights reserved to Alireza Poshtkohi (c) 2019-2022.
%	Email: arp@poshtkohi.info
%	Website: http://www.poshtkohi.info

%------Functions -------%
function [err] = fhngen_rncx_fit(x, opt, w_experimental)
    global rt_fig;
    
    c = Constants; % Gets constants for the model
    k = x;
    v = opt.v;
    
    err = 0.0;
	
    CAi = 45e-9;  % M
    NAi = 8e-3;   %M
    
    %k(1) = 0.502e-3;    % K_NCX_CAi
    %k(2) = 3.04e-4;     % d_NCX
    %k(3) = 0.483;       % gamma
    %k(4) = 1.99;        %g_NCX 
    
    for i=1:1:length(v)
        I_NA_NCX = NA_NCX_Current_Fit(k, v(i));
        I_CA_NCX = -2 * 3^-1 * I_NA_NCX;
        s_p = 243.53 * 1e-12; % the area of processes
        s_b = 74.3 * 1e-12; % the area of cell body
        s_tot = s_p + s_b;
        w(i) = (I_CA_NCX + I_CA_NCX);
    end
	err = loss_function(w, w_experimental, 2);%c.loss_function);
    
    figure(rt_fig);
    clf(rt_fig, 'reset');
    hold on;
    plot(v, w, '-.r*', 'MarkerIndices', 1:30:length(w));
    plot(v, w_experimental, '-bo', 'MarkerIndices', 1:30:length(w_experimental));
    xlabel('Time (s)');
    ylabel('I_{tot} (pA)'); 
    %ylabel('Normalised [O]');
    %xlim([0.0 8.5]);
    %ylim([0.0 0.05]);
    legend('Fitted', 'Experimental', 'FontWeight', 'bold');
    %drawnow update;
    hold off;

    drawnow update;
end
%--------------------%