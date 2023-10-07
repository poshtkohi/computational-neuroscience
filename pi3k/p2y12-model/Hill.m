%% This function implements the Hill function where:
%                 x^2 
%     Hill =   ---------
%              X^2 + K^n
% Sensitivity list: (x,K,n)

function out = Hill(x,K,n)

out = (x^n)/((x^n) + (K^n));

end