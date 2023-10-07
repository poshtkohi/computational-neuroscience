%------Functions -------%
% g(x)=a*exp(-(x-b)^2/2c^2)
function [g] = gaussian_function(a, b, c, x)
    g = a * exp(-(x-b).^2/(2*c^2));
end

%--------------------%