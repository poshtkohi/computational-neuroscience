%------Functions -------%
%{
function [A, G] = periodic_function(t, delay, width, amplitude)
    global t1 t2 g enabled;
    
    % update phase
    if t >= t2
        t1 = t2 + delay;
        t2 = t1 + width;
        enabled = 1;
        %g = 2 * g;
    end
    
    if t > t1 && enabled == 1
        g = 2 * g;
        enabled = 0;
    end
    
    if t >= t1 && t <= t2
        A = amplitude;
    else
        A = 0.0;
    end
    
    G = g;
end
%}
%----------------------%
function [A, G] = periodic_function(t, t1, t2, delay, amplitude)
    % https://math.stackexchange.com/questions/1821507/pulse-wave-formula
    tau = t2 - t1; % width of the signal
    T = tau + delay; % delay between
    t = t - t1;
    A = t / T - floor(t / T) < tau / T;
    %G = (1 + floor(t / T)) / (1 + ceil(t / T));
    %G = 2 * (1 + floor(t / T)) / (1 + ceil(t / T));
    G = 1.5;
    %G = G^0.5;
    A = A * amplitude;
    %G = 1;
end