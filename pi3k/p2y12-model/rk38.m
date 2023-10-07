% Runge-Kutta 3/8 method
function [t, y] = rk38(odefun, y0, ts, tf, dt)
t = ts:dt:tf;
y = zeros(length(y0),length(t));
y(:,1) = y0;

for k = 2 : length(t)
    yn = y(:,k-1);
    tn = t(k-1);
    k1 = dt * odefun(tn, yn);
    k2 = dt * odefun(tn + dt / 3, yn + k1 / 3);
    k3 = dt * odefun(tn + dt * 2 / 3, yn + - k1 / 3 + k2);
    k4 = dt * odefun(tn + dt, yn + k1 - k2 + k3);
    
    y(:,k) = yn + 1 / 8 * k1 + 3 / 8 * k2 + 3 / 8 * k3 + 1 / 8 * k4;
end
end
