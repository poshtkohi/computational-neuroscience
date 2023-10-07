function J_IP3R = ip3r_current(O, CAs, CAi, D, alpha_ip3r)
    %J_IP3R = k(19) * OO * CAs - k(20) * OO * CAi;
    c0 = 2;
    c1 = 0.125;
    h = 1 - D;
    J_IP3R = alpha_ip3r * h * O^2 * (c0 - (1 + c1) * CAi);
end