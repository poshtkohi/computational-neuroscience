function J_LEAK = leak_current(CAs, CAi, alpha_leak)
    %J_LEAK = k(22) * CAi - k(21) * CAs;
    c0 = 2;
    c1 = 0.125;
    J_LEAK = alpha_leak * (c0 - (1 + c1) * CAi);
end