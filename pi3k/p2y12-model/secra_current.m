function J_SECRA = secra_current(CAi, V_SECRA, K_SECRA)
    N_SECRA = 2;
    J_SECRA = V_SECRA * CAi^N_SECRA/(K_SECRA^N_SECRA+CAi^N_SECRA);
end