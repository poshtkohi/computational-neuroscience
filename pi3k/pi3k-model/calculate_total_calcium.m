function cai_tot = calculate_total_calcium(cai_p2y, cai_p2x, k)
    alpha = 1;
	beta = k(27);
    
    cai_tot = alpha * cai_p2y + beta * cai_p2x;
end