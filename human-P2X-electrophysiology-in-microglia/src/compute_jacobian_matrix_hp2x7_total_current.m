%------Functions -------%
function [J] = compute_jacobian_matrix_hp2x7_total_current(t, x)
    global my_k my_ATP;
    k = my_k;
    A = my_ATP;
    
    J(1, 1) = -(k(1) * A) + k(10);
    J(1, 2) = k(3);
    J(1, 3) = k(8);
    J(1, 4) = k(9);

    J(2, 1) = k(1) * A;
    J(2, 2) = -(k(5) + k(3) + k(2) * A);
    J(2, 3) = k(6) * A;
    J(2, 4) = k(4);

    J(3, 1) = 0;
    J(3, 2) = k(5);
    J(3, 3) = -(k(8) + k(6) * A);
    J(3, 4) = k(7);

    J(4, 1) = 0;
    J(4, 2) = k(2) * A;
    J(4, 3) = 0;
    J(4, 4) = (k(10) - k(4) - k(7) - k(9));
    %A
end
%--------------------%