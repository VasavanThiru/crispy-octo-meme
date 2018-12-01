function [d_x, d_l] = pq(Q, A, g, b)
    d_l = -linsolve(A * Q^(-1) * A', A * Q^(-1) * g + b);
    d_x = -linsolve(Q, A' * d_l + g);
end
