function [d_x, d_l] = pq(Q, A, g, b)
    % Q: hessien
    % A: gradient de c la contrainte
    % g: gradient de f le cout
    % b: la contrainte -c
    % d_l = -linsolve(A * Q^(-1) * A', A * Q^(-1) * g + b);
    % d_x = -linsolve(Q, A' * d_l + g);
    % d_l = -linsolve(A * inv(Q) * A', A * inv(Q) * g + b);
    % d_x = -linsolve(Q, A' * d_l + g);
    % L'erreur cachee
    if min(abs(eig(A * Q^(-1) * A'))) < 1e-12
        Q
        A
        g
        b
    end
    d_l = -inv(A * Q^(-1) * A') * (A * Q^(-1) * g + b);
    d_x = -inv(Q) * (A' * d_l + g);
    l = length(d_l);
    x = length(d_x);
    for i = 1:l
        if d_l(i) ~= d_l(i)
            d_l(i) = 0;
        end
    end
    for i = 1:x
        if d_x(i) ~= d_x(i)
            d_x(i) = 0;
        end
    end
end
