function f = f_epsilon(x, l, F, rho)
    f_x = F(x);
    m = length(f_x) - 1;
    f = f_x(1) + rho * norm(f_x(2:m+1), 1);
    % f = f_x(1) + dot(f_x(2:m+1), l) + rho * norm(f_x(2:m+1), 1);
end
