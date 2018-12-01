function [g, j] = Gradient(x, F, delta)
    n = length(x);
    F_x = F(x);
    m = length(F_x) - 1;
    f_x = F_x(1);
    c_x = F_x(2:m+1);
    g = zeros(n, 1);
    j = zeros(n, m);
    for i = 1:n
        h = x(i);
        x(i) = x(i) + delta(i);
        % [f_xi, c_xi] = F(x);
        F_xi = F(x);
        f_xi = F_xi(1);
        c_xi = F_xi(2:m+1);
        g(i) = (f_xi - f_x) / delta(i);
        j(i, 1:m) = (c_xi - c_x) / delta(i);
        x(i) = h;
    end
end

% F = @(x) [0.5 * (A * x)' * x + d * x; [b * x; c * x]];
