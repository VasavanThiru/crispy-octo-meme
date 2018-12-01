function [x, dx] = Globalisation(x0, l0, d, F, c, rho)
    % BLS, the mystery
    %m = length(F(x)) - 1;
    kmax = 1000;
    s = 1;
    n = length(x0);
    m = length(l0);
    subindex = @(A, r, c) A(r:c);
    %f_eps = @(X) F(X)(1) + l0' * F(X)(2:m+1) + rho * norm(F(X)(2:m+1), 1);
    f_eps = @(X) [[1; l0]' * F(X) + rho * norm(subindex(F(X), 2, m + 1), 1)];
    f_0 = f_eps(x0);
    if f_eps(x0 + s * d) >= f_0
        k = 0;
        [g, ~] = Gradient(x0, f_eps, repmat(0.001, 1, n)); % not nice
        y = g' * d;
        while f_eps(x0 + s * d) >= f_0 + c * s * y && k < kmax
            s = s / 2;
            k = k + 1;
        end
    end
    x = x0 + s * d;
    dx = s * d;
end

