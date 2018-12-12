function [x, dx, H] = Globalisation(x_0, l, d, F, g_0, H_i, c, rho)
    % BLS, the mystery
    %m = length(F(x)) - 1;
    kmax = 1000;
    s = 1;
    n = length(x_0);
    m = length(l);
    subindex = @(A, r, c) A(r:c);
    %f_eps = @(X) F(X)(1) + l' * F(X)(2:m+1) + rho * norm(F(X)(2:m+1), 1);
    rho_c_0 = rho * norm(subindex(F(x_0), 2, m + 1), 1);
    f_eps = @(X) [[1; l]' * F(X) + rho * norm(subindex(F(X), 2, m + 1), 1)];
    f_0 = f_eps(x_0);
    if f_eps(x_0 + s * d) >= f_0 % La fonction de merite ne decroit pas
        k = 0;
        %[g, ~] = Gradient(x_0, f_eps, repmat(0.001, 1, n)); % not nice
        y = g_0' * d - rho_c_0;
        while y > 0 && k < m % Recherche d'une direction de descente
            rho = rho + max(l);
            rho_c_0 = rho * norm(subindex(F(x_0), 2, m + 1), 1);
            f_eps = @(X) [[1; l]' * F(X) + rho * norm(subindex(F(X), 2, m + 1), 1)];
            %[g, ~] = Gradient(x_0, f_eps, repmat(0.001, 1, n)); % not nice
            y = g_0' * d - rho_c_0;
            k = k + 1;
        end
        if k == m % Si aucun direction de descente n'est trouvee
            fprintf("Reinitialisation du hessien\n");
            H_i = eye(n);
        end
        k = 0;
        while f_eps(x_0 + s * d) >= f_0 + c * s * y && k < kmax
            s = s / 2;
            k = k + 1;
        end
    end
    x = x_0 + s * d;
    dx = s * d;
    H = H_i;
end

