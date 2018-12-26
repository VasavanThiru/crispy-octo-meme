function [x, dx, H] = Globalisation(x_0, l, d, F, g_0, H_i, c, rho)
    % x_0: point initial, x_k
    % l: lambda k+1
    % d: 
    % F: fonction cout et contraintes
    % g_0: gradient de la fonction cout au temps k
    % H_i: hessien 
    % c: valeur par default 1/10=0.1
    % rho: valeur tres grande de l'ordre 10^9
    % BLS, the mystery
    %m = length(F(x)) - 1;
    kmax = 64;
    s = 1;
    n = length(x_0);
    m = length(l);
    subindex = @(A, r, c) A(r:c);
    rho_c_0 = rho * norm(subindex(F(x_0), 2, m + 1), 1);
    % f_eps = @(X) [[1; l]' * F(X) + rho * norm(subindex(F(X), 2, m + 1), 1)];
    f_eps = @(X) [subindex(F(X), 1, 1) + rho * norm(subindex(F(X), 2, m + 1), 1)];
    f_0 = f_eps(x_0);
    if f_eps(x_0 + s * d) >= f_0 % La fonction de merite ne decroit pas
        k = 0;
        y = g_0' * d - rho_c_0;
        while y > 0 && k < m % Recherche d'une direction de descente f_eps'(x_0) < 0
            rho = rho + max(l);
            rho_c_0 = rho * norm(subindex(F(x_0), 2, m + 1), 1);
            % f_eps = @(X) [[1; l]' * F(X) + rho * norm(subindex(F(X), 2, m + 1), 1)];
            f_eps = @(X) [subindex(F(X), 1, 1) + rho * norm(subindex(F(X), 2, m + 1), 1)];
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
        if k == kmax
            fprintf("Reinitialisation du hessien\n");
            H_i = eye(n);
        end
    end
    x = x_0 + s * d;
    dx = s * d;
    H = H_i;
end

