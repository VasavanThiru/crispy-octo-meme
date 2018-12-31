function [x, g, j, k] = SQP(x_0, F, eps, domain, kiter, choix, to_print)
    g_eps = 1e-5;
    n = length(x_0);
    k = 0;
    rho = 1e9;
    F_0 = F(x_0);
    m = length(F_0) - 1;
    H = eye(n); % differentiel seconde par rapport a x
    c = 0.1;
    [g, j] = Gradient(x_0, F, repmat(g_eps, 1, n));
    [d, l] = pq(H, j', g, -F_0(2:m+1));
    l_0 = l;
    %x = x_0 + d;
    g_0 = g;
    [x, dx, H] = Globalisation(x_0, l_0, d, F, F_0, g_0, H, c, rho);
    d = dx;
    x = Proj(x, domain);
    d = x - x_0;
    l_0 = l;
    g_0 = g + 1 + eps;
    if to_print == 1
        for i = 1:length(x_0)
            fprintf(" %f", x_0(i));
        end
        fprintf("\n");
    end
    % SQP, le grand mystere
    while k < kiter && norm(g, 1) >= eps
    %while k < kiter && norm(g, 1) >= eps && norm(x - x_0, 1) >= eps
    %while k < kiter
        if to_print == 1
            for i = 1:length(x)
                fprintf(" %f", x(i));
            end
            fprintf("\n");
        end
        g_0 = g;
        j_0 = j;
        [g, j] = Gradient(x, F, repmat(g_eps, 1, n));
        H_0 = H;
        H = Qnewton(x, x_0, l, g, g_0, j, j_0, H_0, choix);
        H = H_modification(H, 1e3);
        F_0 = F(x_0);
        x_0 = x;
        [d, l] = pq(H, j', g, -F_0(2:m+1));
        [x, dx, H] = Globalisation(x_0, l_0, d, F, F_0, g_0, H, c, rho);
        d = dx;
        x = Proj(x, domain);
        d = x - x_0;
        l_0 = l;
        k = k + 1;
    end
end
