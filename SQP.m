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
    [x, dx, H] = Globalisation(x_0, l_0, d, F, g_0, H, c, rho);
    d = dx;
    x = Proj(x, domain);
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
        [x, dx, H] = Globalisation(x_0, l_0, d, F, g_0, H, c, rho);
        x = Proj(x, domain);
        d = dx;
        l_0 = l;
        k = k + 1;
    end
end

% F = @(x) [x' * x; sum(x)-1];
%
% x_0 = [-1; 2; 1; -2; -2];
%
% F = @(x) [(x(1)-1)^2+(x(1)-x(2))^2+(x(2)-x(3))^3+(x(3)-x(4))^4+(x(4)-x(5))^4; x(1)+x(2)^2+x(3)^2-3*sqrt(2)-2; x(2)-x(3)^2+x(4)-2*sqrt(2)+2; x(1)*x(5)-2];
% domain = [[-2, 0]; [1, 3]; [0, 2]; [-2, 0]; [-2, 0]];
% xs = [-1.2366; 2.4616; 1.1911; -0.2143; -1.6165];
%
% Ariane test
%
% steps = [[0.1101, 2647.2]; [0.1532, 2922.4]; [0.2154, 4344.3]];
% m_0 = [2e5; 4e4; 1e4];
% m_u = 1700;
% V = 11527;
% m_3 = @(m) steps(3, 2) * log((m_u+(1+steps(3, 1))*m(3))/(m_u+steps(3, 1)*m(3)));
% m_2 = @(m) steps(2, 2) * log((m_u+(1+steps(3, 1))*m(3)+(1+steps(2, 1))*m(2))/(m_u+(1+steps(3, 1))*m(3)+steps(2, 1)*m(2)));
% m_1 = @(m) steps(1, 2) * log((m_u+(1+steps(3, 1))*m(3)+(1+steps(2, 1))*m(2)+(1+steps(1, 1))*m(1))/(m_u+(1+steps(3, 1))*m(3)+(1+steps(2, 1))*m(2)+steps(1, 1)*m(1)));
% F = @(m) [m_u + dot((1 + steps(1:3, 1)), m); -V + m_1(m) + m_2(m) + m_3(m)];
% domain = [[0, 2e5]; [0, 1e5]; [0, 1e5]];
% ms = [145349; 31215; 7933];
