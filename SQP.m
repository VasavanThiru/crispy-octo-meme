function [x, g, j, k] = SQP(x_0, F, eps, b, kiter, choix)
    n = length(x_0);
    k = 0;
    rho = 1 / eps^2;
    F_0 = F(x_0);
    m = length(F_0) - 1;
    H = eye(n); % differentiel seconde par rapport a x
    c = 0.1;
    [g, j] = Gradient(x_0, F, repmat(eps, 1, n))
    [d, l] = pq(H, j', g, -F_0(2:m+1))
    x = x_0 + d
    l_0 = l;
    if norm(x - x_0, 2) > b
        fprintf("Out of bound\n")
    end
    while k < kiter && norm(d, 1) + norm(l, 1) >= eps
    %while k < kiter && norm(x-x_0) >= eps % use gradient, why ?
    %while k < kiter && norm(g) >= eps
        g_0 = g;
        j_0 = j;
        [g, j] = Gradient(x, F, repmat(eps, 1, n));
        H_0 = H;
        H = Qnewton(x, x_0, l, g, g_0, j, j_0, H_0, choix);
        H = H_modification(H);
        F_0 = F(x_0);
        x_0 = x;
        [d, l] = pq(H, j', g, -F_0(2:m+1));
        [x, dx, H] = Globalisation(x_0, l_0, d, F, g_0, H, c, rho);
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
% xs = [-1.2366; 2.4616; 1.1911; -0.2143; -1.6165];
