function [x, g, j, k] = SQP(x0, F, eps, kiter, choix)
    n = length(x0);
    k = 0;
    rho = 1 / eps^2;
    F_0 = F(x0);
    m = length(F_0) - 1;
    H = eye(n); % differentiel seconde par rapport a x
    c = 0.1;
    [g, j] = Gradient(x0, F, repmat(eps, 1, n));
    [d, l] = pq(H, j', g, -F_0(2:m+1));
    x = x0 + d;
    l0 = l;
    while k < kiter
    %while k < kiter && norm(x-x0) >= eps % use gradient, why ?
    %while k < kiter && norm(g) >= eps
        g0 = g;
        j0 = j;
        [g, j] = Gradient(x, F, repmat(eps, 1, n));
        H0 = H;
        H = Qnewton(x, x0, l, g, g0, j, j0, H0, choix);
        H = H_modification(H);
        F_0 = F(x0);
        x0 = x;
        [d, l] = pq(H, j', g, -F_0(2:m+1));
        [x, dx] = Globalisation(x0, l0, d, F, c, rho);
        d = dx;
        l0 = l;
        k = k + 1;
    end
end

% F = @(x) [x' * x; sum(x)-1];
%
% x0 = [-1; 2; 1; -2; -2];
%
 F = @(x) [(x(1)-1)^2+(x(1)-x(2))^2+(x(2)-x(3))^3+(x(3)-x(4))^4+(x(4)-x(5))^4; x(1)+x(2)^2+x(3)^2-3*sqrt(2)-2; x(2)-x(3)^2+x(4)-2*sqrt(2)+2; x(1)*x(5)-2];
