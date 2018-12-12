function [x]=Proj(x_0, domain)
    % x_0: point
    % domain: partie de R^n
    n = length(x_0);
    x = zeros(n, 1);
    for i = 1:n
        if x_0(i) < domain(i, 1)
            x(i) = domain(i, 1);
        elseif x_0(i) > domain(i, 2)
            x(i) = domain(i, 2);
        else
            x(i) = x_0(i);
        end
    end
end
