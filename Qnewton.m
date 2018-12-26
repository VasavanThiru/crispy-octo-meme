function[H] = Qnewton(x, x_0, l, g, g_0, j, j_0, H_0, choix)
    % x: position a l'etape k + 1
    % x_0: position a l'etape k
    % l: lambda a l'etape k
    % g: gradient de f a l'etape k + 1
    % g_0: gradient de f a l'etape k
    % j: gradient de g a l'etape k + 1
    % j_0: gradient de g a l'etape k
    % H_0: hessien a l'etape k
    % choix: 1 pour BFGS, 2 pour SR1
    d_0 = x - x_0;
    y0 = g + j * l - (g_0 + j_0 * l);
    H = H_0;
    if choix == 1 % BFGS
        if y0' * d_0 > 0
            H = H_0 + y0 * y0' / (y0'*d_0) - H_0 * d_0 * d_0' * H_0 / (d_0' * H_0 * d_0);
        end
    elseif choix == 2 % SR1
        %if d_0' * (y0 - H_0 * d_0) ~= 0
        if d_0' * (y0 - H_0 * d_0) > 1e-14 || d_0' * (y0 - H_0 * d_0) < -1e-14
            H = H_0 + (y0 - H_0*d_0) * (y0 - H_0*d_0)' / (d_0'*(y0 - H_0*d_0));
        end
    else
        printf("Error: invalid choix")
    end
end

