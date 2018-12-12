function[H] = Qnewton(x, x0, l, g, g0, j, j0, H0, choix)
    d0 = x - x0;
    y0 = g + j * l - (g0 + j0 * l);
    if choix == 1 % BFGS
        if y0' * d0 > 0
            H = H0 + y0 * y0' / (y0'*d0) - H0 * d0 * d0' * H0 / (d0' * H0 * d0);
        else
            H = H0;
        end
    else % SR1
        %if d0' * (y0 - H0 * d0) ~= 0
        if d0' * (y0 - H0 * d0) > 1e-14 || d0' * (y0 - H0 * d0) < -1e-14
            H = H0 + (y0 - H0*d0) * (y0 - H0*d0)' / (d0'*(y0 - H0*d0));
        else
            H = H0;
        end
    end
end

