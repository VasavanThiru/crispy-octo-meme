function[H] = Qnewton(x, x0, l, g, g0, j, j0, H0, choix)
    d0 = x - x0;
    y0 = g + j * l - (g0 + j0 * l);
    if choix == 1
        if y0' * d0 > 0
            H = H0 + y0 * y0' / (y0'*d0) - H0 * d0 * d0' * H0 / (d0' * H0 * d0);
        else
            H = H0;
        end
    else
        if d0' * (y0 - H0 * d0) ~= 0
            H = H0 + (y0 - H0*d0) * (y0 - H0*d0)' / (d0'*(y0 - H0*d0));
        else
            H = H0;
        end
    end
end

