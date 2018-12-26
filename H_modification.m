function [G] = H_modification(H, tau)
    % H: hessien
    % tau: grande valeur
    v = eig(H);
    v_min = min(v);
    if v_min < 0
        G = H + (tau - v_min) * eye(length(H));
    else
        G = H;
    end
end

