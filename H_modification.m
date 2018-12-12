function[G] = H_modification(H,eps)
    v = eig(H);
    v_min = min(v);
    if v_min < 0
        G = H - (v_min + eps) * eye(length(v_min));
    else
        G = H;
    end
end

