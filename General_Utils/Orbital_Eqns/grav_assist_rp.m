function rp = grav_assist_rp(vA, vD, mu)
    vA = vA(:); vD = vD(:);                % ensure column vectors
    vA_norm = norm(vA);
    vD_norm = norm(vD);

    v_inf = 0.5*(vA_norm + vD_norm); % Take the average in case they dont match perfectly

    c = dot(vA, vD) / max(vA_norm*vD_norm, eps);
    c = max(min(c, 1), -1);
    delta = acos(c);

    s = max(sin(0.5*delta), eps);
    rp = mu / (v_inf^2) * (1/s - 1);
end