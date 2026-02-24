%{
function rp = grav_assist_rp(v_A, v_D, mu)
    DV1 = v_A - v_D;
    DV_norm = norm(DV1);
    v_A_norm = norm(v_A);
    v_D_norm = norm(v_D);
    
    dots = dot(v_A, v_D);
    turn_angle = acos(dots ./ (v_A_norm .* v_D_norm));
    
    rp = mu / DV_norm * (1/sin(turn_angle / 2) - 1);
end
%}
function rp = grav_assist_rp(vA, vD, mu)
    % vA, vD: incoming/outgoing v_infty vectors in planet frame
    % mu: planet GM

    vA = vA(:); vD = vD(:);                % ensure column vectors
    vA_norm = norm(vA);
    vD_norm = norm(vD);

    % For a pure (unpowered) flyby, magnitudes should match:
    v_inf = 0.5*(vA_norm + vD_norm);       % robust to tiny mismatch

    % Turn angle between vA and vD (clamp for numerical safety)
    c = dot(vA, vD) / max(vA_norm*vD_norm, eps);
    c = max(min(c, 1), -1);
    delta = acos(c);

    s = max(sin(0.5*delta), eps);          % avoid divide-by-zero
    rp = mu / (v_inf^2) * (1/s - 1);
end