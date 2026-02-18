function theta = prograde_theta(r1_vec, r2_vec)
    c = cross(r1_vec, r2_vec);
    s = norm(c);
    d = dot(r1_vec, r2_vec);

    theta = atan2(s, d);

    if c(3) < 0
        theta = 2*pi - theta;
    end
end