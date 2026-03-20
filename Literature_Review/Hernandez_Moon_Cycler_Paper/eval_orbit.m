function flybys = eval_orbit(r0, v0, ai, ei, w0, mu, r_all, t0, retrograde)

num_b = length(r_all);
flybys = cell(num_b, 2);

for n = 1:num_b
    r_b = r_all(n);

    [r_b_1, r_b_2, r_sc_1, r_sc_2, tof1, tof2] = ...
        calc_planet_pos_intersection(r0, ai, ei, mu, r_b, w0, t0, retrograde);

    dist1 = norm(r_b_1 - r_sc_1);
    dist2 = norm(r_b_2 - r_sc_2);

    if tof1  == 0
        dist1 = inf;
        tof1 = inf;
    end

    if tof2 == 0
        dist2 = inf;
        tof2 = inf;
    end

    flybys{n, 1} = {r_b_1, r_sc_1, dist1, tof1};
    flybys{n, 2} = {r_b_2, r_sc_2, dist2, tof2};
end
end