function results = generate_guess_flybys(ai, ei, phi0, r0, mu, r_all, t0, retrograde)
    [r0_vec, v0_1, w0_1] = calc_orbit_vecs(ai, ei, phi0, r0, mu, true, retrograde);
    [~,      v0_2, w0_2] = calc_orbit_vecs(ai, ei, phi0, r0, mu, false, retrograde);

    flybys_1 = eval_orbit(r0_vec, v0_1, ai, ei, w0_1, mu, r_all, t0, retrograde);
    flybys_2 = eval_orbit(r0_vec, v0_2, ai, ei, w0_2, mu, r_all, t0, retrograde);

    results = {{flybys_1, r0_vec, v0_1, w0_1}, {flybys_2, r0_vec, v0_2, w0_2}};
end