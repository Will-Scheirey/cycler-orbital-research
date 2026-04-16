function results = generate_guess_flybys(ai, ei, phi0, r0, mu, r_all, t0, retrograde)
    [r0_vec, v0, w0] = calc_orbit_vecs(ai, ei, phi0, r0, mu, true, retrograde);
    [~,      v0_other, w0_other] = calc_orbit_vecs(ai, ei, phi0, r0, mu, false, retrograde);

    flybys = eval_orbit(r0_vec, v0, ai, ei, w0, mu, r_all, t0, retrograde);
    flybys_other = eval_orbit(r0_vec, v0_other, ai, ei, w0_other, mu, r_all, t0, retrograde);

    results = {{flybys, r0_vec, v0, w0}, {flybys_other, r0_vec, v0_other, w0_other}};
end