function [N_max, a_all, vd_all, va_all] = calc_multirev_lambert(r1_vec, r2_vec, mu, tof)

r1 = norm(r1_vec);
r2 = norm(r2_vec);

theta = prograde_theta(r1_vec, r2_vec);

N_max = calc_N_max(r1, r2, theta, tof, mu);
theta = 2*pi - theta;

a_all = nan(N_max, 2);
vd_all = cell(N_max, 2);
va_all = cell(N_max, 2);

for N = 0:N_max
    [~, ~, s, c] = get_lambert_angles(r1, r2, N, theta, true);
    a_m = s/2;

    [tof_eqn_fast, alpha_fast, beta_fast] = lambert_lagrange_eqn(r1, r2, theta, mu, N, true);
    [tof_eqn_slow, alpha_slow, beta_slow] = lambert_lagrange_eqn(r1, r2, theta, mu, N, false);

    tof_target = tof;

    eqn_fast = @(a) tof_eqn_fast(a) - tof_target;
    eqn_slow = @(a) tof_eqn_slow(a) - tof_target;

    a_min = a_m*(1+1e-3);
    a_max = 1e6*a_m;

    a_fast = NaN; a_slow = NaN;

    if sign(eqn_fast(a_min)) ~= sign(eqn_fast(a_max))
        a_fast = fzero(eqn_fast, [a_min, a_max]);
    end

    if sign(eqn_slow(a_min)) ~= sign(eqn_slow(a_max))
        a_slow = fzero(eqn_slow, [a_min, a_max]);
    end

    [vd_fast, va_fast] = calc_lambert_v(mu, a_fast, c, alpha_fast(a_fast), beta_fast(a_fast), r1_vec, r2_vec);
    [vd_slow, va_slow] = calc_lambert_v(mu, a_slow, c, alpha_slow(a_slow), beta_slow(a_slow), r1_vec, r2_vec);

    a_all (N+1, 1) = a_fast;
    vd_all{N+1, 1} = vd_fast;
    va_all{N+1, 1} = va_fast;

    a_all (N+1, 2) = a_slow;
    vd_all{N+1, 2} = vd_slow;
    va_all{N+1, 2} = va_slow;
end

end