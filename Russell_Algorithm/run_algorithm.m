
function [N_max, tof, sol_all] = run_algorithm(tau, p, h, s, n, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec)

sol_all = [];

tof = identical_s_tof(tau, p, h, s) * year2sec;
if tof < 0
    N_max = -1;
    return;
end

theta_b1 = n * tof;

r1 = r_b1 * [cos(theta_b1); sin(theta_b1); 0];
v1 = v_b1 * [-sin(theta_b1); cos(theta_b1); 0];

[N_max, a_all, vd_all, va_all] = calc_multirev_lambert(r0, r1, mu, tof);

if N_max == 0
    return
end

sol_all = struct();
sol_idx = 1;

for i = 1 : N_max+1

    fast_struct = struct();

    a_fast = a_all(i, 1);
    vd_fast = vd_all{i, 1};
    va_fast = va_all{i, 1};

    [v_inf_minus_fast, deltas_fast] = calc_sequence(vd_fast, va_fast, h, s, v0, v1);
    if ~isempty(v_inf_minus_fast)
        ra_fast = ra_from_rv(r0, vd_fast, mu);
        [AR_fast, TR_fast, feasible_fast, max_delta_fast] = is_feasible(a_fast, rp_min, mu_b1, r_b2, norm(v_inf_minus_fast{1}), deltas_fast, TR_min, AR_min, ra_fast);
        fast_struct = struct( ...
            'v_inf', norm(v_inf_minus_fast{1}), ...
            'AR', AR_fast, ...
            'TR', TR_fast, ...
            'feasible', feasible_fast, ...
            'max_delta', max_delta_fast, ...
            'i', i-1, ...
            'va', va_fast, ...
            'vd', vd_fast ...
            );
        fast_struct.v_inf_minus = v_inf_minus_fast;
        fast_struct.deltas = deltas_fast;
    end

    slow_struct = struct();

    a_slow = a_all(i, 2);
    vd_slow = vd_all{i, 2};
    va_slow = va_all{i, 2};

    [v_inf_minus_slow, deltas_slow] = calc_sequence(vd_slow, va_slow, h, s, v0, v1);
    if ~isempty(v_inf_minus_slow)
        ra_slow = ra_from_rv(r0, vd_slow, mu);
        [AR_slow, TR_slow, feasible_slow, max_delta_slow] = is_feasible(a_slow, rp_min, mu_b1, r_b2, norm(v_inf_minus_slow{1}), deltas_slow, TR_min, AR_min, ra_slow);
        slow_struct = struct( ...
            'v_inf', norm(v_inf_minus_slow{1}), ...
            'AR', AR_slow, ...
            'TR', TR_slow, ...
            'feasible', feasible_slow, ...
            'max_delta', max_delta_slow, ...
            'i', i-1, ...
            'va', va_slow, ...
            'vd', vd_slow ...
            );
        slow_struct.v_inf_minus = v_inf_minus_slow;
        slow_struct.deltas = deltas_slow;
    end

    sol_all.(sprintf('sol_%d',sol_idx)) = struct( ...
        'fast', fast_struct, ...
        'slow', slow_struct...
        );

    sol_idx = sol_idx + 1;
end

end