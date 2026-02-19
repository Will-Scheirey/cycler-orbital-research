function [N_max, tof, sol_all] = run_algorithm(tau, p, h, s, mean_motion, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec)

sol_all = [];

tof = identical_s_tof(tau, p, h, s) * year2sec;
if tof < 0
    N_max = -1;
    return;
end

theta_b1 = mean_motion * tof;

r1 = r_b1 * [cos(theta_b1); sin(theta_b1); 0];
v1 = v_b1 * [-sin(theta_b1); cos(theta_b1); 0];

[N_max, a_all, vd_all, va_all, sol_type_all] = calc_multirev_lambert(r0, r1, mu, tof);

if N_max == 0
    return
end

sol_all = struct();
sol_idx = 1;

for i = 1 : N_max+1

    struct_save = struct();
    for z = 1:4
        a = a_all(i, z);
        vd = vd_all{i, z};
        va = va_all{i, z};
    
        [v_inf_minus, deltas, h_all] = calc_sequence(vd, va, h, s, v0, v1);
        if isempty(v_inf_minus)
            continue
        end

        ra = ra_from_rv(r0, vd, mu);
        [AR, TR, feasible, max_delta] = is_feasible(a, rp_min, mu_b1, r_b2, norm(v_inf_minus{1}), deltas, TR_min, AR_min, ra);
        sol_type = sol_type_all{z};

        struct_temp = struct( ...
            'v_inf', norm(v_inf_minus{1}), ...
            'AR', AR, ...
            'TR', TR, ...
            'feasible', feasible, ...
            'max_delta', max_delta, ...
            'i', i-1, ...
            'va', va, ...
            'vd', vd, ...
            'fast', sol_type{1}, ...
            'long', sol_type{2}, ...
            'h_i', h_all ...
            );
        struct_temp.v_inf_minus = v_inf_minus;
        struct_temp.deltas = deltas;

        struct_save.(sprintf('sol_%d', z)) = struct_temp;
    end

    sol_all.(sprintf('rev_%d',sol_idx)) = struct_save;

    sol_idx = sol_idx + 1;
end
end