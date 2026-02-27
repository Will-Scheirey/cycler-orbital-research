function [N_max, a_all, vd_all, va_all, sol_type] = calc_multirev_lambert(r1_vec, r2_vec, mu, tof)

r1 = norm(r1_vec);
r2 = norm(r2_vec);

theta = prograde_theta(r1_vec, r2_vec);

N_max = calc_N_max(r1, r2, theta, tof, mu);

a_all = nan(N_max, 4);
vd_all = cell(N_max, 4);
va_all = cell(N_max, 4);

sol_type = {{true, true}, {false, true}, {true, false}, {false, false}};

tof_target = tof;

for N = 0:N_max
    
    [a1, vd1, va1] = calc_lambert(theta, true, true);
    [a2, vd2, va2] = calc_lambert(theta, false, true);
    [a3, vd3, va3] = calc_lambert(theta, true, false);
    [a4, vd4, va4] = calc_lambert(theta, false, false);

    a_all (N+1, 1) = a1;
    vd_all{N+1, 1} = vd1;
    va_all{N+1, 1} = va1;

    a_all (N+1, 2) = a2;
    vd_all{N+1, 2} = vd2;
    va_all{N+1, 2} = va2;

    a_all (N+1, 3) = a3;
    vd_all{N+1, 3} = vd3;
    va_all{N+1, 3} = va3;

    a_all (N+1, 4) = a4;
    vd_all{N+1, 4} = vd4;
    va_all{N+1, 4} = va4;
end

    function [a, vd, va] = calc_lambert(theta, fast, long)
        
        if long
            theta = 2*pi - theta;
        end

        [~, ~, s, c] = get_lambert_angles(r1, r2, N, theta, fast);
        a_m = s/2;
        a_min = a_m*(1+1e-3);
        a_max = 1e6*a_m;

        [tof_eqn, alpha, beta] = lambert_lagrange_eqn(r1, r2, theta, mu, N, fast);
        eqn = @(a) tof_eqn(a) - tof_target;
        
        a = NaN;
        if sign(eqn(a_min)) ~= sign(eqn(a_max))
            a = fzero(eqn, [a_min, a_max]);
        end

        [vd, va] = calc_lambert_v(mu, a, c, alpha(a), beta(a), r1_vec, r2_vec);
    end

end