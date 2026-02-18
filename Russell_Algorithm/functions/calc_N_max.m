function N_max = calc_N_max(r1, r2, theta, T_max, mu)

N_max_search = 10;
N_max = 0;

options = optimset('Display', 'off');

for N = 0:N_max_search
    [alpha, beta, s] = get_lambert_angles(r1, r2, N, theta, true);
    a_m = s/2;

    tof_eqn_fast = lambert_lagrange_eqn(r1, r2, theta, mu, N, true);

    if N > 0
        f = @(a) (6*N*pi + 3*(alpha(a) - beta(a)) - (sin(alpha(a))) - sin(beta(a))) .* ...
            (sin(alpha(a) - beta(a)) + (sin(alpha(a)) - sin(beta(a)))) - 8*(1 - cos(alpha(a) - beta(a)));

        zero_a = fzero(f, [1, 10]*a_m, options);

        tof_min = tof_eqn_fast(zero_a);

        if tof_min < T_max
            N_max = N;
        else
            break
        end
    else

    end
end

end