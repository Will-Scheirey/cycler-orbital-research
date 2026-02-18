load_params

p_targ = 2;
h_targ = 1;
s_targ = 1;
i_targ = 2;

p = p_targ;
h = h_targ;
s = s_targ;

[N_max, tof, sol_all] = run_algorithm(tau, p, h, s, n, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec);

data = struct('p', p, 'h', h, 's', s, 'N_max', N_max, 'tof', tof, 'solution_data', sol_all);

feasible_solutions = get_feasible_solutions({data}, true);

solution_slow = get_solution(p_targ, h_targ, s_targ, i_targ, '+', feasible_solutions);
if ~isempty(solution_slow)
    plot_orbit(r0, solution_slow.vd, mu, 'Slow'); hold on
    quiver(r0(1), r0(2), solution_slow.vd(1), solution_slow.vd(2), 5e6, 'LineWidth', 2, 'DisplayName', 'V Departure, Slow')
    tof = solution_slow.tof;
    print_solutions(solution_slow)
end

solution_fast = get_solution(p_targ, h_targ, s_targ, i_targ, '-', feasible_solutions);
if ~isempty(solution_fast)
    plot_orbit(r0, solution_fast.vd, mu, 'Fast'); hold on
    quiver(r0(1), r0(2), solution_fast.vd(1), solution_fast.vd(2), 5e6, 'LineWidth', 2, 'DisplayName', 'V Departure, Fast')
    print_solutions(solution_fast)
end

plot_orbit(r0, v0, mu, "Earth")
plot_orbit(r0_2, v0_2, mu, "Mars")

theta_earth = n * tof;
lim = r_b2 * 1.1;

plot(r_b1 * cos(theta_earth), r_b1 * sin(theta_earth), 'k.', 'MarkerSize', 25, 'LineWidth', 2, 'HandleVisibility', 'off')
text(r_b1 * cos(theta_earth) + lim/20, r_b1 * sin(theta_earth), "2")

plot(r0(1), r0(2), 'k.', 'MarkerSize', 25, 'LineWidth', 2, 'HandleVisibility', 'off')
text(r0(1) + lim/20, r0(2), "0")

title(sprintf("Solutions for %d.%d.%d.%d; TOF = %0.2f years", p_targ, h_targ, s_targ, i_targ, tof * sec2year))

xlim([-1,1]*lim)
ylim([-1,1]*lim)
axis square

legend

return
