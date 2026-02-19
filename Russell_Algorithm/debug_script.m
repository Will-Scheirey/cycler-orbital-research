clear; clc; close all
load_params

p_targ = 2;
h_targ = 5;
s_targ = 1;
i_targ = 0;

fast = false;
long = false;

[N_max, tof, rev_all] = run_algorithm(tau, p_targ, h_targ, s_targ, mean_motion, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec);

data = struct('p', p_targ, 'h', h_targ, 's', s_targ, 'N_max', N_max, 'tof', tof, 'rev_data', rev_all);

feasible_solutions = get_feasible_solutions({data}, true);

target_solutions = get_solution(p_targ, h_targ, s_targ, i_targ, feasible_solutions, fast, long);

theta_earth = mean_motion * tof;
r1 = r_b1 * [cos(theta_earth); sin(theta_earth); 0];
v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];

the_plot = plot_orbit(r0, v0, mu, "Earth"); hold on
the_plot.LineWidth = 3;

the_plot = plot_orbit(r0_2, v0_2, mu, "Mars");
the_plot.LineWidth = 3;

vd = target_solutions{1}.vd;

[r, h, e, p_generic, w, nu] = calc_orbit(r0, vd, mu);

theta_intersect = theta_from_r(p_generic, e, w, r_b2);
r2 = r_b2;

if imag(theta_intersect) ~= 0
    ra = ra_from_rv(r0, vd, mu);
    theta_intersect = real(theta_from_r(p_generic, e, w, ra));
    r2 = ra;
end

r_intersect = r2 * [cos(theta_intersect); sin(theta_intersect); 0];

theta_end1 = theta_earth + pi;
r1_end = r_b1 * [cos(theta_end1); sin(theta_end1); 0];
v1_end = v_b1 * [-sin(theta_end1); cos(theta_end1); 0];

for i = 1:length(target_solutions)
    solution = target_solutions{i};
    plot_solution(solution, mu, r1, v1, vd, r1_end, v1_end)
end

lim = r_b2 * 1.1;

theta_end = p_targ * angle_turn;
r_end = r_b1 * [cos(theta_end); sin(theta_end); 0];

plot(r1_end(1), r1_end(2), 'rx', 'MarkerSize', 20, 'LineWidth', 3)

plot_intersection(r0, lim/20, "0")
plot_intersection(r_intersect, lim/20, "1")
plot_intersection(r1, lim/20, "2")
plot_intersection(r_end, lim/20, "End")

title(sprintf("Solutions for %d.%d.%d.%d; TOF = %0.2f years", p_targ, h_targ, s_targ, i_targ, tof * sec2year))

xlim([-1,1]*lim)
ylim([-1,1]*lim)
axis square

legend

function plot_intersection(r, offset_x, text_str)
    plot(r(1), r(2), 'k.', 'MarkerSize', 25, 'LineWidth', 2, 'HandleVisibility', 'off')
    text(r(1) + offset_x, r(2), text_str)
end

function plot_solution(solution, mu, r1_generic, v1_generic, vd_generic, r1_end, v1_end)

    r1_mag = norm(r1_generic);
    r0 = [r1_mag; 0; 0];

    v1_mag = norm(v1_generic);
    v0 = [0; v1_mag; 0];

    v_inf_generic = vd_generic - v0;

    plot_orbit(r0, vd_generic, mu, sprintf("Orbit 1"));

    for k = 2:length(solution.v_inf_minus) - 1
        v1 = v1_generic;
        r1 = r1_generic;

        vd = solution.v_inf_minus{k} + v1;

        plot_orbit(r1, vd, mu, sprintf("Orbit %d", k));
    end

    plot_orbit(r1_end, v1_end + solution.v_inf_minus{end}, mu, "Orbit New");

end