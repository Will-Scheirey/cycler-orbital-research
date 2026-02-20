clear; clc; close all
load_params

figure(1)
clf
hold on

p_targ = 2;
h_targ = 5;
s_targ = 1;
i_targ = 0;

fast = false;
long = false;

[N_max, tof, rev_all] = russell_algorithm(tau, p_targ, h_targ, s_targ, mean_motion, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec);

data = struct('p', p_targ, 'h', h_targ, 's', s_targ, 'N_max', N_max, 'tof', tof, 'rev_data', rev_all);

feasible_solutions = get_feasible_solutions({data}, true);

target_solutions = get_solution(p_targ, h_targ, s_targ, i_targ, feasible_solutions, fast, long);

theta_earth = mean_motion * tof;
r1 = r_b1 * [cos(theta_earth); sin(theta_earth); 0];
v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];

plot_orbit(r0, v0, mu, 'DisplayName', 'Earth', 'LineStyle', '--', 'LineWidth', 2); hold on
plot_orbit(r0_2, v0_2, mu, 'DisplayName', 'Mars', 'LineStyle', '--', 'LineWidth', 2);

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

theta_end1 = theta_earth;

if mod(h, 2) ~= 0
    theta_end1 = theta_end1 + pi;
end
r1_end = r_b1 * [cos(theta_end1); sin(theta_end1); 0];
v1_end = v_b1 * [-sin(theta_end1); cos(theta_end1); 0];

for i = 1:length(target_solutions)
    solution = target_solutions{i};
    plot_solution(solution, mu, r_b1, v_b1, theta_earth)
end

lim = r_b2 * 1.1;

theta_end = p_targ * angle_turn;

plot(r1_end(1), r1_end(2), 'rx', 'MarkerSize', 20, 'LineWidth', 3, 'HandleVisibility', 'off')
plot_intersection(r1_end, lim/20, "End")

plot_intersection(r0, lim/20, "0")
plot_intersection(r_intersect, lim/20, "1")
plot_intersection(r1, lim/20, "2")

title(sprintf("Solutions for %d.%d.%d.%d; TOF = %0.2f years", p_targ, h_targ, s_targ, i_targ, tof * sec2year))

xlim([-1,1]*lim)
ylim([-1,1]*lim)
axis square

legend

function plot_intersection(r, offset_x, text_str)
    plot(r(1), r(2), 'k.', 'MarkerSize', 25, 'LineWidth', 2, 'HandleVisibility', 'off')
    text(r(1) + offset_x, r(2), text_str)
end

function plot_solution(solution, mu, r_b1, v_b1, theta_generic)

    k_all = 1:length(solution.v_inf_minus) - 1;
    num_k = length(k_all);

    theta_earth = theta_generic;
    for k_idx = 1:num_k
        k = k_all(k_idx);

        r1 = r_b1 * [cos(theta_earth); sin(theta_earth); 0];
        v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];

        vd = solution.v_inf_minus{k} + v1;
        
        if k == 1
            name = "Trajectory 0; Generic";
        else
            name = sprintf("Trajectory %d; %0.1f-rev", k-1, solution.dt_years(k-1));
        end

        plot_orbit(r1, vd, mu, 'DisplayName', name, 'LineWidth', 1.5);

        if k > 2
            theta_earth = theta_earth + solution.dt_years(k_idx - 1) * 2*pi;
        end
    end

    % plot_orbit(r1_end, v1_end + solution.v_inf_minus{end}, mu, 'DisplayName', 'Generic Return Reinitiated', 'LineWidth', 1.5);

end