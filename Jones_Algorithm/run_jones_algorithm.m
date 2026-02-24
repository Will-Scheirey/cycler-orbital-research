clear; clc;

[mu, ~, r_b3, r_b1, r_b2] = get_circular_planets();

v_b1 = sqrt(mu / r_b1);

[dv1, dv2, a_t, tof_h] = calc_hohmann(mu, r_b1, r_b2);

[~, ~, ~, synodic_s] = calc_synodic3([r_b1, r_b2, r_b3], mu);

S_12  = calc_synodic(mu, r_b1, r_b2);
S_123 = synodic_s;

S_ratio = S_123 / S_12;

t_wait = calc_hohmann_wait(r_b1, r_b2, mu);

t0_star_all = S_12 * (0:ceil(S_ratio)) + t_wait;

t0_range = year2sec(1/20);
dt0      = year2sec(1/20);

figure(1)
clf
for t_star_idx = 1:1
    t0_star = t0_star_all(t_star_idx);

    t0_min = (t0_star - t0_range);
    t0_max = (t0_star + t0_range);

    tof_min = (tof_h - t0_range);
    tof_max = (tof_h + t0_range);

    t0_all   = t0_min  : dt0 : t0_max;
    tof1_all = tof_min : dt0 : tof_max;

    for t0_idx = 1:length(t0_all)

        t0 = t0_all(t0_idx);
        
        for tof_idx = 1:length(tof1_all)

            tof = tof1_all(tof_idx);
            tf = t0 + tof;

            r0 = get_circular_r(mu, r_b1, t0);
            rf = get_circular_r(mu, r_b2, tf);

            [N_max, a_all, vd_all, va_all, sol_type] = calc_multirev_lambert(r0, rf, mu, tof);

            for n = 1:length(vd_all)
                [r, h, e, p, w, nu] = calc_orbit(r0, vd_all{n}, mu);

                theta_intersect = theta_from_r(p, e, w, r_b2);
                plot_orbit_theta(r0, vd_all{n}, mu, [nu+w, theta_intersect], '-k', 'HandleVisibility', 'off'); hold on
            end
            plot(r0(1), r0(2), 'rx', 'MarkerSize', 10, 'LineWidth', 3, 'HandleVisibility', 'off')
            plot(rf(1), rf(2), 'bx', 'MarkerSize', 10, 'LineWidth', 3, 'HandleVisibility', 'off')

        end
    end
end


theta_plot = linspace(0, 2*pi, 1000);


plot(r_b1 .* cos(theta_plot), r_b1 .* sin(theta_plot), '-b', 'DisplayName', 'Earth', 'LineWidth', 2); hold on
plot(r_b2 .* cos(theta_plot), r_b2 .* sin(theta_plot), '-r', 'DisplayName', 'Mars', 'LineWidth', 2);
plot(r_b3 .* cos(theta_plot), r_b3 .* sin(theta_plot), '-', 'DisplayName', 'Venus', 'LineWidth', 2, 'Color', [0.9, 0.8, 0]);

axis equal
legend

function [nu1, nu2, nu3] = get_true_anomalies(t, r_b1, r_b2, r_b3, mu)
    n1 = mean_motion(r_b1, mu);
    n2 = mean_motion(r_b2, mu);
    n3 = mean_motion(r_b3, mu);

    nu1 = mod(n1*t, 2*pi);
    nu2 = mod(n2*t, 2*pi);
    nu3 = mod(n3*t, 2*pi);
end