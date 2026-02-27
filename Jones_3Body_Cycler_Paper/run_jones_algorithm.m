% Based on: https://ntrs.nasa.gov/citations/20190028464

clear; clc; close all

[mu, ~, r_b3, r_b1, r_b2] = get_circular_planets();

v_b1 = sqrt(mu / r_b1);

[dv1, dv2, a_t, tof_h] = calc_hohmann(mu, r_b1, r_b2);

[~, ~, ~, synodic_s] = calc_synodic3([r_b1, r_b2, r_b3], mu);

S_12  = calc_synodic(mu, r_b1, r_b2);
S_123 = synodic_s;

tof_tot = S_123 * 1;
tof_leg_max = year2sec(2);
tof_leg_min = year2sec(0.2);

S_ratio = S_123 / S_12;

t_wait = calc_hohmann_wait(r_b1, r_b2, mu);

t0_all = S_12 * (0:ceil(S_ratio)) + t_wait;

t0_range = year2sec(1/10);
dt0      = year2sec(1/20);

r_all = [r_b1, r_b2, r_b3];

sequence = [1, 2, 3, 2, 1];
% sequence = [1, 2, 1];
num_flyby = length(sequence);

seq_all = lambert_seq(sequence, r_all, t_wait, tof_h, t0_range, dt0, tof_tot, mu);

this_sequence = seq_all(1, :);
this_leg = this_sequence;

t0 = t_wait;

ri = get_circular_r(mu, r_all(sequence(1)), t0);

figure(1)
clf
hold on

theta_plot = linspace(0, 2*pi, 1000);

plot(r_b1 .* cos(theta_plot), r_b1 .* sin(theta_plot), '-b', 'DisplayName', 'Earth', 'LineWidth', 2); hold on
plot(r_b2 .* cos(theta_plot), r_b2 .* sin(theta_plot), '-r', 'DisplayName', 'Mars', 'LineWidth', 2);
plot(r_b3 .* cos(theta_plot), r_b3 .* sin(theta_plot), '-', 'DisplayName', 'Venus', 'LineWidth', 2, 'Color', [0.9, 0.8, 0]);


r_end = get_circular_r(mu, r_b1, t0 + tof_tot);
plot(r_end(1), r_end(2), 'rx', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'End (Theoretical)')

idx_stop = num_flyby-1;

idx = 1;
while idx <= idx_stop
    vd = this_leg{2};
    leg_tof = this_leg{4};
    rf = this_leg{5};
    ri2 = this_leg{6};

    tf = t0 + leg_tof;

    plot(rf(1), rf(2), 'k.', 'MarkerSize', 20, 'HandleVisibility', 'off');

    plot(ri2(1), ri2(2), 'k.', 'MarkerSize', 20, 'HandleVisibility', 'off');

    stop = false;

    if height(vd) > 1
        vd = vd(1, :);
    end

    for i = 1:numel(vd)
        if ~any(isnan(vd{i}))
            stop = true;
            break
        end
    end

    vd = vd{i};

    if stop
        plot_orbit_from_state(ri, vd, mu, leg_tof, 'HandleVisibility', 'off');
        quiver(ri(1), ri(2), vd(1), vd(2), 3e6, 'LineWidth', 3, 'DisplayName', sprintf('Depart %d', idx))
    end

    if idx == idx_stop
        plot(rf(1), rf(2), 'g.', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'End (Actual)')
        break
    end

    ri = rf;
    this_leg = this_leg{7}(1, :);
    idx = idx + 1;
    t0 = tf;
end

axis equal
legend

function data_out = lambert_seq(sequence, r_all, t0, tof_0, t_range, dt, tof_tot, mu)
    num_seq = length(sequence);

    if num_seq == 1
        data_out = {};
        return
    end

    tof_leg = (tof_0 - t_range) : dt : (tof_0 + t_range);
    num_tof = length(tof_leg);

    if num_seq == 2
        tof_leg = tof_tot;
        num_tof = 1;
    end

    data_out = cell(num_tof, 7);

    b1_idx = sequence(1);
    b2_idx = sequence(2);

    for idx=1:num_tof
        tof = tof_leg(idx);

        t1 = t0;
        t2 = t0 + tof;

        t_remain = tof_tot - tof;

        if t_remain < 0
            continue;
        end

        r1 = get_circular_r(mu, r_all(b1_idx), t1);
        r2 = get_circular_r(mu, r_all(b2_idx), t2);

        [N_max, ~, vd_all, va_all, ~] = calc_multirev_lambert(r1, r2, mu, tof);

        data_out{idx, 1} = 0:N_max;
        data_out{idx, 2} = vd_all;
        data_out{idx, 3} = va_all;
        data_out{idx, 4} = tof;
        data_out{idx, 5} = r2;
        data_out{idx, 6} = r1;
        data_out{idx, 7} = lambert_seq(sequence(2:end), r_all, t2, tof_0, t_range, dt, t_remain, mu);
    end
end