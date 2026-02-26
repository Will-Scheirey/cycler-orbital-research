clear; clc;

[mu, ~, r_b3, r_b1, r_b2] = get_circular_planets();

v_b1 = sqrt(mu / r_b1);

[dv1, dv2, a_t, tof_h] = calc_hohmann(mu, r_b1, r_b2);

[~, ~, ~, synodic_s] = calc_synodic3([r_b1, r_b2, r_b3], mu);

S_12  = calc_synodic(mu, r_b1, r_b2);
S_123 = synodic_s;

tof_tot = S_123;
tof_leg_max = year2sec(2);
tof_leg_min = year2sec(0.2);

S_ratio = S_123 / S_12;

t_wait = calc_hohmann_wait(r_b1, r_b2, mu);

t0_all = S_12 * (0:ceil(S_ratio)) + t_wait;

t0_range = year2sec(1/20);
dt0      = year2sec(1/20);

r_all = [r_b1, r_b2, r_b3];

sequence = [1, 2, 1, 3, 1, 1];
% sequence = [1, 2, 3,];
num_flyby = length(sequence);

seq_all = lambert_seq(sequence, r_all, t_wait, tof_h, t0_range, dt0, tof_tot, mu);

this_sequence = seq_all(1, :);
this_leg = this_sequence;

t0 = t_wait;

ri = get_circular_r(mu, r_all(sequence(1)), t0);

figure(1)
clf
plot(ri(1), ri(2), 'b.', 'MarkerSize', 40); hold on

idx_stop = num_flyby-1;

idx = 1;
while idx <= idx_stop
    leg_tof = this_leg{3};

    tf = t0 + leg_tof;

    rf = this_leg{4};
    plot(rf(1), rf(2), 'k.', 'MarkerSize', 20);

    ri2 = this_leg{5};
    plot(ri2(1), ri2(2), 'k.', 'MarkerSize', 20);

    stop = false;
    vd = this_leg{1};

    for i = 1:numel(vd)
        if ~any(isnan(vd{i}))
            stop = true;
            break
        end
    end

    vd = vd{i};

    if stop
        [r, h, e, p, w, nu] = calc_orbit(ri, vd, mu);

        nu_f = theta_from_tof(ri, vd, mu, leg_tof);

        theta_start = mod(nu   + w, 2*pi);
        theta_end   = mod(nu_f + w, 2*pi);

        plot_orbit_theta(ri, vd, mu, [theta_start, theta_end], '-k');
    end

    if idx == idx_stop
        break
    end

    ri = rf;
    this_leg = this_leg{6}(1, :);
    idx = idx + 1;
    t0 = tf;
end

theta_plot = linspace(0, 2*pi, 1000);

plot(r_b1 .* cos(theta_plot), r_b1 .* sin(theta_plot), '-b', 'DisplayName', 'Earth', 'LineWidth', 2); hold on
plot(r_b2 .* cos(theta_plot), r_b2 .* sin(theta_plot), '-r', 'DisplayName', 'Mars', 'LineWidth', 2);
plot(r_b3 .* cos(theta_plot), r_b3 .* sin(theta_plot), '-', 'DisplayName', 'Venus', 'LineWidth', 2, 'Color', [0.9, 0.8, 0]);

axis equal
legend
%{

tof_seed = tof_h - t0_range : dt0 : tof_h + t0_range;
tof_leg = tof_leg_min : dt0 : tof_leg_max;

tof_leg1 = tof_leg;
tof_leg5 = tof_leg;

tof_seq = cell(1);
tof_seq_idx = 1;

t0 = t_wait;

for seed_idx = 1:length(tof_seed)

    r0 = get_circular_r(mu, r_b1, )
    rf = get_circular_r(mu, r_b2, tof_seed(seed_idx));

    [~, ~, vd_all, va_all, ~] = calc_multirev_lambert(r0, rf, mu, tof);

    for leg1_idx = 1:length(tof_leg1)
        for leg5_idx = 1:length(tof_leg5)
            tof_remain = tof_tot - tof_seed(seed_idx) - tof_leg1(leg1_idx) - tof_leg5(leg5_idx);
            if tof_remain < year2sec(2)
                continue;
            end

            tstep = tof_remain / (3 * 3);
            tof_leg2 = tof_leg_min : tstep : tof_remain/2;

            for leg2_idx = 1:length(tof_leg2)

                tof_remain_2 = tof_remain - tof_leg2(leg2_idx);
                tstep_remain = tof_remain_2 / (3 * 2);
                leg3_tof = tof_leg_min : tstep_remain : tof_remain_2/2;

                for leg3_idx = 1:length(leg3_tof)
                
                    leg4_tof = tof_remain_2 - leg3_tof(leg3_idx);

                    seq = [tof_seed(seed_idx), tof_leg1(leg1_idx), tof_leg2(leg2_idx), leg3_tof(leg3_idx), leg4_tof, tof_leg5(leg5_idx)];

                    if any(seq < 0)
                        continue
                    end

                    tof_seq{tof_seq_idx} = seq;
                                        
                    tof_seq_idx = tof_seq_idx + 1;

                end
            end
        end
    end
end

% [leg1, leg1_t0, leg1_tf] = calc_leg(r_b1, r_b2, t0_all, tof_h, t0_range, dt0, mu);


return

%}

function data_out = lambert_seq(sequence, r_all, t0, tof_0, t_range, dt, tof_tot, mu)
    num_seq = length(sequence);

    if num_seq == 1
        data_out = {};
        return
    end

    tof_leg = (tof_0 - t_range) : dt : (tof_0 + t_range);
    num_tof = length(tof_leg);

    data_out = cell(num_tof, 5);

    for idx=1:num_tof
        tof = tof_leg(idx);

        t1 = t0;
        t2 = t1 + tof;
        t_remain = tof_tot - tof;

        if t_remain < 0
            continue;
        end

        r1 = get_circular_r(mu, r_all(sequence(1)), t1);
        r2 = get_circular_r(mu, r_all(sequence(2)), t2);

        [~, ~, vd_all, va_all, ~] = calc_multirev_lambert(r1, r2, mu, tof);

        data_out{idx, 1} = vd_all;
        data_out{idx, 2} = va_all;
        data_out{idx, 3} = tof;
        data_out{idx, 4} = r2;
        data_out{idx, 5} = r2;
        data_out{idx, 6} = lambert_seq(sequence(2:end), r_all, t2, tof_0, t_range, dt, t_remain, mu);
    end
end

function [legs_all, t0_all, tf_all] = calc_leg(r1, r2, t_initial, tof_0, t_range, dt, mu)
legs_all = cell(1, 1);
leg_idx = 1;

t0_all = [];
tf_all = [];
t_idx = 1;

for t_i_idx = 1:length(t_initial)
    ti = t_initial(t_i_idx);

    t0_min = (ti - t_range);
    t0_max = (ti + t_range);

    tof_min = (tof_0 - t_range);
    tof_max = (tof_0 + t_range);

    t0_all   = t0_min  : dt : t0_max;
    tof1_all = tof_min : dt : tof_max;

    for t0_idx = 1:length(t0_all)

        t0 = t0_all(t0_idx);
        
        for tof_idx = 1:length(tof1_all)

            tof = tof1_all(tof_idx);
            tf = t0 + tof;

            t0_all(t_idx) = t0;
            tf_all(t_idx) = tf;
            t_idx = t_idx + 1;

            r0 = get_circular_r(mu, r1, t0);
            rf = get_circular_r(mu, r2, tf);

            [~, ~, vd_all, ~, ~] = calc_multirev_lambert(r0, rf, mu, tof);

            for n = 1:length(vd_all)
                legs_all{leg_idx} = vd_all{n};
                leg_idx = leg_idx + 1;
            end
        end
    end
end
end