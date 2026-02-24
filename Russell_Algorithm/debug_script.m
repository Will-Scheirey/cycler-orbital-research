clear; clc; close all
load_params

%% Run Algorithm

p_targ = 2;
h_targ = 1;
s_targ = 1;
i_targ = 2;

fast = false;
long = false;

[N_max, tof, rev_all] = russell_algorithm(tau, p_targ, h_targ, s_targ, mean_motion, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec);

data = struct('p', p_targ, 'h', h_targ, 's', s_targ, 'N_max', N_max, 'tof', tof, 'rev_data', rev_all);

feasible_solutions = get_feasible_solutions({data}, true);

target_solutions = get_solution(p_targ, h_targ, s_targ, i_targ, feasible_solutions, fast, long);

theta_b1_turn = mean_motion * tof;

%% CALCULATE INTERSECTIONS
vd = target_solutions{1}.vd;

[theta_intersect, r_intersect] = get_intersection(r0, vd, mu, r_b2);

r_intersect_vec = r_intersect * [cos(theta_intersect); sin(theta_intersect); 0];

theta_end = p_targ * angle_turn;

r_end = r_b1 * [cos(theta_end); sin(theta_end); 0];

%% PLOT

figure(1)
clf
hold on

plot_orbit(r0, v0, mu, 'DisplayName', 'Earth', 'LineStyle', '--', 'LineWidth', 2); hold on
plot_orbit(r0_2, v0_2, mu, 'DisplayName', 'Mars', 'LineStyle', '--', 'LineWidth', 2);

for i = 1:length(target_solutions)
    solution = target_solutions{i};
    plot_solution(solution, mu, vd, r_b1, v_b1, theta_b1_turn)
end

lim = r_b2 * 1.1;

plot_intersection(r0, lim/20, "0")
plot_intersection(r_intersect_vec, lim/20, "1")
plot_intersection(r_end, lim/20, "end")

title(sprintf("Solutions for %d.%d.%d.%d; TOF = %0.2f years", p_targ, h_targ, s_targ, i_targ, tof * sec2year))

xlim([-1,1]*lim)
ylim([-1,1]*lim)
axis square

legend

function plot_intersection(r, offset_x, text_str)
    plot(r(1), r(2), 'k.', 'MarkerSize', 25, 'LineWidth', 2, 'HandleVisibility', 'off')
    text(r(1) + offset_x, r(2), text_str)
end

function plot_solution(solution, mu, vd_generic, r_b1, v_b1, theta_generic)

    num_s = solution.s;
    
    theta_earth = theta_generic;

    total_lines = 0;
    for s = 1:num_s
        total_lines = total_lines + length(solution.v_inf_minus{s});
    end

    colors = copper(total_lines+1);
    line_idx = 2;

    color = colors(1, :);

    r0 = [r_b1; 0; 0];

    plot_orbit(r0, vd_generic, mu, 'DisplayName', 'Initial Generic', 'LineWidth', 1.5, 'Color', color);
    quiver(0, 0, r0(1), r0(2), 0, 'LineWidth', 1.5, 'Color', color, 'HandleVisibility', 'off');

    for s = 1:num_s

        v_inf_minus = solution.v_inf_minus{s};
        dt_years    = solution.dt_years{s};

        k_all = 1:length(v_inf_minus);
        num_k = length(k_all);

        if s > 1
            k_min = 1;
        else
            k_min = 2;
        end

        if s < num_s
            k_max = num_k;
        else
            k_max = num_k-1;
        end
        
        for k_idx = k_min:k_max
            k = k_all(k_idx);
        
            r1 = r_b1 * [cos(theta_earth); sin(theta_earth); 0];
            v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];
        
            vd = v_inf_minus{k} + v1;
            
            if k_idx == 1 || k == num_k
                name = sprintf("Seq %d; Trajectory 0; Generic", s);
            elseif k_idx < num_k 
                name = sprintf("Seq %d; Trajectory %d; %0.1f-rev", s, k_idx, dt_years(k_idx - 1));
            end

            theta_rev = 2*pi;
            if k_idx > 1 && k_idx < num_k
                theta_rev = dt_years(k_idx - 1) * 2*pi;
            end
        
            color = colors(line_idx, :);

            plot_orbit_theta(r1, vd, mu, [theta_earth, theta_earth+theta_rev], 'DisplayName', name, 'LineWidth', 1.5, 'Color', color);
            quiver(0, 0, r1(1), r1(2), 0, 'LineWidth', 1.5, 'Color', color, 'HandleVisibility', 'off');
        
            if k_idx > 1 && k_idx < num_k
                theta_earth = theta_earth + dt_years(k_idx - 1) * 2*pi;
            else
                theta_earth = theta_earth + theta_generic;
            end
            line_idx = line_idx + 1;
        end
        % theta_earth = theta_generic * s;
        % r1 = r_b1 * [cos(theta_earth); sin(theta_earth); 0];
        % v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];
        % plot_orbit(r1, v1 + v_inf_minus{end}, mu, 'DisplayName', 'Generic Return Reinitiated', 'LineWidth', 1.5);
    end
end

function [theta, r2] = get_intersection(r0, vd, mu, r_b2)
[r, h, e, p_generic, w, nu] = calc_orbit(r0, vd, mu);

theta = theta_from_r(p_generic, e, w, r_b2);
r2 = r_b2;

if imag(theta) ~= 0
    ra = ra_from_rv(r0, vd, mu);
    theta = real(theta_from_r(p_generic, e, w, ra));
    r2 = ra;
end
end