clear; clc; close all;

% --- Setup ---

[mu_sun, mu_earth, mu_mars,... 
    r_earth_sun, r_mars_sun,... 
    T_earth, T_mars, S, angle_turn] = get_earth_mars_params();

% --- Lambert Solutions --- 

num_ratios = 100;
num_ns = 100;
max_n = 10;
ns = linspace(1, max_n, num_ns);

num_paths = zeros(num_ns, num_ratios);

zlim = 1200;
numz = 1000;

for n_idx = 1:num_ns
n  = ns(n_idx);
    leg1_t = linspace(T_mars / 5, n*S, num_ratios);
for idx = 1:num_ratios
    fprintf("n=%d, idx=%d\n", n, idx)
    clf
    tof_tot = n*S;
    ratio = leg1_t(idx);
    
    t_leg1 = ratio;
    
    theta_mars1 = t_leg1 / T_mars * 2*pi;
    
    pos_mars1 = r_mars_sun * [cos(theta_mars1); sin(theta_mars1); 0];
    
    R11 = [r_earth_sun; 0; 0];
    R21 = pos_mars1;
    tof1 = t_leg1;
    
    [v_D1, v_A1] = lambert(R11, R21, tof1, mu_sun, zlim, numz);
    v_D1 = cell2mat(v_D1); v_A1 = cell2mat(v_A1);
    
    if isempty(v_D1)
    continue
    end
    
    R12 = R21;
    R22 = r_earth_sun * [cos(2*pi*tof_tot/T_earth); sin(2*pi*tof_tot/T_earth); 0];
    tof2 = tof_tot - tof1;
    
    [v_D2, v_A2] = lambert(R12, R22, tof2, mu_sun, zlim, numz);
    v_D2 = cell2mat(v_D2); v_A2 = cell2mat(v_A2);
    
    if isempty(v_D2)
    continue
    end
    
    %{
    theta_mars2 = (tof_tot + t_leg1) / T_mars * 2*pi;
    
    R13 = R22;
    R23 = r_mars_sun * [cos(theta_mars2); sin(theta_mars2); 0];
    tof3 = t_leg1;
    
    [v_D3, v_A3] = lambert(R13, R23, tof3, mu_sun, zlim, numz);
    v_D3 = cell2mat(v_D3); v_A3 = cell2mat(v_A3);
    
    if isempty(v_D3)
    continue
    end
    %}
    % --- Gravity Assist ---
    
    num_paths(n_idx, idx) = calc_num_paths(v_A1, v_D2, mu_mars, 3396);
    
    %{
    % --- Plotting ---
    
    thetas = linspace(0, 2*pi, 100);
    
    figure(1)
    hold on
    
    plot(r_earth_sun*cos(thetas), r_earth_sun*sin(thetas), 'b', 'LineWidth', 2); hold on
    plot(r_mars_sun*cos(thetas), r_mars_sun*sin(thetas), 'r', 'LineWidth', 2)
    
    plot(R11(1), R11(2), '.b', 'MarkerSize', 20);
    plot(R21(1), R21(2), '.r', 'MarkerSize', 20);
    
    plot(R22(1), R22(2), '.b', 'MarkerSize', 20);
    plot(R23(1), R23(2), '.r', 'MarkerSize', 20);
    
    for k = 1:width(v_D1)
        plot_lambert_arc_time(R11, v_D1(:, k), mu_sun, tof1, 'k')
    end
    
    for k = 1:width(v_D2)
        plot_lambert_arc_time(R12, v_D2(:, k), mu_sun, tof2, 'k')
    end
    
    for k = 1:width(v_D3)
        plot_lambert_arc_time(R13, v_D3(:, k), mu_sun, tof3, 'k')
    end
    
    
    xlim([-1, 1] * r_mars_sun*2)
    ylim([-1, 1] * r_mars_sun*2)
    
    axis square
    drawnow
    pause(0.1)
    %}
end
end

%% Plotting
figure(1)
clf

surf(ns, leg1_t / S, num_paths, 'LineStyle', 'none')

xlabel("n")
ylabel("Leg1 tof (S/n)")
zlabel("Valid Paths")

function num_paths = calc_num_paths(v_A_list, v_D_list, mu, rmin)
    n1 = width(v_A_list);
    n2 = width(v_D_list);

    num_paths = 0;

    for i1 = 1:n1
        for i2 = 1:n2
            rp = grav_assist_rp(v_A_list(:, i1), v_D_list(:, i2), mu);
            if rp > rmin
                num_paths = num_paths + 1;
            end
        end
    end
end

function plot_lambert_arc(R1, R2, V1, mu, color)
    % --- Compute orbital elements ---
    h = cross(R1, V1);
    h_hat = h / norm(h);
    e_vec = cross(V1, h)/mu - R1/norm(R1);
    e = norm(e_vec);

    % Semi-latus rectum and parameter
    p = norm(h)^2 / mu;

    % --- Define orbital plane basis ---
    p_hat = e_vec / e;
    q_hat = cross(h_hat, p_hat); % completes right-handed basis

    % --- Compute true anomalies in that plane ---
    theta1 = atan2(dot(R1, q_hat), dot(R1, p_hat));
    theta2 = atan2(dot(R2, q_hat), dot(R2, p_hat));

    % --- Choose the *short-way* arc ---
    dtheta = mod(theta2 - theta1, 2*pi);
    if dtheta > pi
        dtheta = dtheta - 2*pi;
    end

    % If the Lambert solution is mirrored, force the correct direction:
    if sign(dot(cross(R1, R2), h_hat)) < 0
        h_hat = -h_hat;
        q_hat = -q_hat;
    end

    % --- Generate arc ---
    thetas = linspace(theta1, theta1 + dtheta, 400);
    r = p ./ (1 + e*cos(thetas));

    % Construct inertial coordinates properly
    r_orb = p_hat * (r .* cos(thetas)) + q_hat * (r .* sin(thetas));

    % --- Plot ---
    plot3(r_orb(1,:), r_orb(2,:), r_orb(3,:), color, 'LineWidth', 2);
    plot3([R1(1) R2(1)], [R1(2) R2(2)], [R1(3) R2(3)], 'ko', 'MarkerFaceColor','k');
    axis equal; grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
end