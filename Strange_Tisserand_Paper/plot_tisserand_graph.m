clear; clc; close all

% Reference Paper: Graphical Method for Gravity-Assist Trajectory Design
% by Nathan J. Strange and James M. Longuski

mu_mercury    = 22030;                    % [km^3 s^-2] 
r_mercury     = 2440;                     % [km]
r_mercury_sun = 57.91e6;                  % [km]

mu_venus = 324900;                        % [km^3 s^-2] 
r_venus = 6051.8;                         % [km]
r_venus_sun = 108.2e6;                    % [km]

mu_earth = 398600;                        % [km^3 s^-2]
r_earth = 6378;                           % [km]
r_earth_sun = 149.6e6;                    % [km]

mu_mars = 42828;                          % [km^3 s^-2] 
r_mars = 3396;                            % [km]
r_mars_sun = 227.9e6;                     % [km]

mu_jupiter = 126686000;                   % [km^3 s^-2] 
r_jupiter = 71490;                        % [km]
r_jupiter_sun = 778.6e6;                  % [km]

mu_saturn = 37931000;                     % [km^3 s^-2] 
r_saturn = 60270;                         % [km]
r_saturn_sun = 1.433e9;                   % [km]

mu_sun = 132712000000;                    % [km^3 s^-2]

num_v_inf = 10;
v_infs = linspace(1, 30, num_v_inf);

num_alphas = 20;
alphas = linspace(0, 180, num_alphas);

num_ks = 10;
ks = linspace(0, 50, num_ks);

points_mercury = generate_tisserand_points(v_infs, alphas, ks, mu_sun, r_mercury_sun);
points_venus = generate_tisserand_points(v_infs, alphas, ks, mu_sun, r_venus_sun);
points_earth = generate_tisserand_points(v_infs, alphas, ks, mu_sun, r_earth_sun);
points_mars  = generate_tisserand_points(v_infs, alphas, ks, mu_sun, r_mars_sun);
points_jupiter = generate_tisserand_points(v_infs, alphas, ks, mu_sun, r_jupiter_sun);
points_saturn = generate_tisserand_points(v_infs, alphas, ks, mu_sun, r_saturn_sun);

%% Plot
figure(1)
clf
plot_tisserand(points_mercury, mu_mercury, r_mercury + 500, v_infs, alphas, ks, [0.1, 0.1, 0.1], 'Mercury')
plot_tisserand(points_venus, mu_venus, r_venus + 500, v_infs, alphas, ks, [0.9, 0.75, 0.4], 'Venus')
plot_tisserand(points_earth, mu_earth, r_earth + 500, v_infs, alphas, ks, [0.1, 0.1, 0.9], 'Earth')
plot_tisserand(points_mars, mu_mars, r_mars + 500, v_infs, alphas, ks, [1, 0.1, 0.1], 'Mars')
plot_tisserand(points_jupiter, mu_jupiter, r_jupiter*5, v_infs, alphas, ks, [1, 0.5, 0.1], 'Jupiter')
plot_tisserand(points_saturn, mu_saturn, r_jupiter*5, v_infs, alphas, ks, [0.9, 0.8, 0.7], 'Saturn')

legend
xlabel("Periapsis Radius (AU)")
ylabel("Specific Energy (km^2/s^2)")
zlabel("Inclination (deg)")
title("Tisserand Plots for Multiple Planets")

% Allows us to zoom in farther before things start disappearing visually
gca.CameraViewAngle = 1.5;

function plot_tisserand(points, mu, rp_min, v_infs, alphas, ks, color, name)

    num_points = height(points);

    %{
    for i = 1:num_points
        surf(alphas, ks, squeeze(points(i, :, :, 3)), 'FaceColor', color, 'LineStyle', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off'); hold on
    end
    %}
    
    num_ks = length(ks);
    num_alphas = length(alphas);

    num_points = height(points);
    for i = 1:num_points
        delta_max = delta_from_rp(mu, v_infs(i), rp_min);
        alpha_step = find_closest_idx(alphas, delta_max);
        curr_step = 0;

        for n = 1:num_alphas
            arg1 = squeeze(points(i, n, :, 1)); % rp
            arg2 = squeeze(points(i, n, :, 2)); % Specific Energy
            arg3 = squeeze(points(i, n, :, 3)); % Period
            arg4 = squeeze(points(i, n, :, 4)); % Inclination

            if i == 1 && n == 1
                plot(arg1, arg2, 'Color', color, 'LineStyle', '-', 'DisplayName', name, 'LineWidth', 0.5);
            else
                plot(arg1, arg2, 'Color', color, 'LineStyle', '-','HandleVisibility', 'off', 'LineWidth', 0.5);
            end

            if floor(alpha_step) > curr_step
                plot(arg1, arg2, 'Color', color, 'Marker', '.', 'MarkerSize', 5, 'HandleVisibility', 'off');
            end

            hold on
        end
    end
end

function delta = delta_from_rp(mu, v_inf, rp)
    e = 1 + (rp*v_inf^2)/mu;
    delta = 2*asind(1/e);
end

function rp = rp_from_delta(mu, v_inf, delta_deg)
    e = 1/sind(delta_deg/2);
    rp = (e - 1) * (mu / v_inf^2);
end

function points = generate_tisserand_points(v_infs, alphas, ks, mu_2, r_1_2)
    AU = 149.6e6;
    num_v_inf = length(v_infs);
    num_alphas = length(alphas);
    num_ks = length(ks);

    v1_2 = sqrt(mu_2 / r_1_2);

    v_vec_1_2 = [0, v1_2, 0];

    points = zeros(num_v_inf, num_alphas, num_ks, 4);

    for i=1:num_v_inf
        v_inf_1 = v_infs(i);
        for j=1:num_alphas
            for n = 1:num_ks
                alpha = alphas(j);
                k = ks(n);
    
                r_vec_2 = [r_1_2, 0, 0];
                v_inf_vec_1 = v_inf_1 * [sind(alpha)*sind(k), cosd(alpha)*sind(k), cosd(k)];

                v_vec_2 = v_vec_1_2 + v_inf_vec_1;
                v_2 = norm(v_vec_2);
                
                h_vec_2 = cross(r_vec_2, v_vec_2);
                h_2 = norm(h_vec_2);
        
                epsilon_2 = v_2^2/2 - mu_2 / r_1_2;
                a_2 = -mu_2 / (2*epsilon_2);
                
                e_2 = sqrt(1 - h_2^2/mu_2 * 1/a_2);
            
                rp_2 = h_2^2/mu_2 * 1/(1 + e_2);
                period = 2*pi*sqrt(a_2^3/mu_2);

                if imag(period) ~= 0 || period > 1e9
                    period = 0;
                end

                inclination = acosd(h_vec_2(3) / h_2);

                points(i, j, n, :) = [rp_2 / AU, epsilon_2, period / 86400 / 365, inclination];
            end
        end
    end
end

function idx = find_closest_idx(arr, target)
    [~, idx] = min(abs(arr - target));
end