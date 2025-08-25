clear; clc; close all

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

num_v_inf = 5;
v_infs = linspace(1, 10, num_v_inf);

num_alphas = 100;
alphas = linspace(0, 180, num_alphas);

points_mercury = generate_tisserand_points(v_infs, alphas, mu_sun, r_mercury_sun);
points_venus = generate_tisserand_points(v_infs, alphas, mu_sun, r_venus_sun);
points_earth = generate_tisserand_points(v_infs, alphas, mu_sun, r_earth_sun);
points_mars  = generate_tisserand_points(v_infs, alphas, mu_sun, r_mars_sun);
points_jupiter = generate_tisserand_points(v_infs, alphas, mu_sun, r_jupiter_sun);
points_saturn = generate_tisserand_points(v_infs, alphas, mu_sun, r_saturn_sun);

plot_tisserand(points_mercury, '-k', 'Mercury')
plot_tisserand(points_venus, '-y', 'Venus')
plot_tisserand(points_earth, '-b', 'Earth')
plot_tisserand(points_mars, '-r', 'Mars')
plot_tisserand(points_jupiter, '-m', 'Jupiter')
plot_tisserand(points_saturn, '-g', 'Saturn')

legend
xlabel("Perihelion (AU)")
ylabel("Orbital Specific Energy (km^2/s^2)")
title("Tisserand Plots for Multiple Planets")

function plot_tisserand(points, line, name)
    num_points = height(points);
    semilogx(points(1,:,1), points(1,:,2), line, 'DisplayName', name, 'LineWidth', 2); hold on
    for i = 2:num_points
        rp = squeeze(points(i, :, 1));
        T  = squeeze(points(i, :, 2));
        semilogx(rp, T, line, 'HandleVisibility', 'off', 'LineWidth', 2);
    end
end

function points = generate_tisserand_points(v_infs, alphas, mu_2, r_1_2)
    AU = 149.6e6
    num_v_inf = length(v_infs);
    num_alphas = length(alphas);

    v1_2 = sqrt(mu_2 / r_1_2);

    v_vec_1_2 = [0, v1_2, 0];

    points = zeros(num_v_inf, num_alphas, 2);

    for i=1:num_v_inf
        v_inf_1 = v_infs(i);
    
        for j=1:num_alphas
            alpha = alphas(j);
    
            r_vec_2 = [r_1_2, 0, 0];
            v_inf_vec_1 = v_inf_1 * [sind(alpha), cosd(alpha), 0];
    
            v_vec_2 = v_vec_1_2 + v_inf_vec_1;
            v_2 = norm(v_vec_2);
            
            h_vec_2 = cross(r_vec_2, v_vec_2);
            h_2 = norm(h_vec_2);
    
            epsilon_2 = v_2^2/2 - mu_2 / r_1_2;
            a_2 = -mu_2 / (2*epsilon_2);
            
            e_2 = sqrt(1 - h_2^2/mu_2 * 1/a_2);
        
            rp_2 = h_2^2/mu_2 * 1/(1 + e_2);
            period = 2*pi*sqrt(a_2^3/mu_2);
            % points(i, j, :) = [rp_2 / 149.6e6, period / 86400];
            points(i, j, :) = [rp_2, epsilon_2];
            % points(i,j,1) = points(i,j,1) / AU;
        end
    end
end
