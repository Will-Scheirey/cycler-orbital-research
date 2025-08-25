clear; clc; close all

mu_venus = 324900;     % [km^3 s^-2] 
r_venus = 6051.8;      % [km]
r_venus_sun = 108.2e6; % [km]

mu_earth = 398600;     % [km^3 s^-2]
r_earth = 6378;        % [km]
r_earth_sun = 149.6e6; % [km]

mu_mars = 42828;       % [km^3 s^-2] 
r_mars = 3396;         % [km]
r_mars_sun = 227.9e6;  % [km]

mu_sun = 132712000000; % [km^3 s^-2]

v_venus_sun = sqrt(mu_sun / r_venus_sun); % [km/s]
v_earth_sun = sqrt(mu_sun / r_earth_sun); % [km/s]
v_mars_sun  = sqrt(mu_sun / r_mars_sun);  % [km/s]

num_v_inf = 5;
v_infs = linspace(1, 10, num_v_inf);

num_alphas = 100;
alphas = linspace(0, 180, num_alphas);

points_venus = generate_tisserand_points(v_infs, alphas, mu_sun, v_venus_sun, r_venus_sun);
for i = 1:num_v_inf
    rp = squeeze(points_venus(i, :, 1));
    T  = squeeze(points_venus(i, :, 2));
    plot(rp, T, '-k'); hold on
end

points_earth = generate_tisserand_points(v_infs, alphas, mu_sun, v_earth_sun, r_earth_sun);
for i = 1:num_v_inf
    rp = squeeze(points_earth(i, :, 1));
    T  = squeeze(points_earth(i, :, 2));
    plot(rp, T, '-b'); hold on
end

points_mars = generate_tisserand_points(v_infs, alphas, mu_sun, v_mars_sun, r_mars_sun);
for i = 1:num_v_inf
    rp = squeeze(points_mars(i, :, 1));
    T  = squeeze(points_mars(i, :, 2));
    plot(rp, T, '-r'); hold on
end

function points = generate_tisserand_points(v_infs, alphas, mu_2, v1_2, r_1_2)
    num_v_inf = length(v_infs);
    num_alphas = length(alphas);

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
            points(i, j, :) = [rp_2 / 149.6e6, epsilon_2];
        end
    end
end
