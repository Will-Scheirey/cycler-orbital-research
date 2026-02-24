clear; clc; close all;

mu_sun = 132712000000;                    % [km^3 s^-2]

r_earth_sun = 149.6e6;                    % [km]
v_earth_sun = sqrt(mu_sun / r_earth_sun);  % [km s^-1]

r_mars_sun = 227.9e6;                     % [km]
v_mars_sun = sqrt(mu_sun / r_mars_sun);  % [km s^-1]

r_dep = [r_earth_sun, 0, 0];

num_v = 1000;
min_v = 0;
max_v = 5;
vx = linspace(min_v, max_v, num_v);
vy = linspace(min_v, max_v, num_v);

dv_tot = nan(num_v, num_v);

for i = 1:num_v
    for n = 1:num_v
        dv_dep = norm([vx(i), vy(n)]);

        v_dep = [vx(i), vy(n) + v_earth_sun, 0];
        eps = norm(v_dep)^2/2 - mu_sun/r_earth_sun;

        h_vec = cross(r_dep, v_dep);
        h = norm(h_vec);
        a = -mu_sun/(2*eps);

        e = 1 - h^2/(mu_sun * a);

        if h^2/mu_sun * 1/(1-e) < r_mars_sun, continue; end

        theta_encounter_mars = acos((h^2/(mu_sun*r_mars_sun) - 1)/e);

        vr = mu_sun/h * e * sin(theta_encounter_mars);
        vt = mu_sun/h * (1 + e*cos(theta_encounter_mars));

        v_arr = [vr, vt];
        v_mars = [0, v_mars_sun];

        dv_arr = norm(v_arr - v_mars);
        
        dv_tot(i, n) = dv_arr + dv_dep;
    end
end

contourf(vx, vy, dv_tot)
colorbar