clear; clc; close all

params = solar_orbit_params();

mu_sun = 132712000000;                    % [km^3 s^-2]
num_points = 100;

dep_start_date = datetime(2020, 1, 1);
dep_end_date   = datetime(2021, 1, 1);

arr_start_date = datetime(2020, 9, 1);
arr_end_date   = datetime(2022, 1, 1);

dep_dates = dep_start_date + days(linspace(0, days(dep_end_date - dep_start_date), num_points));
arr_dates = arr_start_date + days(linspace(0, days(arr_end_date - arr_start_date), num_points));

dv_all = nan(num_points, num_points);

z_lim = 1200;
num_z = 250;   % 1000 is expensive inside a 100x100 grid

parfor dep_idx = 1:num_points
    dep_date = dep_dates(dep_idx);
    r_dep_earth = planet_state_from_datetime(params, 'EarthMoonBarycenter', dep_date).r_ecl;
    v_dep_earth = planet_state_from_datetime(params, 'EarthMoonBarycenter', dep_date).v_ecl;

    for arr_idx = 1:num_points
        arr_date = arr_dates(arr_idx);
        r_arr_mars = planet_state_from_datetime(params, 'Mars', arr_date).r_ecl;
        v_arr_mars = planet_state_from_datetime(params, 'Mars', arr_date).v_ecl;

        tof = seconds(arr_date - dep_date); % s
        if tof <= 0
            continue; % skip non-physical pairs
        end

        dv_all(dep_idx, arr_idx) = min_dv_lambert_universal(r_dep_earth/1e3, v_dep_earth/1e3, r_arr_mars/1e3, v_arr_mars/1e3, tof, mu_sun, z_lim, num_z);

        fprintf("%d/%d, %d/%d\n", dep_idx, num_points, arr_idx, num_points);
    end
end

%% PLOT

contour(days(dep_dates - dep_dates(1)), days(arr_dates - arr_dates(1)), dv_all.', 0:1:20);
colorbar