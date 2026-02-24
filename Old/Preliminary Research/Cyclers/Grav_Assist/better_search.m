clear; clc; close all;

[mu_sun, mu_earth, mu_mars,... 
    r_earth_sun, r_mars_sun,... 
    T_earth, T_mars, S, angle_turn] = get_earth_mars_params();

r_mars = 3396; % [km]
r_earth = 6378;

% Two parameters: tof and position of mars

R11 = [r_earth_sun; 0; 0];
R22 = r_earth_sun * [cos(2*pi*S/T_earth); sin(2*pi*S/T_earth); 0];

num_theta = 1000;
theta_mars = linspace(0, 2*pi, num_theta);

num_tof = 1000;
tof_mars = linspace(0, S, num_tof);

% zlim = 1000;
% numz = 100;

rps = zeros(num_theta, num_tof, 2);

orbit_options = {"1A", "2A", "1B", "2B"};
num_options = numel(orbit_options);

for theta_idx=1:num_theta
    theta = theta_mars(theta_idx);
    R21 = r_mars_sun * [cos(theta); sin(theta); 0];
    R12 = R21;

    for tof_idx = 1:num_tof
        tof1 = tof_mars(tof_idx);
        tof2 = S - tof1;

        R13 = R22;
        R23 = r_mars_sun * [cos(2*pi*(S+tof1)/T_mars); sin(2*pi*(S+tof1)/T_mars); 0];

        %{
        [v_D1, v_A1] = lambert(R11, R21, tof1, mu_sun, zlim, numz);
        v_D1 = cell2mat(v_D1); v_A1 = cell2mat(v_A1);
        
        [v_D2, v_A2] = lambert(R12, R22, tof2, mu_sun, zlim, numz);
        v_D2 = cell2mat(v_D2); v_A2 = cell2mat(v_A2);

        if isempty(v_D2) || isempty(v_D1)
            continue
        end

        num1 = width(v_D1);
        num2 = width(v_D2);

        max_rp = 0;

        for idx1 = 1:num1
            for idx2 = 1:num2
                rp = grav_assist_rp(v_A1(:, idx1), v_D2(:, idx2), mu_mars);
                max_rp = max(max_rp, rp);
            end
        end
        max_rps(theta_idx, tof_idx) = max_rp;
        %}

        max_rp = 0;
        for idx1 = 1:num_options
            for idx2 = 1:num_options
                type_1 = orbit_options{idx1};
                type_2 = orbit_options{idx2};
                [~, ~, v_D1, v_A1] = lambert_solver(R11, R21, tof1, mu_sun, type_1);
                [~, ~, v_D2, v_A2] = lambert_solver(R12, R22, tof2, mu_sun, type_2);
                [~, ~, v_D3, v_A3] = lambert_solver(R13, R23, tof1, mu_sun, type_1);

                rp1 = grav_assist_rp(v_A1, v_D2, mu_mars);
                rp2 = grav_assist_rp(v_A2, v_D3, mu_earth);

                if rp1 > r_mars
                    rps(theta_idx, tof_idx, 1) = 1e10;
                end

                if rp2 > r_earth
                    rps(theta_idx, tof_idx, 2) = 1e10;
                end

                % rps(theta_idx, tof_idx, 1) = real(rp1);
                % rps(theta_idx, tof_idx, 2) = real(rp2);
            end
        end

        % max_rps(theta_idx, tof_idx) = max_rp;

        fprintf("Theta Idx: %d, Tof Idx, %d\n", theta_idx, tof_idx)
    end
end


%% Plotting
figure(1)
clf

% max_rps_new = (max_rps > 6378) * 1;

rps_new = min(rps, r_earth * 5);

surf(theta_mars, tof_mars/S, rps_new(:, :, 1), 'LineStyle', 'none')
xlabel("Theta")
ylabel("tof (S)")
zlabel("rp")
zlim([r_mars, r_earth * 5])

figure(2)
clf
surf(theta_mars, tof_mars/S, rps_new(:, :, 2), 'LineStyle', 'none')
xlabel("Theta")
ylabel("tof (S)")
zlabel("rp")
zlim([r_earth, r_earth * 5])
