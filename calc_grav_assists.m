clear; clc; close all

mu = 132712000000;

AU = 1.496e+8;

r_b1 = AU;
v_b1 = sqrt(mu / r_b1);

mu_b2 = 42828;
r_b2 = 2.2737e8;
radius_b2 = 3396;

mu_b2 = 126686000;
r_b2 = 778.6e6;
radius_b2 = 80000;

soi_b2 = r_b2 * (mu_b2 / mu)^(2/5);

min_rp = radius_b2 + 100;
max_rp = soi_b2;

dv_min = sqrt(2*mu / r_b1 - 2*mu / (r_b1 + r_b2)) - sqrt(mu/r_b1);
dv_max = dv_min + 10;

v_b2 = sqrt(mu / r_b2);
h_b2 = r_b2*v_b2;
e_b2 = 0;

[S, angle_turn] = calc_synodic(mu, r_b1, r_b2);

n = 1;
S = S * n;
angle_turn = mod(angle_turn * n, 2*pi);

num_points = 100;

t_all = linspace(-pi, pi, num_points);
v_all = linspace(dv_min, dv_max, num_points);
rp_all = linspace(min_rp, max_rp, num_points);

all_results = cell(num_points, num_points, num_points);
costs = nan(num_points, num_points, num_points);

r0 = [r_b1; 0; 0];

for t_idx = 1:num_points
    t = t_all(t_idx);
    for v_idx = 1:num_points
        v_mag = v_all(v_idx);
        v0x = v_mag * sin(t);
        v0y =  v_b1 + v_mag * cos(t);

        v0 = [v0x; v0y; 0];

        r_out = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, 1e8, 1, 1);

        if isempty(r_out)
            continue;
        end

        for rp_idx = 1:num_points
            rp = rp_all(rp_idx);

            all_results_temp = calc_all_results(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp);
            all_results_temp = calc_all_costs(all_results_temp, S, angle_turn, t, v_mag);

            result = all_results_temp.best;
            costs(t_idx, v_idx, rp_idx)       = result.cost;
            all_results{t_idx, v_idx, rp_idx} = result;
        end
    end
end

%% PLOT

costs_temp = costs(:);
results_temp = all_results(:);

figure(1)
clf
best_costs = min(costs, [], 3, "omitmissing");   % 20x20

surf(v_all, t_all, best_costs, 'LineStyle','none'); % , 'Marker', '.', 'MarkerSize', 5
xlabel("V Mag")
ylabel("t")
zlabel("Cost")

% return

min_cost = min(costs_temp);
min_result = results_temp{min_cost == costs_temp};

for z = 1:length(results_temp)
min_result = results_temp{z};

if isempty(min_result)
    continue
end

if isinf(min_result.cost)
    continue
end

return

figure(2)

t     = min_result.t;
v_mag = min_result.v0;
rp    = min_result.rp;

v0x = v_mag * sin(t);
v0y =  v_b1 + v_mag * cos(t);
v0 = [v0x; v0y; 0];

[r_sc, h_sc, e_sc, p_sc, w_sc, nu_sc] = calc_orbit(r0, v0, mu);

[r_out, h_out, e_out, p_out, w_out, nu_out, theta_intersect, tof1] = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, min_result.sign, min_result.angle);

theta_sc = linspace(0, theta_intersect, 1000);
nu_sc = theta_sc - w_sc;
r1 = p_sc ./ (1 + e_sc*cos(nu_sc));

theta_plot = linspace(0,2*pi,1000);

theta_b1 = calc_intersect(p_out, AU, e_out, w_out, min_result.best_angle);
theta_out = linspace(theta_intersect, theta_b1, 1000);
nu_out = theta_out - w_out;
r2 = p_out ./ (1 + e_out*cos(nu_out));

figure(2)
clf
plot(r2.*cos(theta_out), r2.*sin(theta_out), '-g', 'LineWidth', 1.5); hold on
plot(r1.*cos(theta_sc), r1.*sin(theta_sc), '--k', 'LineWidth', 1.5); hold on

plot(r_b1*cos(theta_plot), r_b1*sin(theta_plot), '-b'); hold on
plot(r_b2*cos(theta_plot), r_b2*sin(theta_plot), '-r')

plot(r_b2*cos(theta_intersect), r_b2*sin(theta_intersect), 'mx')

plot(AU*cos(theta_b1), AU*sin(theta_b1), 'mx')

axis equal
xlim([-1,1]*r_b2 * 4)
ylim([-1,1]*r_b2 * 4)
drawnow;
title(sprintf("Cost = %0.5f, Rp = %0.1f", min_result.cost, min_result.rp))

fprintf("%d / %d\n", z, length(results_temp))
pause(0.05)
end

function all_results = calc_all_costs(all_results, S, angle_turn, t, v_mag)
        if ~isempty(all_results.result1)
            all_results.result1.t    = t;
            all_results.result1.v0   = v_mag;
            [cost, best_angle] = calc_cost(all_results.result1, S, angle_turn);
            all_results.result1.cost = cost;
            all_results.result1.best_angle = best_angle;
        end

        if ~isempty(all_results.result2)
            all_results.result2.t    = t;
            all_results.result2.v0   = v_mag;
            [cost, best_angle] = calc_cost(all_results.result2, S, angle_turn);
            all_results.result2.cost = cost;
            all_results.result2.best_angle = best_angle;        
        end

        if ~isempty(all_results.result3)
            all_results.result3.t    = t;
            all_results.result3.v0   = v_mag;
            [cost, best_angle] = calc_cost(all_results.result3, S, angle_turn);
            all_results.result3.cost = cost;
            all_results.result3.best_angle = best_angle;        
        end

        if ~isempty(all_results.result4)
            all_results.result4.t    = t;
            all_results.result4.v0   = v_mag;
            [cost, best_angle] = calc_cost(all_results.result4, S, angle_turn);
            all_results.result4.cost = cost;
            all_results.result4.best_angle = best_angle;       
        end

        if isempty(all_results.result1)
            all_results.result1.cost = inf;
        end
        best = all_results.result1;
        if ~isempty(all_results.result2), if all_results.result2.cost < best.cost, best = all_results.result2; end; end
        if ~isempty(all_results.result3), if all_results.result3.cost < best.cost, best = all_results.result3; end; end
        if ~isempty(all_results.result4), if all_results.result4.cost < best.cost, best = all_results.result4; end; end

        if isempty(best)
            disp("!")
        end
        all_results.best = best;
end

function all_results = calc_all_results(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp)
    result1 = calc_result(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, 1, 1);
    if ~isempty(result1)
        result1.sign = 1;
        result1.angle = 1;
    end

    result2 = calc_result(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, 1, 2);
    if ~isempty(result2)
        result2.sign = 1;
        result2.angle = 2;
    end

    result3 = calc_result(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, -1, 1);
    if ~isempty(result3)
        result3.sign = -1;
        result3.angle = 1;
    end

    result4 = calc_result(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, -1, 2);
    if ~isempty(result4)
        result4.sign = -1;
        result4.angle = 2;
    end

    all_results = struct('result1', result1, 'result2', result2, 'result3', result3, 'result4', result4);
end

function [cost, best_angle] = calc_cost(result, S, angle_turn)    

    calc = @(tof_tot, theta_b1) (abs(tof_tot - S)/S)^2 + (abs(theta_b1 - angle_turn)/angle_turn)^2;
    
    cost1 = calc(result.tof_tot_1, result.theta_b1_1);
    cost2 = calc(result.tof_tot_2, result.theta_b1_2);

    cost = min(cost1, cost2);
    best_angle = 1;
    if cost2 < cost1, best_angle = 2; end
end

function result = calc_result(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, sign, angle)

result = [];

[r_out, h_out, e_out, p_out, w_out, nu_out, theta_b2, tof_12] = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, sign, angle);

if isempty(r_out) || e_out >= 1
    return;
end

rp_helio = p_out / (1 + e_out);

if rp_helio > norm(r0)
    return
end

theta_b1_1 = calc_intersect(p_out, norm(r0), e_out, w_out, 1);

if isempty(theta_b1_1)
    return
end

tof_21_1 = calc_tof(theta_b2, theta_b1_1, w_out, e_out, p_out, mu);
tof_tot_1 = tof_12 + tof_21_1;

theta_b1_2 = calc_intersect(p_out, norm(r0), e_out, w_out, 2);
tof_21_2 = calc_tof(theta_b2, theta_b1_2, w_out, e_out, p_out, mu);
tof_tot_2 = tof_12 + tof_21_2;

result = struct( ...
    'rp',         rp, ...
    'theta_b2',   theta_b2, ...
    'theta_b1_1', theta_b1_1, ...
    'tof_tot_1',  tof_tot_1, ...
    'theta_b1_2', theta_b1_2, ...
    'tof_tot_2',  tof_tot_2);
end

function [r_out, h_out, e_out, p_out, w_out, nu_out, theta_intersect, tof1] = calc_grav_assist(r0_vec, v0_vec, mu, mu_b2, r_b2, h_b2, e_b2, rp, sign, angle)

[~, h_sc, e_sc, p_sc, w_sc, ~] = calc_orbit(r0_vec, v0_vec, mu);
a_sc = p_sc / (1 - e_sc^2);

r_out = [];
h_out = [];
e_out = [];
p_out = [];
w_out = [];
nu_out = [];
theta_intersect = [];
tof1 = [];

theta_intersect = calc_intersect(p_sc, r_b2, e_sc, w_sc, angle);

if isempty(theta_intersect)
    return
end

nu_int = theta_intersect - w_sc;

% Spacecraft velocity at intersect
vr_sc = mu/h_sc * e_sc * sin(nu_int);
vt_sc = mu/h_sc * (1 + e_sc * cos(nu_int));

% Body 2 velocity at intersect
vr_b2 = mu/h_b2 * e_b2 * sin(theta_intersect);
vt_b2 = mu/h_b2 * (1 + e_b2 * cos(theta_intersect));

% V infinity
v_inf_in = [vr_sc - vr_b2; vt_sc - vt_b2];
v_inf    = norm(v_inf_in);

% Hyperbolic eccentricity around body 2
e  = 1 + rp * v_inf^2 / mu_b2;
% Turn angle
delta = sign * 2 * asin(1 / e);

% Rotate V infinity
R = [cos(delta) -sin(delta); sin(delta) cos(delta)];

v_inf_out = R * v_inf_in;

% Calculate heliocentric velocity
vr_out = vr_b2 + v_inf_out(1);
vt_out = vt_b2 + v_inf_out(2);

v_out = hypot(vr_out, vt_out);

rx_out = r_b2*cos(theta_intersect);
ry_out = r_b2*sin(theta_intersect);

vx_out = vr_out*cos(theta_intersect) - vt_out*sin(theta_intersect);
vy_out = vr_out*sin(theta_intersect) + vt_out*cos(theta_intersect);

r_vec_out = [rx_out, ry_out, 0];
v_vec_out = [vx_out, vy_out, 0];

[r_out, h_out, e_out, p_out, w_out, nu_out] = calc_orbit(r_vec_out, v_vec_out, mu);

% TOF
E = 2 * atan(sqrt((1 - e_sc) / (1 + e_sc)) * tan(nu_int / 2));
Me = E - e_sc * sin(E);
tof1 = sqrt(a_sc^3 / mu) * Me;
end

function theta_intersect = calc_intersect(p_sc, r_b2, e_sc, w_sc, angle)
theta_intersect = [];
% Intersect
c = (p_sc/r_b2 - 1)/e_sc;
if abs(c) > 1, return; end
nu_int = acos(c);

if imag(nu_int) ~= 0, return; end

% Calculate both points
phi_a = w_sc + nu_int;
phi_b = w_sc - nu_int;

phi_a = mod(phi_a, 2*pi);
phi_b = mod(phi_b, 2*pi);

% Choose one
if angle == 1
    theta_intersect = min(phi_a, phi_b);
else
    theta_intersect = max(phi_a, phi_b);
end

end


function tof = calc_tof(theta_b2, theta_b1, w_out, e_out, p_out, mu)

    nu2 = theta_b2 - w_out;
    nu1 = theta_b1 - w_out;

    a = p_out / (1 - e_out^2);

    E2 = 2 * atan( sqrt((1 - e_out) / (1 + e_out)) * tan(nu2/2) );
    E1 = 2 * atan( sqrt((1 - e_out) / (1 + e_out)) * tan(nu1/2) );

    M2 = E2 - e_out*sin(E2);
    M1 = E1 - e_out*sin(E1);

    dM = M1 - M2;
    if dM < 0
        dM = dM + 2*pi;
    end

    tof = sqrt(a^3 / mu) * dM;
end