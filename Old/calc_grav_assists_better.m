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

num_points = 50;

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

        ok = false;

        [valid0, ~, ~, ~, ~, ~] = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, 1e8,  1, 1);
        ok = ok | any(valid0);

        [valid0, ~, ~, ~, ~, ~] = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, 1e8,  1, 2);
        ok = ok | any(valid0);

        [valid0, ~, ~, ~, ~, ~] = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, 1e8, -1, 1);
        ok = ok | any(valid0);

        [valid0, ~, ~, ~, ~, ~] = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, 1e8, -1, 2);
        ok = ok | any(valid0);

        if ~ok
            continue;
        end

        all_results_temp = calc_all_results(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp_all);
        all_results_temp = calc_all_costs(all_results_temp, S, angle_turn, t, v_mag);

        result = all_results_temp.best;
        
       if isempty(result)
            continue;
        end

        costs(t_idx, v_idx, :)       = result.cost;
        for rp_idx = 1:num_points
            if isfinite(result.cost(rp_idx)) && result.cost(rp_idx) < inf
                all_results{t_idx, v_idx, rp_idx} = struct( ...
                    't', t, ...
                    'v0', v_mag, ...
                    'rp', rp_all(rp_idx), ...
                    'cost', result.cost(rp_idx), ...
                    'sign', result.sign(rp_idx), ...
                    'angle', result.angle(rp_idx), ...
                    'best_angle', result.best_angle(rp_idx) );
            else
                all_results{t_idx, v_idx, rp_idx} = [];
            end
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

return

[min_cost, lin_idx] = min(costs_temp, [], "omitnan");
min_result = results_temp{lin_idx};

for z = 1:length(results_temp)
min_result = results_temp{z};

if isempty(min_result)
    continue
end

if isinf(min_result.cost)
    continue
end

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
pause(5)
end

function all_results = calc_all_costs(all_results, S, angle_turn, t, v_mag)

    % --- evaluate each branch (each returns Nx1 cost, Nx1 best_angle) ---
    [c1,a1] = deal([]);
    [c2,a2] = deal([]);
    [c3,a3] = deal([]);
    [c4,a4] = deal([]);

    rp = 0;

    if ~isempty(all_results.result1)
        all_results.result1.t  = t;  all_results.result1.v0 = v_mag;
        [c1,a1] = calc_cost_vec(all_results.result1, S, angle_turn);
        all_results.result1.cost = c1;  all_results.result1.best_angle = a1;
        rp = all_results.result1.rp;
    end
    if ~isempty(all_results.result2)
        all_results.result2.t  = t;  all_results.result2.v0 = v_mag;
        [c2,a2] = calc_cost_vec(all_results.result2, S, angle_turn);
        all_results.result2.cost = c2;  all_results.result2.best_angle = a2;
        rp = all_results.result2.rp;
    end
    if ~isempty(all_results.result3)
        all_results.result3.t  = t;  all_results.result3.v0 = v_mag;
        [c3,a3] = calc_cost_vec(all_results.result3, S, angle_turn);
        all_results.result3.cost = c3;  all_results.result3.best_angle = a3;
        rp = all_results.result3.rp;
    end
    if ~isempty(all_results.result4)
        all_results.result4.t  = t;  all_results.result4.v0 = v_mag;
        [c4,a4] = calc_cost_vec(all_results.result4, S, angle_turn);
        all_results.result4.cost = c4;  all_results.result4.best_angle = a4;
        rp = all_results.result4.rp;
    end

    % --- elementwise best across branches ---
    % Make missing branches = +Inf vectors of the correct length
    N = 0;
    if ~isempty(c1), N = numel(c1);
    elseif ~isempty(c2), N = numel(c2);
    elseif ~isempty(c3), N = numel(c3);
    elseif ~isempty(c4), N = numel(c4);
    end
    if N == 0
        all_results.best = [];
        return;
    end

    if isempty(c1), c1 = inf(N,1); a1 = ones(N,1); end
    if isempty(c2), c2 = inf(N,1); a2 = ones(N,1); end
    if isempty(c3), c3 = inf(N,1); a3 = ones(N,1); end
    if isempty(c4), c4 = inf(N,1); a4 = ones(N,1); end

    c1(~isfinite(c1)) = inf;
    c2(~isfinite(c2)) = inf;
    c3(~isfinite(c3)) = inf;
    c4(~isfinite(c4)) = inf;

    C = [c1 c2 c3 c4];                 % Nx4
    [best_cost, which_branch] = min(C, [], 2, "omitnan");

    no_valid = ~isfinite(best_cost) | (best_cost==inf);
    which_branch(no_valid) = 1;     % dummy to avoid bad indexing
    best_cost(no_valid) = inf;      % keep cost invalid

    % Build a "best" struct that is vector-valued, per-rp
    best = struct();
    best.cost = best_cost;
    best.best_angle = a1;              % placeholder, overwritten below
    best.sign  = nan(N,1);
    best.angle = nan(N,1);

    % sign/angle per branch
    sign_map  = [ 1  1 -1 -1];
    angle_map = [ 1  2  1  2];

    best.sign  = sign_map(which_branch).';
    best.angle = angle_map(which_branch).';

    % pick the per-branch best_angle vector
    best.best_angle = a1;
    best.best_angle(which_branch==2) = a2(which_branch==2);
    best.best_angle(which_branch==3) = a3(which_branch==3);
    best.best_angle(which_branch==4) = a4(which_branch==4);

    if ~isempty(all_results.result1)
        best.rp = all_results.result1.rp;
    elseif ~isempty(all_results.result2)
        best.rp = all_results.result2.rp;
    elseif ~isempty(all_results.result3)
        best.rp = all_results.result3.rp;
    else
        best.rp = all_results.result4.rp;
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

function [cost, best_angle] = calc_cost_vec(result, S, angle_turn)
% result.tof_tot_1, result.theta_b1_1, result.tof_tot_2, result.theta_b1_2
% can be scalars or Nx1 vectors.
% Returns:
%   cost       Nx1 vector
%   best_angle Nx1 vector (1 or 2)

    calc = @(tof_tot, theta_b1) ...
        (abs(tof_tot - S)./S).^2 + ...
        (abs(theta_b1 - angle_turn)./angle_turn).^2;

    cost1 = calc(result.tof_tot_1, result.theta_b1_1);
    cost2 = calc(result.tof_tot_2, result.theta_b1_2);

    % element-wise min + argmin (matches your scalar logic)
    cost = cost1;
    best_angle = ones(size(cost1));

    mask = cost2 < cost1;
    cost(mask) = cost2(mask);
    best_angle(mask) = 2;

    bad1 = ~isfinite(result.tof_tot_1) | ~isfinite(result.theta_b1_1);
    bad2 = ~isfinite(result.tof_tot_2) | ~isfinite(result.theta_b1_2);
    
    both_bad = bad1 & bad2;
    cost(both_bad) = inf;
    best_angle(both_bad) = 1;
end

function result = calc_result(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, sign, angle)

result = [];

[valid, theta_b2, e_out, p_out, w_out, tof_12] = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, sign, angle);

if ~any(valid), return; end

e_out(~valid) = nan;
p_out(~valid) = nan;
w_out(~valid) = nan;

theta_b1_1 = calc_intersect_vec(p_out, norm(r0), e_out, w_out, 1);

if all(~isfinite(theta_b1_1))
    return
end

tof_21_1 = calc_tof_vec(theta_b2, theta_b1_1, w_out, e_out, p_out, mu);
tof_tot_1 = tof_12 + tof_21_1;

theta_b1_2 = calc_intersect_vec(p_out, norm(r0), e_out, w_out, 2);
tof_21_2 = calc_tof_vec(theta_b2, theta_b1_2, w_out, e_out, p_out, mu);
tof_tot_2 = tof_12 + tof_21_2;

bad1 = ~isfinite(tof_21_1) | ~isfinite(theta_b1_1);
bad2 = ~isfinite(tof_21_2) | ~isfinite(theta_b1_2);

tof_tot_1(bad1) = nan;
tof_tot_2(bad2) = nan;

result = struct( ...
    'rp',         rp, ...
    'theta_b2',   theta_b2, ...
    'theta_b1_1', theta_b1_1, ...
    'tof_tot_1',  tof_tot_1, ...
    'theta_b1_2', theta_b1_2, ...
    'tof_tot_2',  tof_tot_2);
end

function [valid, theta_intersect, e_out, p_out, w_out, tof_12] = calc_grav_assist(r0_vec, v0_vec, mu, mu_b2, r_b2, h_b2, e_b2, rp_vec, sign, angle)

    rp_vec = rp_vec(:);
    N = numel(rp_vec);

    e_out = nan(N,1);
    p_out = nan(N,1);
    w_out = nan(N,1);

    % inbound (scalar)
    [~, h_sc, e_sc, p_sc, w_sc, ~] = calc_orbit(r0_vec, v0_vec, mu);
    if isempty(e_sc) || e_sc >= 1
        valid = false(N,1);
        theta_intersect = [];
        tof_12 = nan(N,1);
        return;
    end
    a_sc = p_sc / (1 - e_sc^2);

    % intersection with body2 (scalar)
    theta_intersect = calc_intersect(p_sc, r_b2, e_sc, w_sc, angle);
    if isempty(theta_intersect)
        valid = false(N,1);
        tof_12 = nan(N,1);
        return;
    end

    nu_int = theta_intersect - w_sc;

    % SC vel at intersect (scalar)
    vr_sc = mu/h_sc * e_sc * sin(nu_int);
    vt_sc = mu/h_sc * (1 + e_sc * cos(nu_int));

    % body2 vel at intersect (scalar)
    vr_b2 = mu/h_b2 * e_b2 * sin(theta_intersect);
    vt_b2 = mu/h_b2 * (1 + e_b2 * cos(theta_intersect));

    % v-infinity (scalar 2x1)
    v_inf_in = [vr_sc - vr_b2; vt_sc - vt_b2];
    v_inf    = hypot(v_inf_in(1), v_inf_in(2));

    % tof Earth->body2 (scalar replicated)
    E  = 2 * atan( sqrt((1 - e_sc) / (1 + e_sc)) * tan(nu_int / 2) );
    M  = E - e_sc*sin(E);
    tof_12 = sqrt(a_sc^3 / mu) * M;
    tof_12 = tof_12 + zeros(N,1);

    % ----- vectorized turning over rp -----
    e_hyp = 1 + rp_vec * (v_inf^2) / mu_b2;
    valid = e_hyp > 1;

    delta = sign * 2 * asin(1 ./ e_hyp);

    c = cos(delta); s = sin(delta);
    vix = v_inf_in(1); viy = v_inf_in(2);

    vxo = c*vix - s*viy;   % Nx1
    vyo = s*vix + c*viy;   % Nx1

    vr_out = vr_b2 + vxo;
    vt_out = vt_b2 + vyo;

    ct = cos(theta_intersect);
    st = sin(theta_intersect);

    rx = r_b2 * ct;    % scalar
    ry = r_b2 * st;    % scalar

    vx = vr_out*ct - vt_out*st;   % Nx1
    vy = vr_out*st + vt_out*ct;   % Nx1

    N = numel(vx);   % number of rp samples

    r_mat = [ ...
        rx * ones(N,1), ...
        ry * ones(N,1), ...
        zeros(N,1) ];
    
    v_mat = [ ...
        vx(:), ...
        vy(:), ...
        zeros(N,1) ];
    
    [r_out, h_out, e_out, p_out, w_out, nu_out] = ...
        calc_orbit_vec(r_mat, v_mat, mu);

    % extra constraints
    valid = valid & isfinite(e_out) & (e_out < 1);
    rp_helio = p_out ./ (1 + e_out);
    valid = valid & (rp_helio <= norm(r0_vec));

    e_out(~valid) = nan;
    p_out(~valid) = nan;
    w_out(~valid) = nan;
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

function theta = calc_intersect_vec(p, r_const, e, w, angle)
    theta = nan(size(p));
    c = (p./r_const - 1)./e;

    ok = isfinite(c) & isfinite(e) & isfinite(w) & (abs(c) <= 1);
    if ~any(ok), return; end

    nu = acos(c(ok));
    phi_a = mod(w(ok) + nu, 2*pi);
    phi_b = mod(w(ok) - nu, 2*pi);

    if angle == 1
        theta(ok) = min(phi_a, phi_b);
    else
        theta(ok) = max(phi_a, phi_b);
    end
end

function [r, h, e, p, w, nu] = calc_orbit(r_vec, v_vec, mu)

h_vec = cross(r_vec, v_vec);
e_vec = cross(v_vec, h_vec) / mu - r_vec/norm(r_vec);

h = norm(h_vec);
e = norm(e_vec);

p = h^2/mu;

w = atan2(e_vec(2), e_vec(1));

theta = atan2(r_vec(2), r_vec(1));
nu = theta - w;

r = p ./ (1 + e*cos(nu));
end

function [r, h, e, p, w, nu] = calc_orbit_vec(r_vec, v_vec, mu)
% r_vec: Nx3 (or Nx2)  [x y (0)]
% v_vec: Nx3 (or Nx2)  [vx vy (0)]
% Outputs: Nx1 vectors matching calc_orbit() semantics

    % --- Accept Nx2 too ---
    if size(r_vec,2) == 2
        r_vec = [r_vec, zeros(size(r_vec,1),1)];
    end
    if size(v_vec,2) == 2
        v_vec = [v_vec, zeros(size(v_vec,1),1)];
    end

    % --- Geometry (row-wise) ---
    h_vec = cross(r_vec, v_vec, 2);                         % Nx3
    r_norm = sqrt(sum(r_vec.^2, 2));                        % Nx1

    e_vec = cross(v_vec, h_vec, 2) ./ mu - r_vec ./ r_norm; % Nx3

    % --- Scalars per row ---
    h = sqrt(sum(h_vec.^2, 2));                             % Nx1
    e = sqrt(sum(e_vec.^2, 2));                             % Nx1

    p = (h.^2) ./ mu;                                       % Nx1

    w = atan2(e_vec(:,2), e_vec(:,1));                      % Nx1

    theta = atan2(r_vec(:,2), r_vec(:,1));                  % Nx1
    nu = theta - w;                                         % Nx1

    r = p ./ (1 + e .* cos(nu));                            % Nx1
end

function [S, angle_turn] = calc_synodic(mu, r_b1, r_b2)
    T1 = 2*pi * sqrt(r_b1 ^ 3 / mu);
    T2 = 2*pi * sqrt(r_b2 ^ 3 / mu);
    
    T_ratio = T2 / T1;
    
    S = 1 / abs(1/T1 - 1/T2);
    
    angle_turn = mod(T_ratio * 2*pi, 2*pi);
    if angle_turn > pi
        angle_turn = angle_turn - 2*pi;
    end
    
    angle_turn = -angle_turn;
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

function tof = calc_tof_vec(theta_b2, theta_b1, w_out, e_out, p_out, mu)
    tof = nan(size(theta_b1));

    ok = isfinite(theta_b1) & isfinite(e_out) & (e_out < 1) & isfinite(p_out) & isfinite(w_out);
    if ~any(ok), return; end

    nu2 = theta_b2 - w_out(ok);
    nu1 = theta_b1(ok) - w_out(ok);

    a = p_out(ok) ./ (1 - e_out(ok).^2);

    E2 = 2 * atan( sqrt((1 - e_out(ok)) ./ (1 + e_out(ok))) .* tan(nu2/2) );
    E1 = 2 * atan( sqrt((1 - e_out(ok)) ./ (1 + e_out(ok))) .* tan(nu1/2) );

    M2 = E2 - e_out(ok).*sin(E2);
    M1 = E1 - e_out(ok).*sin(E1);

    dM = M1 - M2;
    dM(dM < 0) = dM(dM < 0) + 2*pi;

    tof(ok) = sqrt(a.^3 / mu) .* dM;
end