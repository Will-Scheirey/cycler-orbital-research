clear; clc; close all

mu = 132712000000;
AU = 1.496e+8;

r_b1 = AU;
v_b1 = sqrt(mu / r_b1);

% === Body 2 (set one block only) ===
% Mars
% mu_b2 = 42828;     r_b2 = 2.2737e8;  radius_b2 = 3396;

% Jupiter
mu_b2 = 126686000;  r_b2 = 778.6e6;   radius_b2 = 80000;

soi_b2 = r_b2 * (mu_b2 / mu)^(2/5);

min_rp = radius_b2 + 100;
max_rp = soi_b2;

dv_min = sqrt(2*mu / r_b1 - 2*mu / (r_b1 + r_b2)) - sqrt(mu/r_b1);
dv_max = dv_min + 10;

v_b2 = sqrt(mu / r_b2);
h_b2 = r_b2 * v_b2;
e_b2 = 0;

[S, angle_turn] = calc_synodic(mu, r_b1, r_b2);
n = 1;
S = S * n;
angle_turn = mod(angle_turn * n, 2*pi);

num_points = 100;

t_all  = linspace(-pi, pi, num_points);
v_all  = linspace(dv_min, dv_max, num_points);
rp_all = linspace(min_rp, max_rp, num_points);

% We now minimize over rp inside the solver, so costs/results are 2D
all_results = cell(num_points, num_points);
costs       = nan(num_points, num_points);

r0 = [r_b1; 0; 0];

% Optional: precompute sin/cos for speed
sin_t = sin(t_all);
cos_t = cos(t_all);

for t_idx = 1:num_points
    t = t_all(t_idx);

    for v_idx = 1:num_points
        v_mag = v_all(v_idx);

        v0 = [v_mag * sin_t(t_idx);
              v_b1  + v_mag * cos_t(t_idx);
              0];

        % quick feasibility gate (intersection exists)
        r_out = calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, 1e8, 1, 1);
        if isempty(r_out)
            continue;
        end

        [best_res, best_cost] = solve_best_over_rp( ...
            r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp_all, S, angle_turn, t, v_mag);

        costs(t_idx, v_idx)       = best_cost;
        all_results{t_idx, v_idx} = best_res;
    end
end

%% PLOT (2D costs now)
figure(1); clf
surf(v_all, t_all, costs, 'LineStyle','none')
xlabel("V Mag")
ylabel("t")
zlabel("Cost")

return

%% ===== Optional debug plotting (same as before, but iterate 2D results) =====
costs_temp   = costs(:);
results_temp = all_results(:);

min_cost = min(costs_temp, [], "omitmissing");
min_result = results_temp{find(costs_temp == min_cost, 1, 'first')};

for z = 1:numel(results_temp)
    min_result = results_temp{z};
    if isempty(min_result), continue; end
    if ~isfield(min_result,'cost') || ~isfinite(min_result.cost), continue; end

    t     = min_result.t;
    v_mag = min_result.v0;
    rp    = min_result.rp;

    v0 = [v_mag * sin(t);
          v_b1  + v_mag * cos(t);
          0];

    [~, ~, e_sc, p_sc, w_sc, ~] = calc_orbit(r0, v0, mu);

    [r_out, ~, e_out, p_out, w_out, ~, theta_intersect, ~] = ...
        calc_grav_assist(r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp, min_result.sign, min_result.angle);

    if isempty(r_out), continue; end

    theta_sc = linspace(0, theta_intersect, 1000);
    nu_sc = theta_sc - w_sc;
    r1 = p_sc ./ (1 + e_sc*cos(nu_sc));

    theta_plot = linspace(0,2*pi,1000);

    theta_b1 = calc_intersect(p_out, AU, e_out, w_out, min_result.best_angle);
    theta_out = linspace(theta_intersect, theta_b1, 1000);
    nu_out = theta_out - w_out;
    r2 = p_out ./ (1 + e_out*cos(nu_out));

    figure(2); clf
    plot(r2.*cos(theta_out), r2.*sin(theta_out), '-g', 'LineWidth', 1.5); hold on
    plot(r1.*cos(theta_sc),  r1.*sin(theta_sc),  '--k', 'LineWidth', 1.5); hold on

    plot(r_b1*cos(theta_plot), r_b1*sin(theta_plot), '-b'); hold on
    plot(r_b2*cos(theta_plot), r_b2*sin(theta_plot), '-r')

    plot(r_b2*cos(theta_intersect), r_b2*sin(theta_intersect), 'mx')
    plot(AU*cos(theta_b1), AU*sin(theta_b1), 'mx')

    axis equal
    xlim([-1,1]*r_b2 * 4)
    ylim([-1,1]*r_b2 * 4)
    title(sprintf("Cost = %0.5f, Rp = %0.1f", min_result.cost, min_result.rp))
    drawnow;
    pause(1)
end

%% =======================================================================
%% ====================== VECTORIZED RP SOLVER ============================
%% =======================================================================

function [best_res, best_cost] = solve_best_over_rp( ...
    r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp_vec, S, angle_turn, t, v_mag)

    best_cost = inf;
    best_res  = [];

    signs  = [ 1,  1, -1, -1];
    angles = [ 1,  2,  1,  2];

    for k = 1:4
        [res_k, cost_k] = solve_branch_over_rp( ...
            r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp_vec, signs(k), angles(k), S, angle_turn, t, v_mag);

        if isempty(res_k), continue; end
        if cost_k < best_cost
            best_cost = cost_k;
            best_res  = res_k;
        end
    end
end

function [best_res, best_cost] = solve_branch_over_rp( ...
    r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp_vec, sign, angle, S, angle_turn, t, v_mag)

    best_res  = [];
    best_cost = inf;

    rp_vec = rp_vec(:);

    [valid, theta_b2, e_out, p_out, w_out, tof_12] = calc_grav_assist_rp_batch( ...
        r0, v0, mu, mu_b2, r_b2, h_b2, e_b2, rp_vec, sign, angle);

    if ~any(valid), return; end

    r_b1 = norm(r0);

    theta_b1_1 = calc_intersect_vec(p_out, r_b1, e_out, w_out, 1);
    theta_b1_2 = calc_intersect_vec(p_out, r_b1, e_out, w_out, 2);

    tof_21_1 = calc_tof_vec(theta_b2, theta_b1_1, w_out, e_out, p_out, mu);
    tof_21_2 = calc_tof_vec(theta_b2, theta_b1_2, w_out, e_out, p_out, mu);

    tof_tot_1 = tof_12 + tof_21_1;
    tof_tot_2 = tof_12 + tof_21_2;

    c1 = calc_cost_vec(tof_tot_1, theta_b1_1, S, angle_turn);
    c2 = calc_cost_vec(tof_tot_2, theta_b1_2, S, angle_turn);

    [cmin, which] = min([c1.'; c2.'], [], 1);
    cmin  = cmin(:);
    which = which(:);

    [best_cost, idx] = min(cmin);

    if ~isfinite(best_cost), best_res = []; best_cost = inf; return; end

    best_res = struct( ...
        'rp',         rp_vec(idx), ...
        'theta_b2',   theta_b2, ...
        'theta_b1_1', theta_b1_1(idx), ...
        'tof_tot_1',  tof_tot_1(idx), ...
        'theta_b1_2', theta_b1_2(idx), ...
        'tof_tot_2',  tof_tot_2(idx), ...
        'sign',       sign, ...
        'angle',      angle, ...
        'best_angle', which(idx), ...
        'cost',       best_cost, ...
        't',          t, ...
        'v0',         v_mag );
end

function [valid, theta_intersect, e_out, p_out, w_out, tof_12] = calc_grav_assist_rp_batch( ...
    r0_vec, v0_vec, mu, mu_b2, r_b2, h_b2, e_b2, rp_vec, sign, angle)

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

    [h_out, e_out, p_out, w_out] = calc_orbit_2d_batch(rx, ry, vx, vy, mu);

    % extra constraints
    valid = valid & isfinite(e_out) & (e_out < 1);
    rp_helio = p_out ./ (1 + e_out);
    valid = valid & (rp_helio <= norm(r0_vec));

    e_out(~valid) = nan;
    p_out(~valid) = nan;
    w_out(~valid) = nan;
end

function [h, e, p, w] = calc_orbit_2d_batch(rx, ry, vx, vy, mu)
    h = rx.*vy - ry.*vx;                % Nx1
    r = hypot(rx, ry);                  % scalar

    ex = (vy.*h)/mu - rx/r;
    ey = (-vx.*h)/mu - ry/r;

    e = hypot(ex, ey);
    p = (h.^2)/mu;
    w = atan2(ey, ex);
end

function theta = calc_intersect_vec(p, r_const, e, w, angle)
    theta = nan(size(p));
    c = (p./r_const - 1)./e;

    ok = isfinite(c) & (abs(c) <= 1);
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

function c = calc_cost_vec(tof_tot, theta_b1, S, angle_turn)
    c = inf(size(tof_tot));
    ok = isfinite(tof_tot) & isfinite(theta_b1);
    if ~any(ok), return; end
    c(ok) = (abs(tof_tot(ok) - S)./S).^2 + (abs(theta_b1(ok) - angle_turn)./angle_turn).^2;
end

%% =======================================================================
%% ====================== ORIGINAL SUPPORT FUNCTIONS ======================
%% =======================================================================

function [r_out, h_out, e_out, p_out, w_out, nu_out, theta_intersect, tof1] = ...
    calc_grav_assist(r0_vec, v0_vec, mu, mu_b2, r_b2, h_b2, e_b2, rp, sign, angle)

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
    if isempty(theta_intersect), return; end

    nu_int = theta_intersect - w_sc;

    vr_sc = mu/h_sc * e_sc * sin(nu_int);
    vt_sc = mu/h_sc * (1 + e_sc * cos(nu_int));

    vr_b2 = mu/h_b2 * e_b2 * sin(theta_intersect);
    vt_b2 = mu/h_b2 * (1 + e_b2 * cos(theta_intersect));

    v_inf_in = [vr_sc - vr_b2; vt_sc - vt_b2];
    v_inf    = norm(v_inf_in);

    e  = 1 + rp * v_inf^2 / mu_b2;
    delta = sign * 2 * asin(1 ./ e);

    R = [cos(delta) -sin(delta); sin(delta) cos(delta)];
    v_inf_out = R * v_inf_in;

    vr_out = vr_b2 + v_inf_out(1);
    vt_out = vt_b2 + v_inf_out(2);

    rx_out = r_b2*cos(theta_intersect);
    ry_out = r_b2*sin(theta_intersect);

    vx_out = vr_out*cos(theta_intersect) - vt_out*sin(theta_intersect);
    vy_out = vr_out*sin(theta_intersect) + vt_out*cos(theta_intersect);

    r_vec_out = [rx_out, ry_out, 0];
    v_vec_out = [vx_out, vy_out, 0];

    [r_out, h_out, e_out, p_out, w_out, nu_out] = calc_orbit(r_vec_out, v_vec_out, mu);

    E = 2 * atan(sqrt((1 - e_sc) / (1 + e_sc)) * tan(nu_int / 2));
    Me = E - e_sc * sin(E);
    tof1 = sqrt(a_sc^3 / mu) * Me;
end

function theta_intersect = calc_intersect(p_sc, r_b2, e_sc, w_sc, angle)
    theta_intersect = [];
    c = (p_sc/r_b2 - 1)/e_sc;
    if abs(c) > 1, return; end
    nu_int = acos(c);
    if imag(nu_int) ~= 0, return; end

    phi_a = mod(w_sc + nu_int, 2*pi);
    phi_b = mod(w_sc - nu_int, 2*pi);

    if angle == 1
        theta_intersect = min(phi_a, phi_b);
    else
        theta_intersect = max(phi_a, phi_b);
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