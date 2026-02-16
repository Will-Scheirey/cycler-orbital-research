clear; clc; close all
mu = 132712000000;

AU = 1.496e+8;

r_b1 = AU; % Earth
v_b1 = sqrt(mu / r_b1); % Velocity of Earth
mu_b1 = 398600;
rp_min = 6378;

r_b2 = 2.2737e8; % Mars
v_b2 = sqrt(mu / r_b2); % Velocity of Mars

S = calc_synodic(mu, r_b1, r_b2);

AR_min = 0.9;
TR_min = 0.85;
p_max = 4;

h_max = 5*p_max;
s_max = 3*p_max;

r0 = [r_b1; 0; 0];

n = sqrt(mu / r_b1^3);

solutions = cell(1,1);

all_sol_idx = 1;

v0 = v_b1 * [0; 1; 0];
%{
for p = 1:p_max
    for h = 0:h_max
        for s = 1:s_max
            
            [N_max, tof, sol_all] = run_algorithm(S, p, h, s, n, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min);

            if isempty(sol_all)
                continue
            end

            data = struct('p', p, 'h', h, 's', s, 'N_max', N_max, 'tof', tof, 'solution_data', sol_all);
            solutions{all_sol_idx} = data;
            all_sol_idx = all_sol_idx + 1;
        end
    end
end
%}
p = 2;
h = 3;
s = 1;

[N_max, tof, sol_all] = run_algorithm(S, p, h, s, n, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min);

data = struct('p', p, 'h', h, 's', s, 'N_max', N_max, 'tof', tof, 'solution_data', sol_all);

feasible_solutions = get_feasible_solutions({data}, true);

print_solutions(feasible_solutions)

return

%% Processing

feasible_solutions = get_feasible_solutions(all_sol_idx, solutions);
print_solutions(feasible_solutions)

%% Functions

function feasible_solutions = get_feasible_solutions(solutions, allow_all)

if nargin < 3
    allow_all = true;
end

feasible_solutions = {};
feasible_idx = 1;

for i = 1:length(solutions)

    solution = solutions{i};
    solution_data = solution.solution_data;

    solution_names = fieldnames(solution_data);

    for n = 1:length(solution_names)
        this_solution = solution_data.(solution_names{n});

        fast_solution = this_solution.fast;

        if isfield(fast_solution, "feasible")
            if fast_solution.feasible(1) || allow_all
                fast_solution.p = solution.p;
                fast_solution.h = solution.h;
                fast_solution.s = solution.s;
                fast_solution.direction = '-';
                feasible_solutions{feasible_idx} = fast_solution;
                feasible_idx = feasible_idx + 1;
            end
        end

        slow_solution = this_solution.slow;

        if isfield(slow_solution, "feasible")
            if slow_solution.feasible(1) || allow_all
                slow_solution.p = solution.p;
                slow_solution.h = solution.h;
                slow_solution.s = solution.s;
                slow_solution.direction = '+';
                feasible_solutions{feasible_idx} = slow_solution;
                feasible_idx = feasible_idx + 1;
            end
        end
    end
end
end

function print_solutions(feasible_solutions)
    num_feasible = length(feasible_solutions);
    
    for i = 1:num_feasible
        sol = feasible_solutions{i};
        fprintf("Feasible solution for %d.%d.%d.%s%d", sol.p, sol.h, sol.s, sol.direction, sol.i);
        fprintf("\t AR: %0.2f; TR: %0.2f; Earth V_inf: %0.1f\n", sol.AR, sol.TR, sol.v_inf);
    end
end

function [N_max, tof, sol_all] = run_algorithm(S, p, h, s, n, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min)

sol_all = [];

tof = identical_s_tof(S, p, h, s);
if tof < 0
    return;
end

theta = n * tof;

r1 = r_b1 * [cos(theta); sin(theta); 0];
v1 = v_b1 * [-sin(theta); cos(theta); 0];

[N_max, a_all, vd_all, va_all] = calc_multirev_lambert(r0, r1, mu, tof);

if N_max == 0
    return
end

sol_all = struct();
sol_idx = 1;

for i = 1 : N_max+1

    fast_struct = struct();

    a_fast = a_all(i, 1);
    vd_fast = vd_all{i, 1};
    va_fast = va_all{i, 1};

    [v_inf_minus_fast, deltas_fast] = calc_sequence(vd_fast, va_fast, h, s, v0, v1);
    if ~isempty(v_inf_minus_fast)
        ra_fast = calc_ra(r0, vd_fast, mu);
        [AR_fast, TR_fast, feasible_fast, max_delta_fast] = is_feasible(a_fast, rp_min, mu_b1, r_b2, norm(v_inf_minus_fast{1}), deltas_fast, TR_min, AR_min, ra_fast);
        fast_struct = struct( ...
            'v_inf', norm(v_inf_minus_fast{1}), ...
            'AR', AR_fast, ...
            'TR', TR_fast, ...
            'feasible', feasible_fast, ...
            'max_delta', max_delta_fast, ...
            'i', i-1);
        fast_struct.v_inf_minus = v_inf_minus_fast;
        fast_struct.deltas = deltas_fast;
    end

    slow_struct = struct();

    a_slow = a_all(i, 2);
    vd_slow = vd_all{i, 2};
    va_slow = va_all{i, 2};

    [v_inf_minus_slow, deltas_slow] = calc_sequence(vd_slow, va_slow, h, s, v0, v1);
    if ~isempty(v_inf_minus_slow)
        ra_slow = calc_ra(r0, vd_slow, mu);
        [AR_slow, TR_slow, feasible_slow, max_delta_slow] = is_feasible(a_slow, rp_min, mu_b1, r_b2, norm(v_inf_minus_slow{1}), deltas_slow, TR_min, AR_min, ra_slow);
        slow_struct = struct( ...
            'v_inf', norm(v_inf_minus_slow{1}), ...
            'AR', AR_slow, ...
            'TR', TR_slow, ...
            'feasible', feasible_slow, ...
            'max_delta', max_delta_slow, ...
            'i', i-1);
        slow_struct.v_inf_minus = v_inf_minus_slow;
        slow_struct.deltas = deltas_slow;
    end

    sol_all.(sprintf('sol_%d',sol_idx)) = struct( ...
        'fast', fast_struct, ...
        'slow', slow_struct...
        );

    sol_idx = sol_idx + 1;
end

end

function [N_max, a_all, vd_all, va_all] = calc_multirev_lambert(r1_vec, r2_vec, mu, tof)

r1 = norm(r1_vec);
r2 = norm(r2_vec);

theta = acos(dot(r1_vec, r2_vec) / (r1*r2));

N_max = calc_N_max(r1, r2, theta, tof, mu);

options = optimset();

a_all = nan(N_max, 2);
vd_all = cell(N_max, 2);
va_all = cell(N_max, 2);

for N = 0:N_max
    [~, ~, s, c] = get_lambert_angles(r1, r2, N, theta, true);
    a_m = s/2;

    [tof_eqn_fast, alpha_fast, beta_fast] = lambert_lagrange_eqn(r1, r2, theta, mu, N, true);
    [tof_eqn_slow, alpha_slow, beta_slow] = lambert_lagrange_eqn(r1, r2, theta, mu, N, false);

    tof_target = tof;

    eqn_fast = @(a) tof_eqn_fast(a) - tof_target;
    eqn_slow = @(a) tof_eqn_slow(a) - tof_target;

    a_min = a_m*(1+1e-3);
    a_max = 1e6*a_m;

    a_fast = NaN; a_slow = NaN;

    if sign(eqn_fast(a_min)) ~= sign(eqn_fast(a_max))
        a_fast = fzero(eqn_fast, [a_min, a_max], options);
    end

    if sign(eqn_slow(a_min)) ~= sign(eqn_slow(a_max))
        a_slow = fzero(eqn_slow, [a_min, a_max], options);
    end

    [vd_fast, va_fast] = calc_v_inf(mu, a_fast, c, alpha_fast(a_fast), beta_fast(a_fast), r1_vec, r2_vec);
    [vd_slow, va_slow] = calc_v_inf(mu, a_slow, c, alpha_slow(a_slow), beta_slow(a_slow), r1_vec, r2_vec);

    a_all (N+1, 1) = a_fast;
    vd_all{N+1, 1} = vd_fast;
    va_all{N+1, 1} = va_fast;

    a_all (N+1, 2) = a_slow;
    vd_all{N+1, 2} = vd_slow;
    va_all{N+1, 2} = va_slow;
end

end

function N_max = calc_N_max(r1, r2, theta, T_max, mu)

N_max_search = 10;
N_max = 0;

options = optimset('Display', 'off');

for N = 0:N_max_search
    [alpha, beta, s] = get_lambert_angles(r1, r2, N, theta, true);
    a_m = s/2;

    tof_eqn_fast = lambert_lagrange_eqn(r1, r2, theta, mu, N, true);


    if N > 0
        f = @(a) (6*N*pi + 3*(alpha(a) - beta(a)) - (sin(alpha(a))) - sin(beta(a))) .* ...
            (sin(alpha(a) - beta(a)) + (sin(alpha(a)) - sin(beta(a)))) - 8*(1 - cos(alpha(a) - beta(a)));

        zero_a = fzero(f, [1, 10]*a_m, options);

        tof_min = tof_eqn_fast(zero_a);

        if tof_min < T_max
            N_max = N;
        else
            break
        end
    else

    end
end

end

function [tof_eqn, alpha, beta] = lambert_lagrange_eqn(r1, r2, theta, mu, N, fast)
[alpha, beta] = get_lambert_angles(r1, r2, N, theta, fast);

tof_eqn = @(a) sqrt(a.^3 / mu) .* (2*N*pi + alpha(a) - beta(a) - sin(alpha(a)) + sin(beta(a)));
end

function [alpha, beta, S, c] = get_lambert_angles(r1, r2, N, theta, fast)
c  = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta));
S  = 0.5*(c + r1 + r2);

alpha0 = @(a) 2 * asin(sqrt(S ./ (2*a)));
beta0  = @(a) 2 * asin(sqrt((S - c) ./ (2*a)));

if fast
    alpha = @(a) alpha0(a);
else
    alpha = @(a) 2*pi - alpha0(a);
end

if (theta - 2*pi*N) < pi
    beta = @(a) beta0(a);
else
    beta = @(a) -beta0(a);
end
end

function [v_D, v_A] = calc_v_inf(mu, a, c, alpha, beta, r1_vec, r2_vec)

A = sqrt(mu / (4*a)) * cot(alpha/2);
B = sqrt(mu / (4*a)) * cot(beta/2);

r1_hat = r1_vec / norm(r1_vec);
r2_hat = r2_vec / norm(r2_vec);
c_hat  = (r2_vec - r1_vec) / c;

v_D = (B + A)*c_hat + (B - A)*r1_hat;
v_A = (B + A)*c_hat - (B - A)*r2_hat;

end

%{
function [phi_fr, phi_gr] = calc_phi(vd_vec, va_vec, ve_vec_dep, ve_vec_arr)
ve = norm(ve_vec_dep);
v_inf_vec   = vd_vec - ve_vec_dep;
v_inf_1_vec = va_vec - ve_vec_arr;

v_inf   = norm(v_inf_vec);
v_inf_1 = norm(v_inf_1_vec);

phi_fr = -asin(v_inf / (2*ve));

%{
phi_gr = pi/2 - acos((v_inf_1_vec' * ve_vec) / (v_inf * ve));
%}

u = (v_inf_1_vec' * ve_vec_arr) / (v_inf_1 * ve);
u = max(-1,min(1,u));
phi_gr = asin(u);

end
%}
function [phi_fr, phi_gr] = calc_phi(vd_vec, va_vec, ve_vec_dep, ve_vec_arr)
ve_dep = norm(ve_vec_dep);
ve_arr = norm(ve_vec_arr);

v_inf_vec   = vd_vec - ve_vec_dep;   % dep v-infinity
v_inf_1_vec = va_vec - ve_vec_arr;   % arr v-infinity

v_inf   = norm(v_inf_vec);
v_inf_1 = norm(v_inf_1_vec);

q = v_inf / (2*ve_dep);
if ~isfinite(q) || q > 1
    phi_fr = NaN; phi_gr = NaN;
    return
end
phi_fr = -asin(q);

u = (v_inf_1_vec' * ve_vec_arr) / (v_inf_1 * ve_arr);
u = max(-1,min(1,u));
phi_gr = asin(u);
end

function [v_inf_minus, deltas, v_inf_local] = calc_sequence(vd, va, h, s, v_e_dep, v_e_arr)

v_inf_minus = [];
v_inf_local = [];
deltas = [];

[phi_fr, phi_gr] = calc_phi(vd, va, v_e_dep, v_e_arr);

if ~isfinite(phi_fr) || ~isfinite(phi_gr)
    return
end

delta_c = pi - 2 * abs(phi_gr);

hj = floor(h/s);

if mod(hj, 2) == 0
    fj = hj/2 + 1;
else
    fj = 2 * floor(hj/4 + 1);
end

lambda_a = pi / (fj - 2);

if fj == 1
    delta_minimax = delta_c;
elseif fj == 2
    delta_minimax = acos(sin(phi_gr) * sin(phi_fr));
elseif fj > 2
    delta_min = acos(cos(phi_fr)*cos(phi_gr) + sin(phi_fr)*sin(phi_gr));
    delta_a   = acos(cos(phi_fr)^2 * cos(lambda_a) + sin(phi_fr)^2);

    if delta_min >= delta_a
        delta_minimax = delta_min;
    else
        lb = @(l) (pi - 2*l) / (fj - 2);
        fcn = @(l) cos(phi_fr)^2*cos(l)*cos(lb(l) + l) - cos(phi_fr)*cos(l)*cos(phi_gr) + ...
            cos(phi_fr)^2*sin(l)*sin(lb(l) + l) + sin(phi_fr)^2 - sin(phi_fr)*sin(phi_gr);
        l = fsolve(fcn, 0, optimset('Display', 'off'));

        delta_minimax = acos(cos(phi_gr)*cos(phi_fr)*cos(l) + sin(phi_gr)*sin(phi_fr));
    end
end

h_all = zeros(s, 1);

if delta_c > delta_minimax
    h_all(1:s-1) = floor(h/s);
    h_all(s)   = floor(h/s) + mod(h, s);
else
    h_all(1) = h;
    h_all(2:s) = 0;
end

% Assemble reference frame

v_inf_1_minus = va - v_e_arr;
v_inf_mag = norm(v_inf_1_minus);

if v_inf_mag == 0
    return
end

z_hat = v_e_arr / norm(v_e_arr);

x_vec = v_inf_1_minus - dot(v_inf_1_minus, z_hat)*z_hat;
x_hat = x_vec / norm(x_vec);

y_hat = cross(z_hat, x_hat);
y_hat = y_hat / norm(y_hat);

C = [x_hat, y_hat, z_hat];          % local->global

if fj == 1
    deltas = delta_c;

    [vx1, vy1, vz1] = sph2cart(0, phi_gr, v_inf_mag);
    v_inf_minus = {[vx1; vy1; vz1]};
    v_global = v_inf_minus{1};
    % v_global = C * [vx1;vy1;vz1];

    v_inf_local = {v_global, v_global};
    v_inf_local{2}(1) = -v_inf_local{2}(1);
    return
end

if fj == 2
    [vx1, vy1, vz1] = sph2cart(-0,      phi_gr, v_inf_mag);
    [vx2, vy2, vz2] = sph2cart(-pi/2,   phi_gr, v_inf_mag);

    v_inf_minus = cell(2,1);
    v_inf_minus{1} = C * [vx1; vy1; vz1];
    v_inf_minus{2} = C * [vx2; vy2; vz2];

    deltas = delta_minimax;
    return
end

v_inf_local = cell(fj+1, 1);

v_inf_minus = cell(fj, 1);
deltas = zeros(fj - 1, 1);

for k=1:fj
    if k == 1
        lat = phi_gr;
        lon = 0;
    else
        lat = phi_fr;
        if delta_min > delta_a
            lon = lambda_a * (k - 2);
        else
            lon = l + lb(l) * (k-2);
        end
    end

    [vx, vy, vz] = sph2cart(-lon, lat, v_inf_mag);

    v_global = C * [vx;vy;vz];

    v_inf_local{k} = [vx; vy; vz];
    v_inf_minus{k} = v_global;

    if k >= 2
        v_prev = v_inf_minus{k-1};
        v_curr = v_inf_minus{k};

        deltas(k-1) = acos(dot(v_prev, v_curr) / (norm(v_prev)*norm(v_curr)));
    end
end

v_inf_local{fj+1} = v_inf_local{1};
v_inf_local{fj+1}(1) = -v_inf_local{fj+1}(1);

end

function [AR, TR, feasible, delta_max] = is_feasible(a, rp_min, mu, r_b2, v_inf_mag, deltas, TR_min, AR_min, ra)

AR = ra / r_b2;

e = 1 + (rp_min * v_inf_mag^2) / mu;
delta_max = 2*asin(1/e);

TR = delta_max / max(deltas);

feasible = TR > TR_min && AR > AR_min;
end

function ra = calc_ra(r_vec, v_vec, mu)
r = norm(r_vec);
v = norm(v_vec);

h_vec = cross(r_vec, v_vec);
h = norm(h_vec);

eps = v^2/2 - mu/r;                 % specific orbital energy
a = -mu/(2*eps);                    % semi-major axis (can be <0 if hyperbolic)

e = sqrt(1 + 2*eps*h^2/mu^2);        % eccentricity (works for all conics)

ra = a*(1+e);                            % aphelion radius (only meaningful if ellipse: e<1 and a>0)
end

function tof = identical_s_tof(S, p, h, s)
tof = (S*p - h/2) / s;
end