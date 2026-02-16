clear; clc; close all

mu_sun = 132712000000;
mu = mu_sun;

AU = 149.6e6;

r0x = AU;
r0y = 0;
v0_0 = sqrt(mu_sun / r0x);

mu_b2 = 42828;
r_b2 = 227.9e6;

%{
mu_b2 = 126686000;
r_b2 = 778.6e6;
%}

v_b2 = sqrt(mu_sun / r_b2);
h_b2 = r_b2*v_b2;
e_b2 = 0;

rp = 279e5;

num_t = 10;
num_v = 10;

t_all = linspace(0, 2*pi, num_t);
v_all = linspace(2, 10, num_v);

valid = nan(num_t, num_v);

figure(1)
clf

for t_idx = 1:num_t
    t = t_all(t_idx);
    for v_idx = 1:num_v
        v_mag = v_all(v_idx);
        v0x = v_mag * sin(t);
        v0y =  v0_0 + v_mag * cos(t);

        [r_sc, h_sc, e_sc, p_sc, w_sc, nu_sc] = calc_orbit([r0x; r0y; 0], [v0x; v0y; 0], mu, 0);

        % [r_out, h_out, e_out, p_out, w_out, nu_out, theta_intersect] = calc_grav_assist(r0x, r0y, v0x, v0y, mu, mu_b2, r_b2, v_b2, h_b2, e_b2, rp, 1, 1);

        rps = linspace(1e0, 1e8, 100);

        [rp_helio] = calc_grav_assist2(r0x, r0y, v0x, v0y, mu, mu_b2, r_b2, v_b2, h_b2, e_b2, rp, rps, 1, 2);

        if ~isempty(rp_helio)
            plot3(ones(100, 1) * v_mag, rps, rp_helio / AU); hold on
        end

        %{
        if isempty(r_out)
            continue;
        end

        % Plotting

        theta_sc = linspace(0, theta_intersect, 1000);
        nu_sc = theta_sc - w_sc;
        r1 = p_sc ./ (1 + e_sc*cos(nu_sc));

        theta_out = linspace(theta_intersect, theta_intersect + 2*pi, 1000);
        nu_out = theta_out - w_out;
        r2 = p_out ./ (1 + e_out*cos(nu_out));

        theta_planet = linspace(0,2*pi,1000);

        clf
        plot(r1.*cos(theta_sc), r1.*sin(theta_sc), '--g', 'LineWidth', 1.5); hold on
        plot(r2.*cos(theta_out), r2.*sin(theta_out), '-g', 'LineWidth', 1.5); hold on

        plot(r0x*cos(theta_planet), r0x*sin(theta_planet), '-b'); hold on
        plot(r_b2*cos(theta_planet), r_b2*sin(theta_planet), '-r')

        % plot(r_intersect*cos(theta_intersect), r_intersect*sin(theta_intersect), 'mx')

        axis equal
        xlim([-1,1]*r_b2 * 4)
        ylim([-1,1]*r_b2 * 4)

        title(sprintf("t: %d, v: %d", t_idx, v_idx))
        %}
    end
end

xlabel("V Mag")
ylabel("Gravity Assist Rp (km)")
zlabel("Resulting Heliocentric Rp (AU)")


function [r_out, h_out, e_out, p_out, w_out, nu_out, theta_intersect] = calc_grav_assist(r0x, r0y, v0x, v0y, mu, mu_b2, r_b2, v_b2, h_b2, e_b2, rp, sign, angle)

[r_sc, h_sc, e_sc, p_sc, w_sc, nu_sc] = calc_orbit([r0x; r0y; 0], [v0x; v0y; 0], mu, 0);
a_sc = p_sc / (1 - e_sc^2);

r_out = [];
h_out = [];
e_out = [];
p_out = [];
w_out = [];
nu_out = [];
theta_intersect = [];

% Intersect
c = (p_sc/r_b2 - 1)/e_sc;
if abs(c) > 1, return; end
nu_int = acos(c);

if imag(nu_int) ~= 0
    return
end

% Calculate both points
phi_a = w_sc + nu_int;
phi_b = w_sc - nu_int;

phi_a = mod(phi_a, 2*pi);
phi_b = mod(phi_b, 2*pi);

% Choose one
if angle == 1
    theta_intersect = phi_a;
else
    theta_intersect = phi_b;
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

% rp = calc_orb(mu, mu_b2, r0x, deg2rad(10), v_inf_in(1), v_inf_in(2), theta_intersect, vr_b2, vt_b2);
% rp = calc_orb_numeric(mu, mu_b2, r0x, deg2rad(90), v_inf_in(1), v_inf_in(2), ...
%                      theta_intersect, vr_b2, vt_b2, sign);

opts = struct('rp_min',1,'rp_max',1e7,'n_grid',10000,'verbose',false);

theta_E_target = deg2rad(10); % inertial polar angle where you want to hit Earth's orbit
sign_turn = +1;               % pick one branch and stick to it
which_intersection = 2;       % choose encounter point

[r_out, h_out, e_out, p_out, w_out, nu_out, theta_intersect, rp_sol, info] = ...
    calc_grav_assist_direct(r0x, r0y, v0x, v0y, mu, mu_b2, r_b2, ...
                            theta_E_target, sign_turn, which_intersection, opts);

if ~isfinite(rp_sol)
    % No valid root in the bracket range
    return
end

rp = rp_sol


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

[r_out, h_out, e_out, p_out, w_out, nu_out] = calc_orbit(r_vec_out, v_vec_out, mu, theta_intersect);

% TOF
E = 2 * atan(sqrt((1 - e_sc) / (1 + e_sc)) * tan(nu_int / 2));
Me = fsolve(@(M) E - e_sc * sin(E) - M, 0, optimset('Display', 'off'));
tof_1 = sqrt(a_sc^3 / mu) * Me;
end

function [rp_helio] = calc_grav_assist2(r0x, r0y, v0x, v0y, mu, mu_b2, r_b2, v_b2, h_b2, e_b2, rp, rps, sign, angle)

rp_helio = [];

[r_sc, h_sc, e_sc, p_sc, w_sc, nu_sc] = calc_orbit([r0x; r0y; 0], [v0x; v0y; 0], mu, 0);
a_sc = p_sc / (1 - e_sc^2);

% Intersect
c = (p_sc/r_b2 - 1)/e_sc;
if abs(c) > 1, return; end
nu_int = acos(c);

if imag(nu_int) ~= 0
    return
end

% Calculate both points
phi_a = w_sc + nu_int;
phi_b = w_sc - nu_int;

phi_a = mod(phi_a, 2*pi);
phi_b = mod(phi_b, 2*pi);

% Choose one
if angle == 1
    theta_intersect = phi_a;
else
    theta_intersect = phi_b;
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

for i = 1:length(rps)
    rp = rps(i);

    rp_helio(i) = calc_helio_rp(mu, mu_b2, rp, deg2rad(25), v_inf_in(1), v_inf_in(2), theta_intersect, vr_b2, vt_b2);
end

end

function [r, h, e, p, w, nu] = calc_orbit(r_vec, v_vec, mu, theta)
h_vec = cross(r_vec, v_vec);
e_vec = cross(v_vec, h_vec) / mu - r_vec/norm(r_vec);

h = norm(h_vec);
e = norm(e_vec);

p = h^2/mu;

w = atan2(e_vec(2), e_vec(1));

nu = theta - w;

r = p ./ (1 + e*cos(nu));
end
