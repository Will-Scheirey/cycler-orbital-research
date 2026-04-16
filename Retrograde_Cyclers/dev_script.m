clear; clc; close all

jupiter_moons = get_jupiter_moon_params();

num_v_inf = 100;
v_inf = linspace(0, 40, num_v_inf);

% Equation : rp = mu / (v_inf^2) * [1/sin(delta/2) - 1];
mu            = jupiter_moons.io.mu;
rp            = jupiter_moons.io.r;
e_io = rp .* v_inf.^2/mu + 1;
delta_io = 2 * asin(1 ./ e_io);

mu            = jupiter_moons.europa.mu;
rp            = jupiter_moons.europa.r;
e_eu = rp .* v_inf.^2/mu + 1;
delta_eu = 2 * asin(1 ./ e_eu);

mu            = jupiter_moons.ganymede.mu;
rp            = jupiter_moons.ganymede.r;
e_ga = rp .* v_inf.^2/mu + 1;
delta_ga = 2 * asin(1 ./ e_ga);

plot(v_inf, delta_io, 'DisplayName', 'Io', 'LineWidth', 1.5); hold on
plot(v_inf, delta_eu, 'DisplayName', 'Eu', 'LineWidth', 1.5); hold on
plot(v_inf, delta_ga, 'DisplayName', 'Ga', 'LineWidth', 1.5); hold on

legend

% Given a velocity vector with v_tangential and v_radial, to enter a
% retrograde orbit, we need to make v_tangential < 0. First, calculate the
% flight path angle:

vt = 5;
vr = 1;

gamma = atan2(vr, vt);

% For a valid retrograde orbit, the most fundamental controllable parameter
% is the periapsis radius wrt to the primary body - we need to make sure we 
% are not violating any constraints. And we can also choose apoapsis radius
% ra based on where we need to intersect in the future.

% Vis-Viva: v^2 = mu * (2/r - 1/a)
% We know r, and we know v = v(delta)

% We know rp = h^2 / mu * 1 / (1 + e)
% And h = r * vt, where we also know vt = vt(delta), so h = h(delta)
% 
% We can write e = (ra - rp) / (ra + rp), and a = (ra + rp)/2, so
% e = 2 * (ra - rp) / a

% rp = h(delta)^2 / mu * 1 / [1 + 2*(ra - rp) / a]
%
% We need to eliminate a from our equations. Rearranging,
%
% 1 + 2 * (ra - rp) / a = h(delta)^2 / (mu * rp)
% a = 2 * (ra - rp) / [h(delta^2) / (mu*rp) - 1]
%
% and
% a = 1 / (2/r - v(delta)^2 / mu)
%
% So 1 / (2/r - v(delta)^2 / mu) = 2 * (ra - rp) / [h(delta^2) / (mu*rp) - 1]
%
% 2/r - v(delta)^2 / mu = [h(delta^2) / (mu*rp) - 1] / [2 * (ra - rp)]
% Solve this eqn for delta

% First, find v(delta)
%
% The incoming v infinity is v_sc - v_b, or
%
% [v_b_r - v_sc_r, v_b_t - v_sc_t]

