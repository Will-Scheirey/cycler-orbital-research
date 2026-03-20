function [r0_vec, v0_vec, w0] = calc_orbit_vecs(ai, ei, phi0, r0, mu, other_branch, retrograde)
if nargin < 6
    other_branch = false;
end

if nargin < 7
    retrograde = false;
end

v0 = sqrt(mu * (2/r0 - 1/ai));
epsilon = v0^2/2 - mu/r0;
h = sqrt(-1/2 * (mu^2/epsilon) * (1 - ei^2));

cos_theta = (ai*(1 - ei^2)/r0 - 1)/ei;
cos_theta = max(-1, min(1, cos_theta));

theta_i = acos(cos_theta);
if other_branch, theta_i = 2*pi - theta_i; end

w0 = phi0 - theta_i;

vr = mu/h *      ei*sin(theta_i);
vt = mu/h * (1 + ei*cos(theta_i));

vx = vr*cos(phi0) - vt*sin(phi0);
vy = vr*sin(phi0) + vt*cos(phi0);

r0_vec = r0 * [cos(phi0 ); sin(phi0); 0];
v0_vec = [vx; vy; 0];

if retrograde
    v0_vec = -v0_vec;
end
end