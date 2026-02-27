function [theta, r2] = get_intersection(r0, vd, mu, r_b2)
[~, ~, e, p_generic, w] = calc_orbit(r0, vd, mu);

theta = theta_from_r(p_generic, e, w, r_b2);
r2 = r_b2;

if imag(theta) ~= 0
    ra = ra_from_rv(r0, vd, mu);
    theta = real(theta_from_r(p_generic, e, w, ra));
    r2 = ra;
end