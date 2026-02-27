function [tof_eqn, alpha, beta] = lambert_lagrange_eqn(r1, r2, theta, mu, N, fast)
[alpha, beta] = get_lambert_angles(r1, r2, N, theta, fast);

tof_eqn = @(a) sqrt(a.^3 / mu) .* (2*N*pi + alpha(a) - beta(a) - sin(alpha(a)) + sin(beta(a)));
end