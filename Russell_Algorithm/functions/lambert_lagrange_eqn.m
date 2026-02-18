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

if theta < pi
    beta = @(a) beta0(a);
else
    beta = @(a) -beta0(a);
end
end