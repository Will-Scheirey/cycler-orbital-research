function [a, p, v_D, v_A] = lambert_solver(r1_vec, r2_vec, tof, mu, transfer_type)

r1 = norm(r1_vec);
r2 = norm(r2_vec);
c  = norm(r2_vec - r1_vec);
s  = 0.5*(c + r1 + r2);

% Transfer angle
cos_dtheta = dot(r1_vec, r2_vec) / (r1*r2);
cos_dtheta = max(-1, min(1, cos_dtheta));
dtheta = acos(cos_dtheta);

% Decide short-way vs long-way from transfer_type:
% "1*" = short-way, "2*" = long-way (your naming)
isLongWay = startsWith(transfer_type, "2");

alpha0 = @(a) 2 * asin( sqrt( s/(2*a) ) );
beta0m = @(a) 2 * asin( sqrt( (s-c)/(2*a) ) ); % magnitude

% Sign convention: long-way uses negative beta in classic formulation
beta0 = @(a) (isLongWay)*(-beta0m(a)) + (~isLongWay)*(beta0m(a));

% Time-of-flight equation (elliptic only here; your H cases omitted for clarity)
if     transfer_type == "1A"
    eqn = @(a) -sqrt(mu)*tof + a^(3/2) * ( (alpha0(a)-sin(alpha0(a))) - (beta0(a)-sin(beta0(a))) );
elseif transfer_type == "1B"
    eqn = @(a) -sqrt(mu)*tof + a^(3/2) * (2*pi - (alpha0(a)-sin(alpha0(a))) - (beta0(a)-sin(beta0(a))) );
elseif transfer_type == "2A"
    eqn = @(a) -sqrt(mu)*tof + a^(3/2) * ( (alpha0(a)-sin(alpha0(a))) + (beta0(a)-sin(beta0(a))) );
elseif transfer_type == "2B"
    eqn = @(a) -sqrt(mu)*tof + a^(3/2) * (2*pi - (alpha0(a)-sin(alpha0(a))) + (beta0(a)-sin(beta0(a))) );
else
    error("Invalid transfer type: " + transfer_type)
end

options = optimoptions('fsolve', 'MaxIterations', 1000, 'Display', 'off');

% IMPORTANT: enforce a >= s/2 to keep asin arguments <= 1
a_guess = s/2 * 1.05;
a = fsolve(@(x) real(eqn(max(x, s/2*1.000001))), a_guess, options);
a = max(a, s/2*1.000001);

alpha = alpha0(a);
beta  = beta0(a);

% p must match the same +/- structure
if isLongWay
    psi = (alpha + beta)/2;
else
    psi = (alpha - beta)/2;
end
p = (4*a*(s-r1)*(s-r2)/c^2) * sin(psi)^2;

A = sqrt(mu / (4*a)) * cot(alpha/2);
B = sqrt(mu / (4*a)) * cot(beta/2);

r1_hat = r1_vec / r1;
r2_hat = r2_vec / r2;
c_hat  = (r2_vec - r1_vec) / c;

v_D = (B + A)*c_hat + (B - A)*r1_hat;
v_A = (B + A)*c_hat - (B - A)*r2_hat;

% Ensure prograde (positive hz) if desired:
h = cross(r1_vec, v_D);
if h(3) < 0
    v_D = -v_D;
    v_A = -v_A;
end

end