function plot_orbit_from_state(r0, v0, mu, line, lineWidth, end_theta)

% constants and vectors
h = cross(r0, v0);
h_norm = norm(h);

h_hat = h / h_norm;

% eccentricity vector
e_vec = cross(v0, h)/mu - r0/norm(r0);
e = norm(e_vec);

% semi-major axis from energy
energy = norm(v0)^2/2 - mu/norm(r0);

if e > 1
    a = h_norm^2/mu * 1/(e^2 - 1);
else
    a = -mu/(2*energy);
end

% periapsis direction: guard for circular orbits
if e > 1e-10
    p_hat = e_vec / e;
else
    % nearly circular: use radial direction for periapsis (arbitrary but continuous)
    p_hat = r0 / norm(r0);
end

% q_hat completes the right-handed triad (perifocal frame)
q_hat = cross(h_hat, p_hat);
q_hat = q_hat / norm(q_hat);

% rotation matrix from perifocal (PQW) to inertial
R_pqw = [p_hat, q_hat, h_hat];   % 3x3

% sample true anomaly and radius in perifocal plane
thetas = linspace(0, 2*pi, 400);            % 1xN
r_pf_mag = a * (1 - e^2) ./ (1 + e .* cos(thetas));  % 1xN

% perifocal coordinates (3xN)
r_pf = [r_pf_mag .* cos(thetas);    % x_P = r cos(theta)
    r_pf_mag .* sin(thetas);    % y_P = r sin(theta)
    zeros(1, numel(thetas))];   % z_P = 0

% transform to inertial: (3x3) * (3xN) -> 3xN
r_eci = R_pqw * r_pf;

% r_eci = r_eci(:, vecnorm(r_eci) > 149.6e6);

% plot (assumes hold is on if you want multi plots)
plot(r_eci(1,:), r_eci(2,:), line, 'LineWidth', lineWidth)
% optionally mark periapsis direction and departure point:
% plot([0 p_hat(1)*max(r_pf_mag)], [0 p_hat(2)*max(r_pf_mag)], '--k');
% plot(r0(1), r0(2), 'og', 'MarkerFaceColor','g');
end