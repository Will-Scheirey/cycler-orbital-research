function plot_orbit_from_state(r0, v0, mu, dt, varargin)
% Plot a Keplerian orbit arc starting from (r0,v0) for time dt.
% Works for elliptical (a>0) and hyperbolic (a<0). Parabolic not supported.

r0 = r0(:);
v0 = v0(:);

r = norm(r0);
v = norm(v0);

h_vec = cross(r0, v0);
h = norm(h_vec);
if h < 1e-12
    error('Degenerate orbit: |r x v| is ~0')
end

eps = 0.5*v^2 - mu/r;
if abs(eps) < 1e-12
    error('Parabolic case not supported')
end

a = -mu/(2*eps);

e_vec = ((v^2 - mu/r)*r0 - dot(r0,v0)*v0) / mu;
e = norm(e_vec);

% For (near) circular, periapsis direction is ill-defined. Still plot using plane basis.
is_circ = (e < 1e-10);

h_hat = h_vec / h;

if ~is_circ
    p_hat = e_vec / e;            % points to periapsis in orbital plane
    q_hat = cross(h_hat, p_hat);  % completes right-handed basis in plane
else
    p_hat = r0 / norm(r0);        % arbitrary in-plane axis (along current r)
    q_hat = cross(h_hat, p_hat);
    q_hat = q_hat / norm(q_hat);
end

% Initial true anomaly (or in-plane angle for circular)
if ~is_circ
    cos_nu0 = dot(e_vec, r0) / (e*r);
    cos_nu0 = max(-1, min(1, cos_nu0));
    nu0 = acos(cos_nu0);
    if dot(r0, v0) < 0
        nu0 = 2*pi - nu0;
    end
else
    % Angle in (p_hat,q_hat) plane
    nu0 = atan2(dot(r0, q_hat), dot(r0, p_hat));
end

% "Mean motion"-like scale
if a > 0
    n = sqrt(mu / a^3);
else
    n = sqrt(mu / (-a)^3);
end

% Sample times
N = 1200;
t_vec = linspace(0, dt, N);
r_hist = zeros(3, N);

% Precompute for ellipse
if a > 0 && ~is_circ
    E0 = atan2(sqrt(1 - e^2)*sin(nu0), e + cos(nu0));
    if E0 < 0, E0 = E0 + 2*pi; end
    M0 = E0 - e*sin(E0);
elseif a < 0 && ~is_circ
    % Hyperbolic anomaly at t0
    F0 = asinh( sin(nu0)*sqrt(e^2 - 1)/(1 + e*cos(nu0)) );
    M0 = e*sinh(F0) - F0;
else
    M0 = 0; % circular: just advance angle with n*t
end

for k = 1:N
    t = t_vec(k);

    if a > 0
        % Elliptical
        if is_circ
            nu = nu0 + n*t;
            p = a; % for circular, radius is constant = a
            r_pf = p * [cos(nu); sin(nu); 0];
        else
            M = M0 + n*t;
            E = solve_kepler_elliptic(M, e);
            nu = atan2( sqrt(1 - e^2)*sin(E), cos(E) - e );
            p = a*(1 - e^2);
            r_pf = (p / (1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
        end

    else
        % Hyperbolic (non-circular always)
        M = M0 + n*t;
        F = solve_kepler_hyperbolic(M, e);
        nu = atan2( sqrt(e^2 - 1)*sinh(F), e - cosh(F) );
        p = a*(1 - e^2); % note: a<0, (1-e^2)<0 so p>0
        r_pf = (p / (1 + e*cos(nu))) * [cos(nu); sin(nu); 0];
    end

    r_hist(:,k) = p_hat*r_pf(1) + q_hat*r_pf(2);
end

% Plot
plot3(r_hist(1,:), r_hist(2,:), r_hist(3,:), 'k', 'LineWidth', 1.5, varargin{:});
hold on
% plot3(r0(1), r0(2), r0(3), 'go', 'MarkerFaceColor','g');               % start
% plot3(r_hist(1,end), r_hist(2,end), r_hist(3,end), 'ro', 'MarkerFaceColor','r'); % end
% plot3(0,0,0,'bo','MarkerFaceColor','b');                               % central body

end

% --------- local functions (same file) ----------

function E = solve_kepler_elliptic(M, e)
M = mod(M, 2*pi);
E = M + sign(sin(M))*0.85*e; % better than E=M for higher e

for k = 1:60
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = -f/fp;
    E  = E + dE;
    if abs(dE) < 1e-12
        break
    end
end
end

function F = solve_kepler_hyperbolic(M, e)
% Solve e*sinh(F) - F = M
F = asinh(M/e); % decent initial guess

for k = 1:60
    f  = e*sinh(F) - F - M;
    fp = e*cosh(F) - 1;
    dF = -f/fp;
    F  = F + dF;
    if abs(dF) < 1e-12
        break
    end
end
end