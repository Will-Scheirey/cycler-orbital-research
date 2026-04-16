function dt = tof_from_theta(r0, v0, mu, nu_target)

r0 = r0(:);
v0 = v0(:);

r = norm(r0);
v = norm(v0);

eps = 0.5*v^2 - mu/r;
a = -mu/(2*eps);
if a <= 0
    error('Orbit is not elliptical')
end

e_vec = ((v^2 - mu/r)*r0 - dot(r0,v0)*v0) / mu;
e = norm(e_vec);

cos_nu0 = dot(e_vec, r0) / (e*r);
cos_nu0 = max(-1, min(1, cos_nu0));

nu0 = acos(cos_nu0);
if dot(r0, v0) < 0
    nu0 = 2*pi - nu0;
end

% --- initial eccentric anomaly ---
E0 = atan2( sqrt(1 - e^2)*sin(nu0), e + cos(nu0) );
if E0 < 0, E0 = E0 + 2*pi; end

M0 = E0 - e*sin(E0);

% --- target eccentric anomaly ---
E = atan2( sqrt(1 - e^2)*sin(nu_target), e + cos(nu_target) );
if E < 0, E = E + 2*pi; end

M = E - e*sin(E);

% --- mean motion ---
n = sqrt(mu / a^3);

% --- Δt (handle wrapping forward in time) ---
dM = M - M0;
if dM < 0
    dM = dM + 2*pi;
end

dt = dM / n;

end