function nu = theta_from_tof(r0, v0, mu, dt)

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

E0 = atan2( sqrt(1 - e^2)*sin(nu0), e + cos(nu0) );
if E0 < 0, E0 = E0 + 2*pi; end

M0 = E0 - e*sin(E0);

n = sqrt(mu / a^3);
M = M0 + n*dt;
M = mod(M, 2*pi);

E = solve_kepler_elliptic(M, e);

nu = 2 * atan2( sqrt(1 + e)*sin(E/2), sqrt(1 - e)*cos(E/2) );
nu = mod(nu, 2*pi);

end

function E = solve_kepler_elliptic(M, e)
E = M;
for k = 1:50
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = -f/fp;
    E  = E + dE;
    if abs(dE) < 1e-12
        break
    end
end
end