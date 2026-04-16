function r = get_planet_pos(r0, v0, mu, t)

r0 = r0(:);
v0 = v0(:);

r0n = norm(r0);
v0n = norm(v0);

eps = 0.5*v0n^2 - mu/r0n;
a = -mu/(2*eps);

if a <= 0
    error('Orbit is not elliptical')
end

h = cross(r0, v0);
h_hat = h / norm(h);

e_vec = ((v0n^2 - mu/r0n)*r0 - dot(r0,v0)*v0)/mu;
e = norm(e_vec);

% Treat nearly circular orbits as circular
if e < 1e-6
    p_hat = r0 / norm(r0);
    q_hat = cross(h_hat, p_hat);
    q_hat = q_hat / norm(q_hat);

    n = sqrt(mu / a^3);
    theta = n*t;

    r = a*cos(theta)*p_hat + a*sin(theta)*q_hat;
    return
end

p_hat = e_vec / e;
q_hat = cross(h_hat, p_hat);
q_hat = q_hat / norm(q_hat);

dotrv0 = dot(r0, v0);

cosE0 = (1 - r0n/a) / e;
sinE0 = dotrv0 / (e * sqrt(mu*a));
E0 = atan2(sinE0, cosE0);
M0 = E0 - e*sin(E0);

n = sqrt(mu / a^3);
M = M0 + n*t;
M = mod(M, 2*pi);

E = M;
for k = 1:50
    fval = E - e*sin(E) - M;
    fp   = 1 - e*cos(E);
    dE   = -fval/fp;
    E    = E + dE;
    if abs(dE) < 1e-13
        break
    end
end

x_pf = a*(cos(E) - e);
y_pf = a*sqrt(1 - e^2)*sin(E);

r = x_pf*p_hat + y_pf*q_hat;
end