function [x,y,z,times] = get_planet_pos_time(r0, v0, mu, t_lim)

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

times = linspace(t_lim(1), t_lim(2), 1e4);

x = zeros(size(times));
y = zeros(size(times));
z = zeros(size(times));

% circular orbit
if e < 1e-10
    p_hat = r0 / norm(r0);
    q_hat = cross(h_hat, p_hat);

    n = sqrt(mu / r0n^3);
    theta0 = 0;

    for i = 1:numel(times)
        theta = theta0 + n*times(i);
        r = r0n * (cos(theta)*p_hat + sin(theta)*q_hat);

        x(i) = r(1);
        y(i) = r(2);
        z(i) = r(3);
    end
    return
end

% non-circular elliptic orbit
dotrv0 = dot(r0, v0);

cosE0 = (1 - r0n/a) / e;
sinE0 = dotrv0 / (e * sqrt(mu*a));
E0 = atan2(sinE0, cosE0);
M0 = E0 - e*sin(E0);

n = sqrt(mu/a^3);

for i = 1:numel(times)
    dt = times(i);
    M = M0 + n*dt;

    E = M;
    for k = 1:100
        fval = E - e*sin(E) - M;
        fp   = 1 - e*cos(E);
        dE   = -fval/fp;
        E    = E + dE;

        if abs(dE) < 1e-12
            break
        end
    end

    dE = E - E0;

    f = 1 - a/r0n * (1 - cos(dE));
    g = dt + sqrt(a^3/mu) * (sin(dE) - dE);

    r = f*r0 + g*v0;

    x(i) = r(1);
    y(i) = r(2);
    z(i) = r(3);
end
end