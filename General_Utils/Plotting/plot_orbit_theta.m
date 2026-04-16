function hPlot = plot_orbit_theta(r0, v0, mu, theta_lim, varargin)

r0 = r0(:);
v0 = v0(:);

h = cross(r0, v0);
h_hat = h / norm(h);

e_vec = ((norm(v0)^2 - mu/norm(r0)) * r0 - dot(r0, v0) * v0) / mu;
e = norm(e_vec);

if e < 1e-12
    p_hat = r0 / norm(r0);
else
    p_hat = e_vec / e;
end

q_hat = cross(h_hat, p_hat);

p_orbit = norm(h)^2 / mu;

t0 = theta_lim(1);
t1 = theta_lim(2);

if t1 < t0
    t1 = t1 + 2*pi;
end

thetas = linspace(t0, t1, 1000);

r = p_orbit ./ (1 + e*cos(thetas));

R = p_hat .* (r .* cos(thetas)) + q_hat .* (r .* sin(thetas));

x = R(1, :);
y = R(2, :);
z = R(3, :);

if isempty(varargin)
    varargin = {'LineWidth', 1.5};
end

hPlot = plot3(x, y, z, varargin{:});

end