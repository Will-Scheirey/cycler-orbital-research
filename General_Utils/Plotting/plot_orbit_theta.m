function hPlot = plot_orbit_theta(r0, v0, mu, theta_lim, varargin)

[~, ~, e, p_orbit, w, ~] = calc_orbit(r0, v0, mu);

t0 = theta_lim(1);
t1 = theta_lim(2);

t0 = mod(t0, 2*pi);
t1 = mod(t1, 2*pi);

% force forward sweep (short-ish way in the forward direction)
if t1 < t0
    t1 = t1 + 2*pi;
end

thetas = linspace(t0, t1, 1000);

% conic in polar form: r(nu) = p / (1 + e cos nu)
% here nu = theta - w (theta is inertial angle)
r = p_orbit ./ (1 + e*cos(thetas - w));

x = r .* cos(thetas);
y = r .* sin(thetas);

if isempty(varargin)
    varargin = {'LineWidth', 1.5};
end

hPlot = plot(x, y, varargin{:});

end