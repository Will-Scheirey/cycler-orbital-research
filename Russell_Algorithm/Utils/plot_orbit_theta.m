function p = plot_orbit_theta(r0, v0, mu, theta_lim, varargin)
[~, ~, e, p, w, ~] = calc_orbit(r0, v0, mu);

thetas = linspace(theta_lim(1), theta_lim(2), 1000);

r = p ./ (1 + e*cos(thetas - w));

r_vec = r .* [cos(thetas); sin(thetas)];

if nargin < 4
    varargin = {'LineWidth', 1.5};
end

p = plot(r_vec(1, :), r_vec(2, :), varargin{:});
end