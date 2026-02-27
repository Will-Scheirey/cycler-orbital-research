function p = plot_orbit(r0, v0, mu, varargin)
    p = plot_orbit_theta(r0, v0, mu, [0, 2*pi], varargin{:});
end