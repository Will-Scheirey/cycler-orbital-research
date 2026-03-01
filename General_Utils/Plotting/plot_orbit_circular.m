function p = plot_orbit_circular(r, mu, varargin)
    v0 = sqrt(mu / r);

    v0_vec = [0; v0; 0];
    r0_vec = [r; 0; 0];

    p = plot_orbit(r0_vec, v0_vec, mu, varargin{:});
end