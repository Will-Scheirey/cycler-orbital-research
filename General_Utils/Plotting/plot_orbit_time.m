function [hPlot, x, y, z, times] = plot_orbit_time(r0, v0, mu, t_lim, varargin)

[x,y,z,times] = get_planet_pos_time(r0, v0, mu, t_lim);
hPlot = plot3(x, y, z, varargin{:});

end