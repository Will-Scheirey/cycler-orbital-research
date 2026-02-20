function [r, h, e, p, w, nu] = calc_orbit(r_vec, v_vec, mu)

h_vec = cross(r_vec, v_vec);
e_vec = cross(v_vec, h_vec) / mu - r_vec/norm(r_vec);

h = norm(h_vec);
e = norm(e_vec);

p = h^2/mu;

w = atan2(e_vec(2), e_vec(1));

theta = atan2(r_vec(2), r_vec(1));
nu = theta - w;

r = p ./ (1 + e*cos(nu));
end