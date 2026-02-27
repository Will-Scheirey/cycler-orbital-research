function ra = ra_from_rv(r_vec, v_vec, mu)
r = norm(r_vec);
v = norm(v_vec);

h_vec = cross(r_vec, v_vec);
h = norm(h_vec);

eps = v^2/2 - mu/r;
a = -mu/(2*eps);

e = sqrt(1 + 2*eps*h^2/mu^2);

ra = a*(1+e);
end