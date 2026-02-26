function r = r_from_tof(r0, v0, mu, dt)

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
h2 = dot(h,h);

e_vec = ((v0n^2 - mu/r0n)*r0 - dot(r0,v0)*v0) / mu;
e = norm(e_vec);

nu = true_anomaly_after_dt(r0, v0, mu, dt);

p = h2 / mu;
rmag = p / (1 + e*cos(nu));

e_hat = e_vec / e;
h_hat = h / sqrt(h2);
q_hat = cross(h_hat, e_hat);

r = rmag * (cos(nu)*e_hat + sin(nu)*q_hat);

end