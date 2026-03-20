function dt = tof_between_thetas(theta0, theta, e, n, retrograde)

M0 = theta_to_M(theta0, e);
M  = theta_to_M(theta,  e);

if retrograde
    dM = mod(M0 - M, 2*pi);
else
    dM = mod(M - M0, 2*pi);
end

dt = dM / n;
end