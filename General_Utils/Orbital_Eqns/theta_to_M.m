function M = theta_to_M(theta, e)

E = 2*atan2( sqrt(1-e)*sin(theta/2), sqrt(1+e)*cos(theta/2) );
E = mod(E, 2*pi);

M = mod(E - e*sin(E), 2*pi);

end