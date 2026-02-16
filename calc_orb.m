function r_p_sol = calc_orb(mu_1, mu_2, r1, theta_r, v_inf_in_r, v_inf_in_t, theta, Vr_b2, Vt_b2)

syms r_p;

v_inf = sqrt(v_inf_in_r^2 + v_inf_in_t^2);

beta = 1 / (1 + (r_p*v_inf^2)/mu_2);
v_r = Vr_b2 + v_inf_in_r * (1 - 2*beta^2)            - v_inf_in_t * (2*beta*sqrt(1 - beta^2));
v_t = Vt_b2 + v_inf_in_r * (2*beta*sqrt(1 - beta^2)) + v_inf_in_t * (1 - 2*beta^2);

A = v_t*sin(theta) - v_r*cos(theta);

term = mu_1*sin(theta)^2 - r1^2*A^2 - r1*A*v_r * cos(theta_r);

eqn = term;
f = matlabFunction(eqn, 'Vars', {r_p});

% r_p_sol = fsolve(f, 1e3, optimset('Display', 'off'));
r_p_sol = fzero(f, 1e3);
end