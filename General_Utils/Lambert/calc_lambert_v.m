function [v_D, v_A] = calc_lambert_v(mu, a, c, alpha, beta, r1_vec, r2_vec)

A = sqrt(mu / (4*a)) * cot(alpha/2);
B = sqrt(mu / (4*a)) * cot(beta/2);

r1_hat = r1_vec / norm(r1_vec);
r2_hat = r2_vec / norm(r2_vec);
c_hat  = (r2_vec - r1_vec) / c;

v_D = (B + A)*c_hat + (B - A)*r1_hat;
v_A = (B + A)*c_hat - (B - A)*r2_hat;
%{
h = cross(r1_vec, v_D);

if h(3) < 0
    v_D = -v_D;
    v_A = -v_A;
end
%}
end