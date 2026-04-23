function orbits = get_orbits(solution, mu, vd_generic, r_b1, v_b1, theta_generic)

num_s = solution.s;

r0 = [r_b1; 0;    0];
v0 = [0;    v_b1; 0];

[~, ~, ~, ~, w, ~] = calc_orbit(r0, vd_generic, mu);
dt = tof_from_theta(r0, vd_generic, mu, theta_generic - w) + year2sec(floor(theta_generic / (2*pi)));

orb_idx = 1;
orbits{orb_idx} = {r0, vd_generic, dt, v0};
orb_idx = orb_idx + 1;

for s = 1:num_s

    v_inf_minus_s = solution.v_inf_minus{s};
    v_local_s     = solution.v_local{s};
    dt_years_s    = solution.dt_years{s};

    if s < num_s
        dt_years_s = [dt_years_s; dt];
    end

    theta_s       = solution.theta_earth_all{s};

    num_k = length(v_inf_minus_s);

    k_min = 2;
    k_max = num_k - 1;

    for k = k_min:k_max
        theta_earth = theta_s(k);

        r1 = r_b1 * [cos(theta_earth);  sin(theta_earth); 0];
        v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];

        dt_leg = dt_years_s(k - 1);

        v_local = v_local_s{k};
        v_inf_global = local_vinf_to_global( ...
            v_local, ...
            solution.rotm, ...
            theta_earth, ...
            theta_generic);

        vd = v1 + v_inf_global;

        orbits{orb_idx} = {r1, vd, year2sec(dt_leg), v1};
        orb_idx = orb_idx + 1;
    end

    theta_earth = theta_s(num_k);

    r1 = r_b1 * [cos(theta_earth);  sin(theta_earth); 0];
    v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];

    Rz = [cos(theta_earth), -sin(theta_earth), 0;
        sin(theta_earth),  cos(theta_earth), 0;
        0,                 0,                1];

    vd = Rz * vd_generic;

    orbits{orb_idx} = {r1, vd, dt, v1};
    orb_idx = orb_idx + 1;
end
end

function tf = is_odd_half_year(dt_year)
n_half = round(2*dt_year);
tf = (mod(n_half, 2) == 1);
end

function v_global = local_vinf_to_global(v_local, C0, theta_earth, theta_generic)
dtheta = theta_earth - theta_generic;

Rz = [cos(dtheta), -sin(dtheta), 0;
    sin(dtheta),  cos(dtheta), 0;
    0,            0,           1];

v_global = Rz * C0 * v_local;
end