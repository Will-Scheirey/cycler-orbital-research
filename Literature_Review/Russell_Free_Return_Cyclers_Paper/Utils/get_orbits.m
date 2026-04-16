function orbits = get_orbits(solution, mu, vd_generic, r_b1, v_b1, theta_generic)

num_s = solution.s;

theta_earth = theta_generic;

r0 = [r_b1; 0;    0];
v0 = [0;    v_b1; 0];

[~, ~, ~, ~, w, ~] = calc_orbit(r0, vd_generic, mu);

dt = tof_from_theta(r0, vd_generic, mu, theta_generic - w) + year2sec(floor(theta_generic / (2*pi)));

orb_idx = 1;
orbits{orb_idx} = {r0, vd_generic, dt, v0};

orb_idx = orb_idx + 1;

for s = 1:num_s

    v_inf_minus = solution.v_inf_minus{s};
    dt_years    = solution.dt_years{s};

    k_all = 1:length(v_inf_minus);
    num_k = length(k_all);

    if s > 1
        k_min = 1;
    else
        k_min = 2;
    end

    if s < num_s
        k_max = num_k;
    else
        k_max = num_k-1;
    end

    for k_idx = k_min:k_max
        k = k_all(k_idx);

        theta_earth = solution.theta_earth_all(k_idx);

        r1 = r_b1 * [cos(theta_earth);  sin(theta_earth); 0];
        v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];

        dt_leg = dt_years(k_idx - 1);

        v_local = solution.v_local{k};
        v_inf_global = local_vinf_to_global( ...
            v_local, ...
            solution.rotm, ...
            theta_earth, ...
            theta_generic);

        vd = v1 + v_inf_global;

        orbits{orb_idx} = {r1, vd, year2sec(dt_leg), v1};
        orb_idx = orb_idx + 1;
    end

    theta_earth = solution.theta_earth_all(num_k);

    r1 = r_b1 * [cos(theta_earth);  sin(theta_earth); 0];
    v1 = v_b1 * [-sin(theta_earth); cos(theta_earth); 0];

    Rz = [cos(theta_earth), -sin(theta_earth), 0;
        sin(theta_earth),  cos(theta_earth), 0;
        0,                 0,                1];

    vd = Rz * vd_generic;

    orbits{orb_idx} = {r1, vd, dt, v1};
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