% Based on: https://ssd.jpl.nasa.gov/planets/approx_pos.html

function out = planet_state_from_jd(params, planetName, JD_TDB_like)
%PLANET_STATE_FROM_JD
% Uses your struct fields:
%   a,e,i,L,w,W,omega plus secular rates *_dot
% where:
%   w = longitude of perihelion (varpi)
%   W = longitude of ascending node (Omega)
%
% Returns:
%   r_orb,v_orb (perifocal/orbital plane)
%   r_ecl,v_ecl (J2000 ecliptic)
%   r_eq, v_eq  (J2000 equatorial)

    % Sun GM (m^3/s^2)
    mu_sun = 1.32712440018e20;

    p0 = params.(planetName);

    wrap2pi = @(x) mod(x, 2*pi);
    wrapToPi = @(x) mod(x + pi, 2*pi) - pi;

    % 1) Julian centuries since J2000
    T = (JD_TDB_like - 2451545.0) / 36525.0;

    % 2) time-update the elements (secular)
    a     = p0.a + p0.a_dot * T;    % m
    e     = p0.e + p0.e_dot * T;    % -
    I     = p0.i + p0.i_dot * T;    % rad
    L     = p0.L + p0.L_dot * T;    % rad
    varpi = p0.w + p0.w_dot * T;    % rad (longitude of perihelion)
    Omega = p0.W + p0.W_dot * T;    % rad (longitude of node)

    % Argument of periapsis
    omega = wrap2pi(varpi - Omega);

    % Mean anomaly (no extra bT^2/cos/sin terms included)
    M = wrapToPi(L - varpi);

    % 3) Solve Kepler in radians: M = E - e sin E
    E = kepler_E_rad(M, e);

    % True anomaly
    nu = 2 * atan2( sqrt(1+e) * sin(E/2), sqrt(1-e) * cos(E/2) );
    nu = wrap2pi(nu);

    % Distance
    rmag = a * (1 - e*cos(E));          % m
    p    = a * (1 - e^2);               % semi-latus rectum, m

    % 4) Orbital plane (PQW / perifocal) position & velocity
    r_orb = [rmag*cos(nu); rmag*sin(nu); 0];

    v_orb = sqrt(mu_sun / p) * [-sin(nu); e + cos(nu); 0];  % m/s

    % 5) Rotate to J2000 ecliptic (column-vector active rotation)
    % r = Rz(Omega)*Rx(I)*Rz(omega)*r_orb
    R = Rz(Omega) * Rx(I) * Rz(omega);
    r_ecl = R * r_orb;
    v_ecl = R * v_orb;

    % 6) Ecliptic -> Equatorial using J2000 obliquity
    eps = deg2rad(23.43928);
    R_eps = Rx(eps);

    r_eq = R_eps * r_ecl;
    v_eq = R_eps * v_ecl;

    % Output
    out = struct();
    out.T = T;

    out.a = a;
    out.e = e;
    out.i = I;
    out.L = L;
    out.w = varpi;    % longitude of perihelion at date
    out.W = Omega;    % longitude of node at date
    out.omega = omega;

    out.M  = M;
    out.E  = E;
    out.nu = nu;

    out.r_orb = r_orb;
    out.v_orb = v_orb;

    out.r_ecl = r_ecl;
    out.v_ecl = v_ecl;

    out.r_eq  = r_eq;
    out.v_eq  = v_eq;
end

% ---------- helpers ----------

function E = kepler_E_rad(M, e)
% Newton solve, radians
    E = M; % good default
    for k = 1:50
        f  = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        dE = -f / fp;
        E  = E + dE;
        if abs(dE) < 1e-12
            break;
        end
    end
end

function R = Rz(a)
    ca = cos(a); sa = sin(a);
    R = [ ca -sa  0;
          sa  ca  0;
           0   0  1 ];
end

function R = Rx(a)
    ca = cos(a); sa = sin(a);
    R = [ 1  0   0;
          0  ca -sa;
          0  sa  ca ];
end