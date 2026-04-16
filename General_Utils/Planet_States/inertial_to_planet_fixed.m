function r_sc_rot = inertial_to_planet_fixed(r1, r2, r_sc, v1, v2, origin_type, m1, m2)

if nargin < 6 || isempty(origin_type)
    origin_type = "planet1";
end

if origin_type == "planet1"
    r_origin = r1;
elseif origin_type == "barycenter"
    r_origin = (m1*r1 + m2*r2) / (m1 + m2);
else
    error("origin_type must be 'planet1' or 'barycenter'");
end

xhat = (r2 - r1) / norm(r2 - r1);

if nargin >= 5 && ~isempty(v1) && ~isempty(v2)
    h = cross(r2 - r1, v2 - v1);
    zhat = h / norm(h);
else
    zhat = [0; 0; 1];
    zhat = zhat / norm(zhat);
end

yhat = cross(zhat, xhat);
yhat = yhat / norm(yhat);

% Re-orthogonalize zhat in case of numerical error
zhat = cross(xhat, yhat);
zhat = zhat / norm(zhat);

C_RI = [xhat, yhat, zhat];
C_IR = C_RI';

r_rel = r_sc - r_origin;
r_sc_rot = C_IR * r_rel;

end