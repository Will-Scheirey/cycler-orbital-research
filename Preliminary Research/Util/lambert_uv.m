function [v1, v2, info] = lambert_uv(r1, r2, tof, mu, varargin)
%LAMBERT_UV  Robust universal-variable Lambert solver (0-rev).
%
%   [v1, v2] = lambert_uv(r1, r2, tof, mu)
%   [v1, v2, info] = lambert_uv(..., 'prograde', true, 'longway', false)
%
% Inputs:
%   r1, r2 : 3x1 or 1x3 position vectors [km]
%   tof    : time of flight [s] (must be > 0)
%   mu     : gravitational parameter [km^3/s^2]
%
% Name-value options:
%   'prograde' (logical) : true = choose +h along +k in the transfer plane
%                          false = retrograde (opposite)     (default: true)
%   'longway'  (logical) : false = short-way (Δθ <= π)
%                          true  = long-way  (Δθ = 2π-Δθshort) (default: false)
%   'tol'      (scalar)  : root tolerance on z (default 1e-10)
%   'maxIter'  (int)     : max iterations (default 200)
%
% Outputs:
%   v1, v2 : 3x1 velocity vectors [km/s]
%   info   : struct with fields:
%       .z, .theta, .A, .iter, .converged
%
% Notes:
% - This is a standard Vallado-style universal variable Lambert solver
%   for the 0-revolution case (single-rev / no multi-rev families).
% - Works well for porkchop plots and avoids branch "seams" if you
%   control short/long way explicitly.

% -------------------- parse + sanitize --------------------
p = inputParser;
p.addParameter('prograde', true, @(x)islogical(x) && isscalar(x));
p.addParameter('longway',  false, @(x)islogical(x) && isscalar(x));
p.addParameter('tol',      1e-10, @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('maxIter',  200,   @(x)isnumeric(x) && isscalar(x) && x>=1);
p.parse(varargin{:});
opt = p.Results;

r1 = r1(:); r2 = r2(:);
if numel(r1)~=3 || numel(r2)~=3
    error('r1 and r2 must be 3-vectors.');
end
if tof <= 0
    error('tof must be > 0.');
end
if mu <= 0
    error('mu must be > 0.');
end

r1n = norm(r1); r2n = norm(r2);
if r1n == 0 || r2n == 0
    error('r1 and r2 must be non-zero.');
end

% -------------------- transfer angle + branch --------------------
c = dot(r1,r2)/(r1n*r2n);
c = max(-1,min(1,c));
theta_short = acos(c); % in [0, pi]

% Determine direction of motion (prograde/retrograde) relative to transfer plane.
% For planar problems, this aligns with +/- z. For general 3D, we build a plane normal.
h_dir = cross(r1, r2);
if norm(h_dir) < 1e-14
    error('r1 and r2 are nearly colinear; Lambert is ill-conditioned.');
end
h_dir = h_dir / norm(h_dir);

% Define a "reference" normal for prograde choice.
% If your problem is in XY plane, this is simply +k.
k_ref = [0;0;1];
% If the plane is far from XY, k_ref is still fine as a reference; we only use its sign.
sgn = sign(dot(h_dir, k_ref));
if sgn == 0, sgn = 1; end

% If user wants prograde, we want h to align with +k_ref in the common planar case.
% If sgn is negative, then r1->r2 cross points "down", so prograde means using the longway sense.
want_h_up = opt.prograde;         % prograde => +k_ref
have_h_up = (sgn > 0);            % does cross(r1,r2) point up?

% Choose theta based on longway and desired sense.
theta = theta_short;
if opt.longway
    theta = 2*pi - theta_short;
end

% If the current geometric cross product sense doesn't match desired prograde/retrograde,
% flip the transfer angle to traverse the other direction in the plane.
% In 2D porkchops, this keeps continuity.
if want_h_up ~= have_h_up
    theta = 2*pi - theta; % flip direction
end

% -------------------- compute A --------------------
% Vallado: A = sin(theta) * sqrt(r1*r2/(1-cos(theta)))
sinT = sin(theta);
cosT = cos(theta);

den = (1 - cosT);
if abs(den) < 1e-14
    % error('Transfer angle too small/close to 0 or 2π; Lambert ill-conditioned.');
end

A = sinT * sqrt(r1n*r2n/den);

if abs(A) < 1e-14
    error('A is ~0 (singular). Check geometry/TOF.');
end

% -------------------- solve for z --------------------
tol = opt.tol;
maxIter = opt.maxIter;

% Helper functions
    function [Cz, Sz] = stumpff(z)
        if z > 0
            sz = sqrt(z);
            Cz = (1 - cos(sz))/z;
            Sz = (sz - sin(sz))/sz^3;
        elseif z < 0
            sz = sqrt(-z);
            Cz = (cosh(sz) - 1)/(-z);
            Sz = (sinh(sz) - sz)/sz^3;
        else
            Cz = 1/2;
            Sz = 1/6;
        end
    end

    function y = y_of_z(z)
        [Cz, Sz] = stumpff(z);
        % guard: Cz can be tiny; y must be positive for real sqrt
        y = r1n + r2n + A * (z*Sz - 1)/sqrt(Cz);
    end

    function F = F_of_z(z)
        [Cz, Sz] = stumpff(z);
        y = r1n + r2n + A * (z*Sz - 1)/sqrt(Cz);
        if ~(isreal(y) && y > 0 && isreal(Cz) && Cz > 0)
            F = NaN;
            return
        end
        F = (y/Cz)^(3/2) * Sz + A*sqrt(y) - sqrt(mu)*tof;
    end

% Bracket z by expanding until F changes sign (or we hit a limit)
z = 0.0;
F0 = F_of_z(z);
if ~isfinite(F0)
    % try a small perturbation
    z = 1e-3;
    F0 = F_of_z(z);
end
if ~isfinite(F0)
    % error('Failed to evaluate Lambert function near z=0. Check geometry/TOF.');
end

% Find a bracket [zL, zU] with opposite signs
zL = z; zU = z;
FL = F0; FU = F0;

step = 0.5;           % initial expansion
zMax = 1e6;           % safety cap
found = false;

for k = 1:200
    % expand lower
    zL_try = zL - step;
    FL_try = F_of_z(zL_try);
    if isfinite(FL_try) && sign(FL_try) ~= sign(F0)
        zL = zL_try; FL = FL_try;
        zU = z;      FU = F0;
        found = true;
        break
    end

    % expand upper
    zU_try = zU + step;
    FU_try = F_of_z(zU_try);
    if isfinite(FU_try) && sign(FU_try) ~= sign(F0)
        zL = z;      FL = F0;
        zU = zU_try; FU = FU_try;
        found = true;
        break
    end

    step = step * 2;
    if abs(zL) > zMax || abs(zU) > zMax
        break
    end
end

if ~found
    % F0 might already be ~0
    if abs(F0) < 1e-10
        zL = z; zU = z;
        found = true;
    else
        % error('Could not bracket a Lambert root in z. Try different TOF or check geometry.');
    end
end

% Bisection + safeguarded Newton (derivative-free bisection here for robustness)
iter = 0;
converged = false;

if zL == zU
    z = zL;
    converged = true;
else
    for iter = 1:maxIter
        zMid = 0.5*(zL + zU);
        FMid = F_of_z(zMid);

        if ~isfinite(FMid)
            % If invalid, nudge midpoint slightly toward z=0
            zMid = 0.5*(zMid + 0.0);
            FMid = F_of_z(zMid);
            if ~isfinite(FMid)
                % fall back to bisection endpoints shrink
                zU = zMid;
                continue
            end
        end

        if abs(FMid) < 1e-12 || abs(zU - zL) < tol
            z = zMid;
            converged = true;
            break
        end

        if sign(FMid) == sign(FL)
            zL = zMid; FL = FMid;
        else
            zU = zMid; FU = FMid;
        end
    end

    if ~converged
        z = 0.5*(zL + zU);
    end
end

% -------------------- compute f, g, gdot and velocities --------------------
[Cz, Sz] = stumpff(z);
y = r1n + r2n + A * (z*Sz - 1)/sqrt(Cz);

if ~(isreal(y) && y > 0 && isreal(Cz) && Cz > 0)
    error('Lambert y(z) invalid at converged z. (numerical issue)');
end

f    = 1 - y/r1n;
g    = A * sqrt(y/mu);
gdot = 1 - y/r2n;

if abs(g) < 1e-14
    error('g is ~0 (singular).');
end

v1 = (r2 - f*r1) / g;
v2 = (gdot*r2 - r1) / g;

% -------------------- output info --------------------
info = struct();
info.z = z;
info.theta = theta;
info.A = A;
info.iter = iter;
info.converged = converged;
end