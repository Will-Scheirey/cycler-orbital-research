% Adapted from 'Orbital Mechanics for Engineering Students' by Howard
% Curtis; Appendix D
function [V1_list, V2_list] = lambert_curtis(R1, R2, t, mu, z_lim, num_z)
%
%{
This function solves Lambert's problem.
mu            - gravitational parameter (km^3/s^2)
R1, R2        - initial and final position vectors (km)
r1, r2        - magnitudes of R1 and R2
t             - the time of flight from R1 to R2 (s)
V1, V2        - initial and final velocity vectors (km/s)
c12           - cross product of R1 into R2
theta         - angle between R1 and R2
string        - 'pro' if the orbit is prograde
                'retro' if the orbit is retrograde
A             - a constant given by Equation 5.35
z             - alpha*x^2, where alpha is the reciprocal of the
                semimajor axis and x is the universal anomaly
y(z)          - a function of z given by Equation 5.38
F(z,t)        - a function of the variable z and constant t,
                given by Equation 5.40
dFdz(z)       - the derivative of F(z,t), given by Equation 5.43
ratio         - F/dFdz
tol           - tolerance on precision of convergence
nmax          - maximum number of iterations of Newton's procedure
f, g          - Lagrange coefficients
gdot          - time derivative of g
C(z), S(z)    - Stumpff functions
dum           - a dummy variable
User M-functions required: stumpC and stumpS
%}
% ––––––––––––––––––––––––––––––––––––––––––––––

%...Magnitudes of R1 and R2:
r1 = norm(R1);
r2 = norm(R2);

c12 = cross(R1, R2);
theta0 = acos(dot(R1,R2)/r1/r2);
thetas = [theta0, 2*pi - theta0];  % short and long transfer

V1_list = {};
V2_list = {};

T_list = [];

for th = thetas
    theta = th;

    %...Equation 5.35:
    A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));

    %{
    num_z = 1000;
    z_lim = 1200;
    %}
    zs = linspace(-z_lim, z_lim, num_z);

    options = optimset('Display','off');

    for zz = 1:num_z
        %...Determine approximately where F(z,t) changes sign, and
        %...use that value of z as the starting value for Equation 5.45:
        z0 = zs(zz);
        if ~isreal(F(z0,t)), continue; end
        try
            z = fzero(@(z)F(z,t), z0, options);
            % store zsol if new
        catch
            continue;
        end
        if ~isfinite(z), continue; end

        %...Set an error tolerance and a limit on the number of iterations:
        tol = 1.e-8;
        nmax = 5000;

        %...Iterate on Equation 5.45 until z is determined to within the
        %...error tolerance:
        options = optimset('TolX', tol, 'Display', 'off', 'MaxIter', nmax);
        z = fzero(@(z) F(z,t), z, options);

        %...Equation 5.46a:
        f = 1 - y(z)/r1;
        %...Equation 5.46b:
        g = A*sqrt(y(z)/mu);
        %...Equation 5.46d:
        gdot = 1 - y(z)/r2;
        if abs(g) < eps, continue; end

        %...Equation 5.28:
        V1 = 1/g*(R2 - f*R1);
        %...Equation 5.29:
        V2 = 1/g*(gdot*R2 - R1);

        if true
            k_hat = [0; 0; 1];
            h_vec = cross(R1, V1);
            if dot(h_vec, k_hat) < 0
                % retrograde -> skip
                continue
            end
        end

        T = orbital_period(R1, V1, mu);

        if any(abs(T_list - T) < t/1000000)
            continue
        end

        T_list(end+1) = T;

        % collect both branches
        V1_list{end+1} = V1;
        V2_list{end+1} = V2;

        % V1_list{end+1} = -V1;
        % V2_list{end+1} = -V2;
    end
end
%
% Subfunctions used in the main body:
%
%...Equation 5.38:
    function dum = y(z)
        dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
    end
%...Equation 5.40:
    function dum = F(z,t)
        dum = (y(z)/C(z))^ 1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
    end
%...Equation 5.43:
    function dum = dFdz(z)
        if z == 0
            dum = sqrt(2)/40*y(0)^ 1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
        else
            dum = (y(z)/C(z))^ 1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
                + 3*S(z)^ 2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...
                + A*sqrt(C(z)/y(z)));
        end
    end
%...Stumpff functions:
    function dum = C(z)
        dum = stumpC(z);
    end

    function dum = S(z)
        dum = stumpS(z);
    end

end %lambert
%

%
function c = stumpC(z)
%
%{
This function evaluates the Stumpff function C(z) according
to Equation 3.53.
z - input argument
c - value of C(z)
User M-functions required: none
%}
% ––––––––––––––––––––––––––––––––––––––––––––––
if z > 0
    c = (1 - cos(sqrt(z)))/z;

elseif z < 0
    c = (cosh(sqrt(-z)) - 1)/(-z);

else
    c = 1/2;
end
end
%

%
function s = stumpS(z)
%
%{
This function evaluates the Stumpff function S(z) according
to Equation 3.52.
z - input argument
s - value of S(z)
User M-functions required: none
%}
% ––––––––––––––––––––––––––––––––––––––––––––––
if z > 0
    s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^ 3;

elseif z < 0
    s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^ 3;
else
    s = 1/6;
end
%
end


function T = orbital_period(r_vec, v_vec, mu)
r = norm(r_vec);
v = norm(v_vec);
epsilon = v^2/2 - mu/r;
a = -mu/(2*epsilon);
T = 2*pi*sqrt(a^3/mu);
end