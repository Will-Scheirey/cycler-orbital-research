function [T12,T13,T23,T3,err] = calc_synodic3(a, mu, maxYears, tolRad, dt, tMinDays)

if nargin < 3 || isempty(maxYears), maxYears = 200; end
if nargin < 4 || isempty(tolRad),   tolRad   = 1*pi/180; end
if nargin < 5 || isempty(dt),       dt       = 86400; end
if nargin < 6 || isempty(tMinDays), tMinDays = 30; end   % ignore first 30 days

n = sqrt(mu ./ a.^3);

T12 = 2*pi / abs(n(1)-n(2));
T13 = 2*pi / abs(n(1)-n(3));
T23 = 2*pi / abs(n(2)-n(3));

dn12 = n(1) - n(2);
dn13 = n(1) - n(3);

tMax = maxYears * 365.25 * 86400;
tMin = tMinDays * 86400;

bestE = inf; bestT = NaN; bestE12 = NaN; bestE13 = NaN;

for t = tMin:dt:tMax
    e12 = mod(dn12*t, 2*pi); e12 = min(e12, 2*pi - e12);
    e13 = mod(dn13*t, 2*pi); e13 = min(e13, 2*pi - e13);
    e   = max(e12, e13);

    if e < bestE
        bestE = e; bestT = t; bestE12 = e12; bestE13 = e13;
    end
end

T3  = bestT;
err = [bestE bestE12 bestE13];

end