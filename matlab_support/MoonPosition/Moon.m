%--------------------------------------------------------------------------
%
% Moon: Computes the Moon's geocentric position using a low precision
%       analytical series
%
% Input:
%   Mjd_TT    Terrestrial Time (Modified Julian Date)
%
% Output:
%   rMoon     Lunar position vector [m] with respect to the mean equator 
%             and equinox of J2000 (EME2000, ICRF)
%
% Last modified:   2022/12/04   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function rMoon = Moon(Mjd_TT)

global const

% Constants
ep = 23.43929111*const.Rad; % Obliquity of J2000 ecliptic

T = (Mjd_TT-const.MJD_J2000)/36525; % Julian cent. since J2000

% Mean elements of lunar orbit
L_0 =           Frac(0.606433 + 1336.851344*T); % Mean longitude [rev] w.r.t. J2000 equinox
l   = const.pi2*Frac(0.374897 + 1325.552410*T); % Moon's mean anomaly [rad]
lp  = const.pi2*Frac(0.993133 +   99.997361*T); % Sun's mean anomaly [rad]
D   = const.pi2*Frac(0.827361 + 1236.853086*T); % Diff. long. Moon-Sun [rad]
F   = const.pi2*Frac(0.259086 + 1342.227825*T); % Argument of latitude 


% Ecliptic longitude (w.r.t. equinox of J2000)
dL = +22640*sin(l) - 4586*sin(l-2*D) + 2370*sin(2*D) +  769*sin(2*l)...
     -668*sin(lp) - 412*sin(2*F) - 212*sin(2*l-2*D) - 206*sin(l+lp-2*D)...
     +192*sin(l+2*D) - 165*sin(lp-2*D) - 125*sin(D) - 110*sin(l+lp)...
     +148*sin(l-lp) - 55*sin(2*F-2*D);

L = const.pi2 * Frac( L_0 + dL/1296.0e3 ); % [rad]

% Ecliptic latitude
S  = F + (dL+412*sin(2*F)+541*sin(lp))/const.Arcs;
h  = F-2*D;
N  = -526*sin(h) + 44*sin(l+h) - 31*sin(-l+h) - 23*sin(lp+h)...
     +11*sin(-lp+h) - 25*sin(-2*l+F) + 21*sin(-l+F);

B = ( 18520.0*sin(S) + N )/const.Arcs; % [rad]

cosB = cos(B);

% Distance [m]
R = 385000e3 - 20905e3*cos(l) - 3699e3*cos(2*D-l) - 2956e3*cos(2*D)...
    -570e3*cos(2*l) + 246e3*cos(2*l-2*D) - 205e3*cos(lp-2*D)...
    -171e3*cos(l+2*D) - 152e3*cos(l+lp-2*D);

% Equatorial coordinates
rMoon = R_x(-ep) * [R*cos(L)*cosB, R*sin(L)*cosB, R*sin(B)]';

