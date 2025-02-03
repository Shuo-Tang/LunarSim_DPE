%--------------------------------------------------------------------------
%
% MoonPos: Computes the Moon's equatorial position using Brown's theory
%          (Improved Lunar Ephemeris)
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
% 
% Output:
%   rMoon     Geocentric equatorial position of the Moon (in [m]), referred
%             to the mean equator and equinox of J2000 (EME2000, ICRF)
%
% Notes: Light-time is already taken into account
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function rMoon = MoonBrown(Mjd_TT)

pi2       = 2*pi;        % 2pi
Arcs      = 3600*180/pi; % Arcseconds per radian
MJD_J2000 = 51544.5;     % Modif. Julian Date of J2000.0
R_Earth   = 6378.137e3;  % Earth's radius [m]; WGS-84

global Dlam gam1C DS sinPi N Cos Sin

T   = (Mjd_TT-MJD_J2000)/36525; % Julian cent. since J2000
o   = 6;                        % Index offset
dim = 2*o+1;                    % Work array dimension
T2  = T*T;                      % Time

% Cosine and sine of mean arguments
Cos = zeros(dim,4);
Sin = zeros(dim,4);

% Reset perturbations
Dlam=0; DS=0; gam1C=0; sinPi=3422.7000; N=0;

% Longperiodic perturbations
S1 = Sine (0.19833+0.05611*T); S2 = Sine (0.27869+0.04508*T);
S3 = Sine (0.16827-0.36903*T); S4 = Sine (0.34734-5.37261*T);
S5 = Sine (0.10498-5.37899*T); S6 = Sine (0.42681-0.41855*T);
S7 = Sine (0.14943-5.37511*T); 

dL0 = 0.84*S1+0.31*S2+14.27*S3+ 7.26*S4+ 0.28*S5+0.24*S6;
dl  = 2.94*S1+0.31*S2+14.27*S3+ 9.34*S4+ 1.12*S5+0.83*S6;
dls =-6.40*S1                                   -1.89*S6;
dF  = 0.21*S1+0.31*S2+14.27*S3-88.70*S4-15.30*S5+0.24*S6-1.86*S7;
dD  = dL0-dls;

Dgam = -3332e-9 * Sine (0.59734-5.37261*T)...
       -539e-9 * Sine (0.35498-5.37899*T)...
       -64e-9 * Sine (0.39943-5.37511*T);

% Mean arguments of the lunar orbit (incl. longperiodic corrections)
% L0 mean longitude of the Moon
% l  mean anomaly of the Moon     l' mean anomaly of the Sun      
% F  mean distance from the node  D  mean elongation from the Sun 

L0 = pi2*Frac(0.60643382+1336.85522467*T-0.00000313*T2) + dL0/Arcs;
l  = pi2*Frac(0.37489701+1325.55240982*T+0.00002565*T2) + dl /Arcs;
ls = pi2*Frac(0.99312619+  99.99735956*T-0.00000044*T2) + dls/Arcs;
F  = pi2*Frac(0.25909118+1342.22782980*T-0.00000892*T2) + dF /Arcs;
D  = pi2*Frac(0.82736186+1236.85308708*T-0.00000397*T2) + dD /Arcs;


% Cosine and sine of multiples of mean arguments 
% incl. secular correction
for i=0:3

switch i
    case 0
        arg=l;  max=4; fac=1.000002208;
    case 1
        arg=ls; max=3; fac=0.997504612-0.002495388*T;
    case 2
        arg=F;  max=4; fac=1.000002708+139.978*Dgam;
    case 3
        arg=D;  max=6; fac=1;
end

Cos(o+1,i+1) = 1;  Cos(o+2,i+1) = cos(arg)*fac;  Cos(o,i+1) = +Cos(o+2,i+1);
Sin(o+1,i+1) = 0;  Sin(o+2,i+1) = sin(arg)*fac;  Sin(o,i+1) = -Sin(o+2,i+1);

for j=2:max
    [Cos(o+j+1,i+1),Sin(o+j+1,i+1)] = AddThe(Cos(o+j,i+1),Sin(o+j,i+1),Cos(o+2,i+1),Sin(o+2,i+1));
    Cos(o-j+1,i+1) = + Cos(o+j+1,i+1);
    Sin(o-j+1,i+1) = - Sin(o+j+1,i+1);
end

end

  % Solar perturbations
  PertAddSol (   13.902,   14.06,-0.001,   0.2607,0, 0, 0, 4);
  PertAddSol (    0.403,   -4.01,+0.394,   0.0023,0, 0, 0, 3);
  PertAddSol ( 2369.912, 2373.36,+0.601,  28.2333,0, 0, 0, 2);
  PertAddSol ( -125.154, -112.79,-0.725,  -0.9781,0, 0, 0, 1);
  PertAddSol (    1.979,    6.98,-0.445,   0.0433,1, 0, 0, 4);
  PertAddSol (  191.953,  192.72,+0.029,   3.0861,1, 0, 0, 2);
  PertAddSol (   -8.466,  -13.51,+0.455,  -0.1093,1, 0, 0, 1);
  PertAddSol (22639.500,22609.07,+0.079, 186.5398,1, 0, 0, 0);
  PertAddSol (   18.609,    3.59,-0.094,   0.0118,1, 0, 0,-1);
  PertAddSol (-4586.465,-4578.13,-0.077,  34.3117,1, 0, 0,-2);
  PertAddSol (   +3.215,    5.44,+0.192,  -0.0386,1, 0, 0,-3);
  PertAddSol (  -38.428,  -38.64,+0.001,   0.6008,1, 0, 0,-4);
  PertAddSol (   -0.393,   -1.43,-0.092,   0.0086,1, 0, 0,-6);
  PertAddSol (   -0.289,   -1.59,+0.123,  -0.0053,0, 1, 0, 4);
  PertAddSol (  -24.420,  -25.10,+0.040,  -0.3000,0, 1, 0, 2);
  PertAddSol (   18.023,   17.93,+0.007,   0.1494,0, 1, 0, 1);
  PertAddSol ( -668.146, -126.98,-1.302,  -0.3997,0, 1, 0, 0);
  PertAddSol (    0.560,    0.32,-0.001,  -0.0037,0, 1, 0,-1);
  PertAddSol ( -165.145, -165.06,+0.054,   1.9178,0, 1, 0,-2);
  PertAddSol (   -1.877,   -6.46,-0.416,   0.0339,0, 1, 0,-4);
  PertAddSol (    0.213,    1.02,-0.074,   0.0054,2, 0, 0, 4);
  PertAddSol (   14.387,   14.78,-0.017,   0.2833,2, 0, 0, 2);
  PertAddSol (   -0.586,   -1.20,+0.054,  -0.0100,2, 0, 0, 1);
  PertAddSol (  769.016,  767.96,+0.107,  10.1657,2, 0, 0, 0);
  PertAddSol (   +1.750,    2.01,-0.018,   0.0155,2, 0, 0,-1);
  PertAddSol ( -211.656, -152.53,+5.679,  -0.3039,2, 0, 0,-2);
  PertAddSol (   +1.225,    0.91,-0.030,  -0.0088,2, 0, 0,-3);
  PertAddSol (  -30.773,  -34.07,-0.308,   0.3722,2, 0, 0,-4);
  PertAddSol (   -0.570,   -1.40,-0.074,   0.0109,2, 0, 0,-6);
  PertAddSol (   -2.921,  -11.75,+0.787,  -0.0484,1, 1, 0, 2);
  PertAddSol (   +1.267,    1.52,-0.022,   0.0164,1, 1, 0, 1);
  PertAddSol ( -109.673, -115.18,+0.461,  -0.9490,1, 1, 0, 0);
  PertAddSol ( -205.962, -182.36,+2.056,  +1.4437,1, 1, 0,-2);
  PertAddSol (    0.233,    0.36, 0.012,  -0.0025,1, 1, 0,-3);
  PertAddSol (   -4.391,   -9.66,-0.471,   0.0673,1, 1, 0,-4);
  PertAddSol (    0.283,    1.53,-0.111,  +0.0060,1,-1, 0,+4);
  PertAddSol (   14.577,   31.70,-1.540,  +0.2302,1,-1, 0, 2);
  PertAddSol (  147.687,  138.76,+0.679,  +1.1528,1,-1, 0, 0);
  PertAddSol (   -1.089,    0.55,+0.021,   0.0   ,1,-1, 0,-1);
  PertAddSol (   28.475,   23.59,-0.443,  -0.2257,1,-1, 0,-2);
  PertAddSol (   -0.276,   -0.38,-0.006,  -0.0036,1,-1, 0,-3);
  PertAddSol (    0.636,    2.27,+0.146,  -0.0102,1,-1, 0,-4);
  PertAddSol (   -0.189,   -1.68,+0.131,  -0.0028,0, 2, 0, 2);
  PertAddSol (   -7.486,   -0.66,-0.037,  -0.0086,0, 2, 0, 0);
  PertAddSol (   -8.096,  -16.35,-0.740,   0.0918,0, 2, 0,-2);
  PertAddSol (   -5.741,   -0.04, 0.0  ,  -0.0009,0, 0, 2, 2);
  PertAddSol (    0.255,    0.0 , 0.0  ,   0.0   ,0, 0, 2, 1);
  PertAddSol ( -411.608,   -0.20, 0.0  ,  -0.0124,0, 0, 2, 0);
  PertAddSol (    0.584,    0.84, 0.0  ,  +0.0071,0, 0, 2,-1);
  PertAddSol (  -55.173,  -52.14, 0.0  ,  -0.1052,0, 0, 2,-2);
  PertAddSol (    0.254,    0.25, 0.0  ,  -0.0017,0, 0, 2,-3);
  PertAddSol (   +0.025,   -1.67, 0.0  ,  +0.0031,0, 0, 2,-4);
  PertAddSol (    1.060,    2.96,-0.166,   0.0243,3, 0, 0,+2);
  PertAddSol (   36.124,   50.64,-1.300,   0.6215,3, 0, 0, 0);
  PertAddSol (  -13.193,  -16.40,+0.258,  -0.1187,3, 0, 0,-2);
  PertAddSol (   -1.187,   -0.74,+0.042,   0.0074,3, 0, 0,-4);
  PertAddSol (   -0.293,   -0.31,-0.002,   0.0046,3, 0, 0,-6);
  PertAddSol (   -0.290,   -1.45,+0.116,  -0.0051,2, 1, 0, 2);
  PertAddSol (   -7.649,  -10.56,+0.259,  -0.1038,2, 1, 0, 0);
  PertAddSol (   -8.627,   -7.59,+0.078,  -0.0192,2, 1, 0,-2);
  PertAddSol (   -2.740,   -2.54,+0.022,   0.0324,2, 1, 0,-4);
  PertAddSol (    1.181,    3.32,-0.212,   0.0213,2,-1, 0,+2);
  PertAddSol (    9.703,   11.67,-0.151,   0.1268,2,-1, 0, 0);
  PertAddSol (   -0.352,   -0.37,+0.001,  -0.0028,2,-1, 0,-1);
  PertAddSol (   -2.494,   -1.17,-0.003,  -0.0017,2,-1, 0,-2);
  PertAddSol (    0.360,    0.20,-0.012,  -0.0043,2,-1, 0,-4);
  PertAddSol (   -1.167,   -1.25,+0.008,  -0.0106,1, 2, 0, 0);
  PertAddSol (   -7.412,   -6.12,+0.117,   0.0484,1, 2, 0,-2);
  PertAddSol (   -0.311,   -0.65,-0.032,   0.0044,1, 2, 0,-4);
  PertAddSol (   +0.757,    1.82,-0.105,   0.0112,1,-2, 0, 2);
  PertAddSol (   +2.580,    2.32,+0.027,   0.0196,1,-2, 0, 0);
  PertAddSol (   +2.533,    2.40,-0.014,  -0.0212,1,-2, 0,-2);
  PertAddSol (   -0.344,   -0.57,-0.025,  +0.0036,0, 3, 0,-2);
  PertAddSol (   -0.992,   -0.02, 0.0  ,   0.0   ,1, 0, 2, 2);
  PertAddSol (  -45.099,   -0.02, 0.0  ,  -0.0010,1, 0, 2, 0);
  PertAddSol (   -0.179,   -9.52, 0.0  ,  -0.0833,1, 0, 2,-2);
  PertAddSol (   -0.301,   -0.33, 0.0  ,   0.0014,1, 0, 2,-4);
  PertAddSol (   -6.382,   -3.37, 0.0  ,  -0.0481,1, 0,-2, 2);
  PertAddSol (   39.528,   85.13, 0.0  ,  -0.7136,1, 0,-2, 0);
  PertAddSol (    9.366,    0.71, 0.0  ,  -0.0112,1, 0,-2,-2);
  PertAddSol (    0.202,    0.02, 0.0  ,   0.0   ,1, 0,-2,-4);
  PertAddSol (    0.415,    0.10, 0.0  ,   0.0013,0, 1, 2, 0);
  PertAddSol (   -2.152,   -2.26, 0.0  ,  -0.0066,0, 1, 2,-2);
  PertAddSol (   -1.440,   -1.30, 0.0  ,  +0.0014,0, 1,-2, 2);
  PertAddSol (    0.384,   -0.04, 0.0  ,   0.0   ,0, 1,-2,-2);
  PertAddSol (   +1.938,   +3.60,-0.145,  +0.0401,4, 0, 0, 0);
  PertAddSol (   -0.952,   -1.58,+0.052,  -0.0130,4, 0, 0,-2);
  PertAddSol (   -0.551,   -0.94,+0.032,  -0.0097,3, 1, 0, 0);
  PertAddSol (   -0.482,   -0.57,+0.005,  -0.0045,3, 1, 0,-2);
  PertAddSol (    0.681,    0.96,-0.026,   0.0115,3,-1, 0, 0);
  PertAddSol (   -0.297,   -0.27, 0.002,  -0.0009,2, 2, 0,-2);
  PertAddSol (    0.254,   +0.21,-0.003,   0.0   ,2,-2, 0,-2);
  PertAddSol (   -0.250,   -0.22, 0.004,   0.0014,1, 3, 0,-2);
  PertAddSol (   -3.996,    0.0 , 0.0  ,  +0.0004,2, 0, 2, 0);
  PertAddSol (    0.557,   -0.75, 0.0  ,  -0.0090,2, 0, 2,-2);
  PertAddSol (   -0.459,   -0.38, 0.0  ,  -0.0053,2, 0,-2, 2);
  PertAddSol (   -1.298,    0.74, 0.0  ,  +0.0004,2, 0,-2, 0);
  PertAddSol (    0.538,    1.14, 0.0  ,  -0.0141,2, 0,-2,-2);
  PertAddSol (    0.263,    0.02, 0.0  ,   0.0   ,1, 1, 2, 0);
  PertAddSol (    0.426,   +0.07, 0.0  ,  -0.0006,1, 1,-2,-2);
  PertAddSol (   -0.304,   +0.03, 0.0  ,  +0.0003,1,-1, 2, 0);
  PertAddSol (   -0.372,   -0.19, 0.0  ,  -0.0027,1,-1,-2, 2);
  PertAddSol (   +0.418,    0.0 , 0.0  ,   0.0   ,0, 0, 4, 0);
  PertAddSol (   -0.330,   -0.04, 0.0  ,   0.0   ,3, 0, 2, 0);
  
  % Solar perturbations in latitude
  PertAddN(-526.069, 0, 0,1,-2);   PertAddN(  -3.352, 0, 0,1,-4);
  PertAddN( +44.297,+1, 0,1,-2);   PertAddN(  -6.000,+1, 0,1,-4);
  PertAddN( +20.599,-1, 0,1, 0);   PertAddN( -30.598,-1, 0,1,-2);
  PertAddN( -24.649,-2, 0,1, 0);   PertAddN(  -2.000,-2, 0,1,-2);
  PertAddN( -22.571, 0,+1,1,-2);   PertAddN( +10.985, 0,-1,1,-2);
  
  % Planetary perturbations
  PertPlanetary(T); 
  
  % Position vector
  S   = F + DS/Arcs;
  fac = 1.000002708+139.978*Dgam;
  
  rMoon = Polar(mod(L0+Dlam/Arcs,pi2),(fac*(18518.511+1.189+gam1C)*...
                sin(S)-6.24*sin(3*S)+N)/Arcs,R_Earth*Arcs/(sinPi*0.999953253));
  
  rMoon = (EclMatrix(Mjd_TT)*PrecMatrix(MJD_J2000,Mjd_TT))'*rMoon;
  
end

%
% AddSol: Summation of solar perturbations
%
function [] = PertAddSol(coeffl, coeffS, coeffg, coeffP, p, q, r, s)

global Dlam gam1C DS sinPi

[x, y] = PertTerm (p,q,r,s);
Dlam   = Dlam  + coeffl*y; DS    = DS    + coeffS*y;
gam1C  = gam1C + coeffg*x; sinPi = sinPi + coeffP*x;

end

%
% AddN: Summation of perturbation in latitude
%
function [] = PertAddN(coeffN, p, q, r, s)

global N

[~,y] = PertTerm(p,q,r,s); 
N = N + coeffN*y;

end

%
% Term: calculates x=cos(p*l+q*ls+r*F+s*D) and y=sin(p*l+q*ls+r*F+s*D)
%
%
function [x, y] = PertTerm(p, q, r, s)

global Cos Sin

o    = 6;    % Index offset
i = zeros(4,1);
i(1)=p; i(2)=q; i(3)=r; i(4)=s;  x=1.0; y=0.0;

for k=1:4
    if (i(k)~=0)
        [x,y] = AddThe(x,y,Cos(o+1+i(k),k),Sin(o+1+i(k),k));
    end
end

end

%
% Planetary: Perturbations of ecliptic latitude by Venus and Jupiter
%
%
function [] = PertPlanetary(T)

global Dlam

Dlam = Dlam+0.82*Sine(0.7736  -62.5512*T)+0.31*Sine(0.0466 -125.1025*T)...
           +0.35*Sine(0.5785  -25.1042*T)+0.66*Sine(0.4591+1335.8075*T)...
           +0.64*Sine(0.3130  -91.5680*T)+1.14*Sine(0.1480+1331.2898*T)...
           +0.21*Sine(0.5918+1056.5859*T)+0.44*Sine(0.5784+1322.8595*T)...
           +0.24*Sine(0.2275   -5.7374*T)+0.28*Sine(0.2965   +2.6929*T)...
           +0.33*Sine(0.3132   +6.3368*T);
end

%
% Sine
%
function out = Sine(x)
    out = sin(2*pi*Frac(x));
end

