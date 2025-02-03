%--------------------------------------------------------------------------
%
% test_MoonPosition: Computes lunar ephemeris
%
% References:
% Montenbruck O., Gill E.; Satellite Orbits: Models, Methods and 
% Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
%
% Montenbruck O., Pfleger T.; Astronomy on the Personal Computer; Springer 
% Verlag, Heidelberg; 4th edition (2000).
%
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill;
% New York; 4th edition (2013).
%
% https://ssd.jpl.nasa.gov/planets/eph_export.html
%
% https://celestrak.org/SpaceData/EOP-All.txt
% 
% Last modified:   2024/08/22   Meysam Mahooti
%
%--------------------------------------------------------------------------
clc
clear all
format long g

global PC    % Planetary Coefficients
global const % Astronomical Constants
SAT_Const

load DE440Coeff
PC = DE440Coeff;

% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k ~= 0)
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);
        for i=numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

MJD_UTC = Mjday(2024,8,22,18,59,0);
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
MJD_TT = MJD_UTC + TT_UTC/86400;
MJD_TDB = Mjday_TDB(MJD_TT);

fprintf('\nlunar coordinates derived from the Chebyshev coefficients of ');
fprintf('the Development Ephemeris DE440 [km]\n');
fprintf('Barycentric Dynamical Time (TDB) is used for JPL ephemerides computations\n');
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_TDB+2400000.5);
1e-3*r_Moon'

% Difference between ephemeris time and universal time
[year, month, day, hour, minute, sec] = invjday(MJD_UTC);
days = finddays(year, month, day, hour, minute, sec);
ET_UT = ETminUT(year+days/365.25);
MJD_ET = MJD_UTC+ET_UT/86400;

fprintf('\nlunar coordinates derived from the Chebyshev coefficients of ');
fprintf('the Development Ephemeris DE440 [km]\n');
fprintf('Ephemeris Time (ET) is used for JPL ephemerides computations\n');
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_ET+2400000.5);
1e-3*r_Moon'

fprintf('very accurate ELP2000-82 lunar coordinates [km]');
JD_TT = MJD_TT+2400000.5;
[X,Y,Z] = moonpos(JD_TT);
rMoon = [X,Y,Z]

fprintf('high-precision analytic lunar coordinates [km]');
rMoon = 1e-3*MoonBrown(MJD_TT)'

fprintf('low-precision analytic lunar coordinates [km]');
rMoon = 1e-3*Moon(MJD_TT)'

fprintf('Simpson analytic lunar coordinates [km]');
[rMoon,vMoon] = MoonSimpson(MJD_TT);
rMoon

