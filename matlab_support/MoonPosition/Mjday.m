%------------------------------------------------------------------------------
%
% Mjday: Modified Julian Date from calendar date and time
%
% Inputs:
%   Year      Calendar date components
%   Month
%   Day
%   Hour      Time components
%   Min
%   Sec
% Output:
%   Modified Julian Date
%
% Last modified:   2022/09/24   Meysam Mahooti
%
% Reference:
% Montenbruck O., Gill E.; Satellite Orbits: Models, Methods and 
% Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
%
%------------------------------------------------------------------------------
function Mjd = Mjday(Year, Month, Day, Hour, Minute, Sec)

if (nargin < 4)
    Hour   = 0;
    Minute = 0;
    Sec    = 0;
end

if (Month<=2)
    Month = Month+12;
    Year  = Year-1;
end

if ( (10000*Year+100*Month+Day) <= 15821004 )
    b = (-2 + fix((Year+4716)/4) - 1179);        % Julian calendar
else
    b = fix(Year/400)-fix(Year/100)+fix(Year/4); % Gregorian calendar
end    

MjdMidnight = 365*Year - 679004 + b + fix(30.6001*(Month+1)) + Day;
FracOfDay   = (Hour+Minute/60+Sec/3600)/24;

Mjd = MjdMidnight + FracOfDay;

