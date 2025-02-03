%------------------------------------------------------------------------------
%
% invjday: Calendar date and time from Modified Julian Date
%
% Input:
%   Mjd       Modified Julian Date
% Outputs:
%   Year      Calendar date components
%   Month
%   Day
%   Hour      Time components
%   Min
%   Sec
%
% Last modified:   2022/05/27   Meysam Mahooti
%
% Reference:
% Montenbruck O., Gill E.; Satellite Orbits: Models, Methods and 
% Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
%
%------------------------------------------------------------------------------
function [Year, Month, Day, Hour, Minute, Sec] = invjday(Mjd)

% Convert Julian day number to calendar date
a = fix(Mjd+2400001.0);

if ( a < 2299161 )  % Julian calendar
    c = a + 1524;
else                % Gregorian calendar
    b = fix((a-1867216.25)/36524.25);
    c = a +  b - fix(b/4) + 1525;
end

d     = fix ( (c-122.1)/365.25 );
e     = 365*d + fix(d/4);
f     = fix ( (c-e)/30.6001 );

Day   = c - e - fix(30.6001*f);
Month = f - 1 - 12*fix(f/14);
Year  = d - 4715 - fix((7+Month)/10);

Hours = 24.0*(Mjd-floor(Mjd));

Hour = fix(Hours);
x = (Hours-Hour)*60.0;
Minute = fix(x);
Sec = (x-Minute)*60.0;

