%--------------------------------------------------------------------------
% Az  : azimuth of vector
% Elev: altitude of vector
% r   : norm of vector
%
% Last modified:   2018/01/27   Meysam Mahooti
%--------------------------------------------------------------------------
function Vec = Polar(Az, Elev, r)

if (nargin<3)
    r = 1;
end
Vec = zeros(3,1);
cosEl = cos(Elev);
Vec(1) = r * cos(Az) * cosEl;
Vec(2) = r * sin(Az) * cosEl;
Vec(3) = r * sin(Elev);