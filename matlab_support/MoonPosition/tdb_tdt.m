function Diff=tdb_tdt(JD)
%-----------------------------------------------------------------
% tdb_tdt function             Convert TDT to TDB or vice versa.
% Input  : - Vector of julian days.
% Output : - TDB-TDT, [seconds].
% Tested : Matlab 5.3
%     By : Eran O. Ofek      June 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------
RAD = 57.29577951308232;
if (nargin==1),
   Option = 'b';
end

G = 357.53 + 0.9856003.*(JD - 2451545.0);
G = G./RAD;
Diff = 0.001658.*sin(G) + 0.000014.*sin(2.*G);  %sec.

