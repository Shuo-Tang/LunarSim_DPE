function ArcSec=dms2as(D,M,S)
%-----------------------------------------------------------------
% dms2as function             Convert DD:MM:SS to arcsec.
% Input  : - Deg.
%          - Min.
%          - Sec.
% Output : - arcsec.
% Tested : Matlab 5.3
%     By : Eran O. Ofek      June 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------

ArcSec = D.*3600.0 + M.*60.0 + S;
