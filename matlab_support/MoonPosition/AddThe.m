%--------------------------------------------------------------------------
%
% AddThe: Calculates cos(alpha+beta) and sin(alpha+beta) using addition 
%         theorems
%
% Inputs:
%   c1        cos(alpha)
%   s1        sin(alpha)
%   c2        cos(beta)
%   s2        sin(beta)
%
% Outputs:
%   c         cos(alpha+beta)
%   s         sin(alpha+beta)
% 
% Last modified:   2018/01/27   Meysam Mahooti
%
%--------------------------------------------------------------------------
function [c, s] = AddThe(c1, s1, c2, s2)

c = c1 * c2 - s1 * s2; 
s = s1 * c2 + c1 * s2;

