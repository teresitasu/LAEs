%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%*******************************************************************
%*  [D, dlnDdt] = cdenGrowth(aexp)                                 *
%*******************************************************************
%*******************************************************************
%
% ARGUMENTS
%  aexp            Cosmological expansion factor aexp = 1/ (1 + z).
%
% RETURNS
%  D               Growth factor, normalised arbitrarily to initial a_i = 0.001.
%  dlnDdt          d(lnD)/dt for growth factor D, normalised
%                  arbitrarily to initial a_i = 0.001.
% Based on Peebles LSSU Eq.(10.12).
% If aexp is an array, it must be ordered from small to large values and
% with a spacing sufficient for good accuracy using the trapezoidal rule.
%
% COMPATIBILITY: Matlab
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  22 01 10 Creation date.
%  05 10 11 Modified to include velocity growth term.
%  16 11 17 Modified to allow for vector aexp.
%
function [D, dlnDdt] = cdenGrowth(aexp);
global omega;
aexp0 = 0.001;
gint = 0;
lena = length(aexp);
disp(lena)
if(lena > 1)
  daexp = gradient(aexp);
  Y = daexp.*cdenGrowthFunc(aexp);
  gint = cumtrapz(Y);
end
dfdap = @cdenGrowthFunc;
g1 = quad(dfdap,aexp0,aexp(1));
gint = g1 + gint;
D = gint.*(omega(1)./aexp.^3 + omega(2) + omega(3)./aexp.^2).^0.5;
dlnDdt = -1.5*omega(1)./aexp.^3 - omega(3)./aexp.^2;
dlnDdt = dlnDdt./ (omega(1)./aexp.^3 + omega(2) + omega(3)./aexp.^2).^0.5;
dlnDdt = dlnDdt + 1./ (gint.*(omega(1)./aexp + omega(2).*aexp.^2 + omega(3)));

