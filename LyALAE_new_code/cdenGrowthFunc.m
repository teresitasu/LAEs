%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to be called by cdenGrowth.
%
function dfdap = cdenGrowthFunc(aexp);
global omega;
dfdap = 1./ (omega(1)./ aexp + omega(2).*aexp.*aexp + omega(3)).^1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
