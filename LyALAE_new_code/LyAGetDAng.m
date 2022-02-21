%
%*******************************
%*  [DAng,h] = LyAGetDAng(zred)
%*******************************
%*******************************
%
% ARGUMENTS
%  zred         Redshift.
%
%
% RETURNS
%  DAng         Proper angular distance (Mpc).
%  h            H0/(100 km/s/Mpc)
%
% COMPATIBILITY: Matlab
%
% REQUIREMENTS: 
%	         cdenCosparamInit.m called previously
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  05 10 21 Creation date. (After SNRGetDLum.m)
%
function [DAng,h] = LyAGetDAng(zred);
[om_m,om_v,om_bh2,h,an,sigma8] = cdenCosparamInit;
z_min = 0;
nz = length(zred);
tol = 1.e-5;
for i = 1:nz
  z_max = zred(i);
%  yfnc = quad("SNRdyfncdz",z_min,z_max,tol); %for octave
  yfnc = integral(@SNRdyfncdz,z_min,z_max,'RelTol',tol);
  D(i) = (2.99792458e3/ h)*yfnc;
end
DAng = D./ (1 + zred);
