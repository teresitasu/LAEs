%
%*********************************************
% Phi = LyAGetQSOLF_KWH19(M1450,zred,imod
%*********************************************
%*********************************************
%
%
% Returns QSO luminosity function
% Uses Kulkarni et al. (2019) MNRAS 488:1035
%
% ARGUMENTS
% M1450  Absolute magnitude at 1450A
% zred     List of redshifts (CURRENTLY FOR A SINGLE REDSHIFT)
% imod    1, 2 or 3
%
% RETURNS
%  Phi       2D array containing Phi(M1450,z) (row,col)
%             (Phi(M1450,z)dM1450 = comoving number density in Mpc^-3 of QSOs at z
%	        with M1450 absolute magnitude between M1450 and M1450 + dM1450)
%  
%
% COMPATIBILITY: Matlab, Octave
%
% REQUIREMENTS
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  11 09 19 Creation date. (From LyAGetQSOLF.m)
%
%
function Phi = LyAGetQSOLF_KWH19(M1450,zred,imod)
M1450 = M1450(:);
lenM = length(M1450);
lenz = length(zred);
M1450mat = repmat(M1450,1,lenz);
Phi = zeros(lenM,lenz);
zp1 = 1 + zred;
% Define Chebyshev poloynomials of first kind, orders 0 to 3.
T0 = ones(1,lenz);
T1 = zp1;
T2 = 2*zp1.*zp1 - 1;
T3 = 4*zp1.*zp1.*zp1 - 3*zp1;
if(imod<1 || imod>3)
  disp('LyAGetQSOLF_KWH19: invalid imod');
  imod
  return
end
if(imod==1)
  c0 = [-7.798,1.128,-0.120];
  c1 = [-17.163,-5.512,0.593,-0.024];
  c2 = [-3.223,-0.258];
  c3 = [-2.312,0.559,3.773,141.884,-0.171];
  phis = 10.^(c0*[T0',T1',T2']');
  Ms = c1*[T0',T1',T2',T3']';
  Msmat = repmat(Ms,lenM,1);
  alpha = c2*[T0',T1']';
  zeta = log10(zp1/ (1 + c3(3)));
  beta = c3(1) + c3(2)./ (10.^(c3(4)*zeta) + 10.^(c3(5)*zeta));
  alphap1 = 1 + alpha;
  betap1 = 1 + beta;
  x = 10.^(0.4*(M1450mat - Msmat));
  phismat = repmat(phis,lenM,1);
  alphap1mat = repmat(alphap1,lenM,1);
  betap1mat = repmat(betap1,lenM,1);
  Phi = phismat./(x.^alphap1mat + x.^betap1mat);
end
if(imod==2)
  c0 = [-7.432,0.953,-0.112];
  c1 = [-15.412,-6.869,0.778,-0.032];
  c2 = [-2.959,-0.351];
  c3 = [-2.264,0.530,2.379,12.527,-0.229];
  phis = 10.^(c0*[T0',T1',T2']');
  Ms = c1*[T0',T1',T2',T3']';
  Msmat = repmat(Ms,lenM,1);
  alpha = c2*[T0',T1']';
  zeta = log10(zp1/ (1 + c3(3)));
  beta = c3(1) + c3(2)./ (10.^(c3(4)*zeta) + 10.^(c3(5)*zeta));
  alphap1 = 1 + alpha;
  betap1 = 1 + beta;
  x = 10.^(0.4*(M1450mat - Msmat));
  phismat = repmat(phis,lenM,1);
  alphap1mat = repmat(alphap1,lenM,1);
  betap1mat = repmat(betap1,lenM,1);
  Phi = phismat./(x.^alphap1mat + x.^betap1mat);
end
if(imod==3)
  c0 = [-6.942,0.629,-0.086];
  c1 = [-15.038,-7.046,0.772,-0.030];
  c2 = [-2.888,-0.383];
  c3 = [-1.602,-0.082];
  phis = 10.^(c0*[T0',T1',T2']');
  Ms = c1*[T0',T1',T2',T3']';
  Msmat = repmat(Ms,lenM,1);
  %Msmat
  %Ms
  alpha = c2*[T0',T1']';
  beta = c3(1) + c3(2)*zp1;
  alphap1 = 1 + alpha;
  betap1 = 1 + beta;
  x = 10.^(0.4*(M1450mat - Msmat));
  phismat = repmat(phis,lenM,1);
  alphap1mat = repmat(alphap1,lenM,1);
  betap1mat = repmat(betap1,lenM,1);
  Phi = phismat./(x.^alphap1mat + x.^betap1mat);
end

