%
%******************************************
%* [fk,Pk,PkH,PkHG,Pk_alpha_nsn,Pk_alpha] = LyATransmissionPkl(zred_out,b_delta,b_Gamma,lorder);
%******************************************
%******************************************
%
% Returns 3D power spectrum of LyA transmission, including correlations in
% metagalactic photoionization rate. Allows for QSO bias from Laurent et al. (2017)
% and galaxy bias.
%
% ARGUMENTS
% zred_out    Output redshifts used for dGammakCorr arrays (low to high)
% b_delta       Gas density bias factor
% b_Gamma   Ionization rate bias factor
% lorder        Legendre order l
%
% RETURNS
% fk               Wavenumber (h/ Mpc)
% Pk              Dark matter power spectrum (zred,fk)
% PkH                 LyA transmission fluctuation power spectrum from density alone
% PkHG               LyA transmission fluctuation power spectrum from density and
%                        density-Gamma cross term
% Pk_alpha_nsn   LyA transmission fluctuation power spectrum without shot noise
% Pk_alpha         LyA transmission fluctuation power spectrum with shot noise
%
% COMPATIBILITY: Matlab(?), Octave
%
% REQUIREMENTS:
%	         cdenCosparamInit.m called previously
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  24 04 19 Creation date. (Adapted from LyASolvedGammakCorr.m.)
%
function [fk,Pk,PkH,PkHG,Pk_alpha_nsn,Pk_alpha] = LyATransmissionPkl(zred_out,b_delta,b_Gamma,lorder);
global omega; %from. eg, cdenCosparamInit.m
om_m = omega(1);
om_v = omega(2);
om_k = omega(3);
lok = 0;
if(lorder==0)
  lok = 1;
end
if(lorder==2)
  lok = 1;
end
if(lorder==4)
  lok = 1;
end
if(lok==0)
  disp('must use lorder = 0, 2 or 4');
  return;
end
if(exist('LyASolvedGammakCorr.mat')==2)
  disp('using existing LyASolvedGammakCorr.mat file');
  load('LyASolvedGammakCorr.mat');
else
  disp('no LyASolvedGammakCorr.mat file');
  return;
end
lenk  = length(fk);
dGammak_nsn = dGammakCorr_nsn.^0.5;
deltk_0 = PS.^0.5;
gg_0 = cdenGrowth(1.);
zp1 = 1 + zred_out;
zp1_3 = zp1.*zp1.*zp1;
Hfac = (om_m*zp1_3 + om_v + om_k*zp1.*zp1).^0.5;
aexp_flip = 1./ fliplr(zp1);
[gg_flip,gv_flip] = cdenGrowth(aexp_flip);
gg = fliplr(gg_flip);
gv = fliplr(gv_flip);
deltk = gg'*deltk_0/ gg_0;
Pk = deltk.*deltk;
ddGammak_nsn = deltk.*dGammak_nsn;
beta = gv./ Hfac;
beta2 = beta.*beta;
betafac = repmat(beta',1,lenk);
beta2fac = repmat(beta2',1,lenk);
PkH = b_delta.*deltk;
PkH = PkH.*PkH;
PkHG = 2*b_delta*b_Gamma.*ddGammak_nsn;
if(lorder==0)
  PkHl = (1 + betafac/ 1.5 + beta2fac/ 5.).*PkH;
  PkHGl = (1 + betafac/ 3.).*PkHG;
end
if(lorder==2)
  PkHl = 4*betafac.*(1./ 3 + betafac/ 7.).*PkH;
  PkHGl = (2*betafac/ 3.).*PkHG;
  b_Gamma = 0;
end
if(lorder==4)
  PkHl = 8*beta2fac.*PkH/ 35.;
  PkHGl = 0;
  b_Gamma = 0;
end
Pk_alpha_nsn = PkHl + PkHGl + b_Gamma*b_Gamma*dGammakCorr_nsn;
Pk_alpha = PkHl + PkHGl + b_Gamma*b_Gamma*dGammakCorr;
