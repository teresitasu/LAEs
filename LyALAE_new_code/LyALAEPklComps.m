%
%******************************************
%* [fk,Pk,PkLAEl,PkLAEGl,Pk_LAE_nsn,Pk_LAE] = LyALAEPklComps(zred_out,b_LAE,b_delta,b_Gamma,a_LAE,f_LAE,f2_LAE,tau_eff,lorder);
%******************************************
%******************************************
%
% Returns 3D power spectrum of LyA emitters, including correlations in
% metagalactic photoionization rate. Allows for QSO bias from Laurent et al. (2017)
% and galaxy bias.
%
% ARGUMENTS
% zred_out    Output redshifts used for dGammakCorr arrays (low to high)
% b_LAE       LyA emitter bias factor
% b_delta     Gas density bias factor of tau_eff
% b_Gamma     Ionization rate bias factor of tau_eff
% a_LAE       LAE ew power law: dN/dw ~ w^a_LAE (see Ouchi et al. 2020)
% f_LAE       Overlap integral for nu-nu_alpha < 0 in emission line profile
% f2_LAE      Overlap integral sqr for nu-nu_alpha < 0 in emission line profile
% tau_eff     Median tau_eff of IGM
% lorder      Legendre order l
%
% RETURNS
% fk          Wavenumber (h/ Mpc)
% Pk          Dark matter power spectrum (zred,fk)
% PkLAEl      Legendre component of LAE power spectrum from density alone
% PkLAEGl     Legendre component of LAE flux power spectrum from density and
%                        density-Gamma cross term
% Pk_LAE_nsn  LyA emitter power spectrum without shot noise
% Pk_LAE      LyA emitter power spectrum with shot noise
%
% COMPATIBILITY: Matlab, Octave
%
% REQUIREMENTS:
%	         cdenCosparamInit.m called previously
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  01 10 21 Creation date. (Adapted from LyATransmissionPklComps.m.)
%  13 12 21 Modified to allow for a_LAE, f_LAE and f2_LAE
%				%
function [fk,Pk,PkLAEl,PkLAEGl,Pk_LAE_nsn,Pk_LAE] = LyALAEPklComps(zred_out,b_LAE,b_delta,b_Gamma,a_LAE,f_LAE,f2_LAE,tau_eff,lorder);
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
PkLAEa = repmat(b_LAE',1,lenk).*deltk;
PkLAEb = repmat(((a_LAE+1)*f_LAE.*tau_eff.*b_delta)',1,lenk).*deltk;
PkLAEb2 = repmat(((a_LAE+1)*f2_LAE.*tau_eff.*b_delta)',1,lenk).*deltk;
PkLAE = PkLAEa.*PkLAEa + 2*PkLAEa.*PkLAEb + PkLAEb2.*PkLAEb2;
PkLAEG = 2*(a_LAE+1)*repmat(((b_LAE*f_LAE + (a_LAE+1)*f2_LAE.*tau_eff.*b_delta).*(tau_eff.*b_Gamma))',1,lenk).*ddGammak_nsn;
if(lorder==0)
  PkLAEl = (1 + betafac/ 1.5 + beta2fac/ 5.).*PkLAE;
  PkLAEGl = (1 + betafac/ 3.).*PkLAEG;
end
if(lorder==2)
  PkLAEl = 4*betafac.*(1./ 3 + betafac/ 7.).*PkLAE;
  PkLAEGl = (2*betafac/ 3.).*PkLAEG;
  b_Gamma = 0;
end
if(lorder==4)
  PkLAEl = 8*beta2fac.*PkLAE/ 35.;
  PkLAEGl = 0*PkLAEl;
  b_Gamma = 0;
end
Pk_LAE_nsn = PkLAEl + PkLAEGl + repmat(((a_LAE+1)^2*f2_LAE.*(tau_eff.*tau_eff).*(b_Gamma.*b_Gamma))',1,lenk).*dGammakCorr_nsn;
Pk_LAE = PkLAEl + PkLAEGl + repmat(((a_LAE+1)^2*f2_LAE.*(tau_eff.*tau_eff).*(b_Gamma.*b_Gamma))',1,lenk).*dGammakCorr;
