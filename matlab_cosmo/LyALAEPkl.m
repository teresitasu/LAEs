%
%******************************************
%* [fk,Pk,PkLAE,PkLAEG,Pk_LAE_nsn,Pk_LAE] = LyALAEPkl(zred_out,b_LAE,b_delta,b_Gamma,tau_eff,lorder);
%******************************************
%******************************************
%
% Returns 3D power spectrum of LyA emitters, including correlations in
% metagalactic photoionization rate. Allows for QSO bias from Laurent et al.
% (2017) and galaxy bias.
%
% ARGUMENTS
% zred_out      Output redshifts used for dGammakCorr arrays (low to high)
% b_LAE         Lya emitter bias
% b_delta       Gas density bias factor (including tau_eff)
% b_Gamma       Ionization rate bias factor (from tau_eff)
% tau_eff       Median tau_eff
% lorder        Legendre order l
%
% RETURNS
% fk            Wavenumber (comoving h/ Mpc)
% Pk            Dark matter power spectrum (zred,fk)
% PkLAE         Lya emitter power spectrum from density alone
% PkLAEG        Lya emitter power spectrum from density and
%                        density-Gamma cross term
% Pk_LAE_nsn    Lya emitter power spectrum without shot noise
% Pk_LAE        Lya emitter power spectrum with shot noise
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
%  27 09 21 Creation date. (Adapted from LyATransmissionPkl.m.)
%
function [fk,Pk,PkLAE,PkLAEG,Pk_LAE_nsn,Pk_LAE] = LyALAEPkl(zred_out,b_LAE,b_delta,b_Gamma,tau_eff,lorder);
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
%if(exist('LyASolvedGammakCorr.mat')==2)
%  disp('using existing LyASolvedGammakCorr.mat file');
%  load('LyASolvedGammakCorr.mat');
  
if(exist('LyASolvedGammakCorr.txt')==2)
  %disp('using existing LyASolvedGammakCorr.mat file');
  %load('LyASolvedGammakCorr.mat');
  data_t = load('LyASolvedGammakCorr.txt');
  fk = transpose(data_t(:,1)); % co lumn 1 of the data text file is assigned the variable x
  dGammakCorr_nsn = transpose(data_t(:,2)); % column 2 is assigned the variable y
  dGammakCorr = transpose(data_t(:,3));
  
  data_t = load('cdenPowsp_omm_0_3_omv_0_7_ombh2_0_0_h_0_7_an_1_0_s8_0_8_ips_3.out');
  PS = transpose(data_t(:,2)); % column 1 of the data text file is assigned the variable x

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
PkLAE = repmat((b_LAE - tau_eff.*b_delta)',1,lenk).*deltk;
fprintf('Pk_LAE before %d \n',size(PkLAE)');
PkLAE = PkLAE.*PkLAE;
fprintf('Pk_LAE after %d \n',size(PkLAE)');

PkLAEG = -2*repmat(((b_LAE - tau_eff.*b_delta).*(tau_eff.*b_Gamma))',1,lenk).*ddGammak_nsn;
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
  PkLAEGl = 0;
  b_Gamma = 0;
end
Pk_LAE_nsn = PkLAEl + PkLAEGl + repmat(((tau_eff.*tau_eff).*b_Gamma.*b_Gamma)',1,lenk).*dGammakCorr_nsn;
Pk_LAE = PkLAEl + PkLAEGl + repmat(((tau_eff.*tau_eff).*b_Gamma.*b_Gamma)',1,lenk).*dGammakCorr;

fprintf('Pk_LAE_nsn %d \n',size(Pk_LAE_nsn)');
fprintf('Pk_LAE %d \n',size(Pk_LAE)');
fprintf('dGammakCorr %d \n',size(dGammakCorr)');

