%
%******************************************
%* [fk,Pk,PkLAEl,PkLAEGl,Pk_LAE_nsn,Pk_LAE] = LyALAEPklComps(zred_out,b_LAE,b_delta,b_Gamma,tau_eff,lorder);
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
%
function [fk,Pk,PkLAEl,PkLAEGl,Pk_LAE_nsn,Pk_LAE] = LyALAEPklComps(zred_out,b_LAE,b_delta,b_Gamma,tau_eff,lorder);
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
if(exist('LyASolvedGammaCorrAAL_W14B21_bG3_tG100_tQ1_z6.0.out')==2)
  %disp('using existing LyASolvedGammakCorrSS_z6.out file');
  %load('LyASolvedGammakCorrSS_z6.out');
  data_t = load('LyASolvedGammaCorrAAL_W14B21_bG3_tG100_tQ1_z6.0.out');
  fk = transpose(data_t(:,1)); % co lumn 1 of the data text file is assigned the variable x
  dGammakCorr_nsn = transpose(data_t(:,3)); % column 2 is assigned the variable y
  dGammakCorr = transpose(data_t(:,4));
  PS = transpose(data_t(:,5));
  
  %data_t = load('cdenPowsp_omm_0_3_omv_0_7_ombh2_0_0_h_0_7_an_1_0_s8_0_8_ips_3.out');
  %PS = transpose(data_t(:,2)); % column 1 of the data text file is assigned the variable x

else
  disp('no LyASolvedGammakCorrSS_z6.out file');
  return;
end

lenk  = length(fk);
dGammak_nsn = dGammakCorr_nsn.^0.5;
deltk_0 = PS.^0.5;
gg_0 = cdenGrowth(1.);
fprintf('gg_0 %f \n',gg_0');
zp1 = 1 + zred_out;
zp1_3 = zp1.*zp1.*zp1;
Hfac = (om_m*zp1_3 + om_v + om_k*zp1.*zp1).^0.5;
aexp_flip = 1./ fliplr(zp1);
fprintf('aexp_flip %f \n',aexp_flip');
[gg_flip,gv_flip] = cdenGrowth(aexp_flip);
gg = fliplr(gg_flip);
gv = fliplr(gv_flip);
fprintf('gg %f, gv %f %f\n',gg,gv);
deltk = gg'*deltk_0/ gg_0;
Pk = deltk.*deltk;
ddGammak_nsn = deltk.*dGammak_nsn;
beta = gv./ Hfac;
beta2 = beta.*beta;
betafac = repmat(beta',1,lenk);
beta2fac = repmat(beta2',1,lenk);
PkLAE = repmat((b_LAE - tau_eff.*b_delta)',1,lenk).*deltk;
PkLAE = PkLAE.*PkLAE;
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
  PkLAEGl = 0*PkLAEl;
  b_Gamma = 0;
end
Pk_LAE_nsn = PkLAEl + PkLAEGl + repmat(((tau_eff.*tau_eff).*(b_Gamma.*b_Gamma))',1,lenk).*dGammakCorr_nsn;
Pk_LAE = PkLAEl + PkLAEGl + repmat(((tau_eff.*tau_eff).*(b_Gamma.*b_Gamma))',1,lenk).*dGammakCorr;
