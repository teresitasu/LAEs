%
%******************************************
%*  [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrSS(Gammai,aj,aS,om_m,hubb,zi,zf,k,Pk,lzred_MG,bfrac,bet,cfrac,iQmod,lgQbolmin/M1450min,lgQbolmax/M1450max,bj,bnH,baA,bLLS,bG)  *
%******************************************
%******************************************
%
% Returns steady-state power-spectrum in fluctuations in photoionization rate
 % Gamma, in comoving units (Mpc/ h)^3. Allows for QSO bias from Laurent et al. (2017) and galaxy bias.
%
% ARGUMENTS
% Gammai      Initial guess on ionization rate (for LyASolveEmissMG.m)
% aj              Source comoving emissivity spectral power-law exponent
% aS             Source comoving emissivity evolution power-law exponent
% om_m        Omega_m for all matter
% hubb          H0/ (100 km/s/Mpc)
% zi              Ionization redshift
% zf              Final redshift
% k               Single comoving wavenumber (h/ Mpc)
% Pk             Value of power spectrum P(k) at k at z = 0 (comoving (Mpc/h)^3)
% lzred_MG   Length of redshift array (when required)
% bfrac         Fraction of baryons in diffuse component of IGM
% bet            Exponent for power law HI distribution of Lyman Limits Systems
% cfrac          Fraction of intergalactic attenuation due to Lyman Limit Systems
% iQmod       QSO evolution model
%% lgQbolmin   Log10 minimum QSO bolometric luminosity
%% lgQbolmax   Log10 maximum QSO bolometric luminosity
% M1450min     Minimum QSO absolute magnitude at 1450A
% M1450max     Maximum QSO absolute magnitude at 1450A
% bj              Source bias (set to <=0 to override bQ and bG and use -bj instead)
% bnH           Hydrogen density bias
% baA           Radiative recombination rate bias
% bLLS          Bias factor for LLS contribution to attenuation
% bG             Galaxy bias factor
%
% RETURNS
%  e24            Comoving emissivity coefficient (1e24 erg/s/Hz/Mpc^3)
%  zred           Array of redshifts from zf to zi
  %  dGammaCorr Fluctuation power spectrum of photoionization rate Gamma from zf to zi (units comoving (Mpc/h)^3)
%  dGammaCorr_nsn Fluctuation power spectrum of photoionization rate Gamma from zf to zi without shotnoise (units comoving (Mpc/h)^3)
%  dGammaCorr_sn Fluctuation power spectrum in photoionization rate Gamma from zf to zi for shotnoise term alone (units comoving (Mpc/h)^3)
%  Gamma     Mean metagalactic photoionization rate
%  aeff_d       Array of effective diffuse attenuation coefficient, matching zred array
%  aeff_LLS    Array of effective LLS attenuation coefficient, matching zred array
%  S              Source q*bj - (2*bnH + baA) array, matching zred array
%
% COMPATIBILITY: Matlab(?), Octave
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  16 02 18 Creation date. (Adapted from LyASolvedGammakCorr.m.)
%  16 03 18 Fix missing hubb factor in x
%  05 07 18 Add lgQbolmin and lgQbolmax arguments
%  01 08 18 Modify for phi .ne. 1
%  24 08 18 Add factor convcomov to convert 1/ neff to comoving (Mpc/ h)^3
%  31 08 21 Adapt to allow for Kulkarni et al. (2019) QSO LF models
%
%function [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrSS(Gammai,aj,aS,om_m,hubb,zi,zf,k,Pk,lzred_MG,bfrac,bet,cfrac,iQmod,lgQbolmin,lgQbolmax,bj,bnH,baA,bLLS,bG);
function [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrSS(Gammai,aj,aS,om_m,hubb,zi,zf,k,Pk,lzred_MG,bfrac,bet,cfrac,iQmod,M1450min,M1450max,bj,bnH,baA,bLLS,bG);
mpc = 3.0856e24;
lmin = 0.01;
%bGamLLS = 0;
bGamLLS = 1 - bet;
om_v = 1 - om_m; %assume flat universe
deltk_0 = Pk^0.5;
gg_0 = cdenGrowth(1.);
% Get ionization rate
if(exist('EmissMG.mat')==2)
  disp('using existing EmissMG.mat');
  load('EmissMG.mat');
  zp1_MG = 1 + zred_MG;
else
  disp('generating new EmissMG.mat');
  zred_MG = linspace(zf,zi,lzred_MG);
  zp1_MG = 1 + zred_MG;
  Gamma_i = Gammai*ones(1,lzred_MG);
  [e24,Gamma_MG,xHI,aeffd_MG,aeffLLS_MG,otsfac_MG,err] = LyASolveEmissMG(Gamma_i,aj,aS,om_m,om_v,hubb,bfrac,bet,cfrac,zred_MG);
save('EmissMG.mat','e24','Gamma_MG','xHI','aeffd_MG','aeffLLS_MG','otsfac_MG','zred_MG','err');
end
aeff_MG = aeffd_MG + aeffLLS_MG;
zred = zred_MG;
lzred = length(zred);
zp1 = 1 + zred;
zp1_3 = zp1.*zp1.*zp1;
%factor to convert shot noise from 1/ neff to comoving (Mpc/h)^3:
convcomov = hubb*hubb*hubb*zp1_3;
aexp_flip = 1./ fliplr(zp1);
gg_flip = cdenGrowth(aexp_flip);
gg = fliplr(gg_flip);
deltk = gg*deltk_0/ gg_0;
bQ = 0.278*zp1.*zp1 + 0.57;
%bQ = 8./zp1.^0.5; %TEST CASE
%bQ = zp1; %TEST CASE
%bQ = 0.25.*zp1.*zp1; %TEST CASE
%bQ = 4; %TEST CASE
Gamma = Gamma_MG;
aeff = aeff_MG;
aeff_d = aeffd_MG;
aeff_LLS = aeffLLS_MG;
otsfac = otsfac_MG;
Hdc = hubb*(om_m*zp1.*zp1.*zp1+om_v).^0.5/(2.99792458e3*mpc);
%[neff_Q,avgL_Q,avgL2_Q] = LyAGetNeffQSO(zred,lgQbolmin,lgQbolmax,iQmod);
[neff_Q,avgL_Q,avgL2_Q] = LyAGetNeffQSO_KWH19(zred,M1450min,M1450max,iQmod);
if(bG > 0)
  [neff_G,avgL_G,avgL2_G] = LyAGetNeffGal(zred,lmin);
  avgL = avgL_Q + avgL_G;
else
  avgL = avgL_Q;
end
if(bj <= 0)
  bj = -bj;
else
  if(bG > 0)
    bj = (bQ.*avgL_Q + bG.*avgL_G)./ avgL;
  else
    bj = bQ;
  end
end
q = 0.999e-16*e24*zp1.^(3 - aS)/ (3 + aj);
shot = zeros(lzred);
shot_Q = zeros(lzred);
avgavgL2_Q = zeros(lzred);
if(bG > 0)
  shot_G = zeros(lzred,lzred);
  avgavgL2_G = zeros(lzred);
end
Delta_t0 = 6523/ (hubb*om_m^0.5);
jdf = q./ (Gamma*mpc);
aTot = jdf;
q = jdf./ aTot;
S = bj.*q - (2*bnH + baA)*(aeff_d./ aTot) - bLLS*(aeff_LLS./ aTot);
S = S.*deltk;
S2ns = S.*S;
% Add shotnoise terms
avgavgL2_Q = zp1_3.*avgL2_Q;
shot_Q = avgavgL2_Q./ (zp1_3.*avgL.*zp1_3.*avgL);
shot_Q = shot_Q.*q.*q;
if(bG > 0)
    avgavgL2_G =zp1_3.*avgL2_G;
    shot_G = avgavgL2_G./ (zp1_3.*avgL.*zp1_3.*avgL);
    shot_G = shot_G.*q.*q;
end
% include conversion to comoving (Mpc/h)^3 through factor convcomov
if(bG> 0)
    shot = convcomov.*(shot_Q + shot_G);
else
    shot = convcomov.*shot_Q;
end
S2 = S2ns + shot;
x = zp1.*k*hubb./ (mpc*aTot);
denom = atan(x)./ x;
denom = 1./ denom - (aeff_d - bGamLLS.*aeff_LLS)./ aTot;
denom = denom.*denom;
dGammaCorr = S2./ denom;
dGammaCorr_nsn = S2ns./ denom;
dGammaCorr_sn = shot./ denom;
