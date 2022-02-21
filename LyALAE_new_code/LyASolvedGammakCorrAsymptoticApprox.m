%
%******************************************
%*  [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrAsymptoticApprox(Gammai,aj,aS,ab,om_m,h,zi,zf,k,Pk,bfrac,bet,cfrac,iQmod,bj,bnH,baA,bLLS,lgQbolmin,lgQbolmax,tau_Q,bG,tau_G)  *
%******************************************
%******************************************
%
% Returns deltaGammaCorr based on asympytotic limits alone: power spectrum of
% photoionization rate Gamma fluctuations, proportional to dark matter density
% power spectrum. Adopts conventional units of the dark matter power spectrum
% as comoving volume units (Mpc/h)^3. Allows for QSO bias from
% Laurent et al. (2017).
%
% ARGUMENTS
% Gammai      Initial guess on ionization rate (for LyASolveEmissMG.m)
% aj          Source comoving emissivity spectral power-law exponent
% aS          Source comoving emissivity evolution power-law exponent
% ab          Value of a_b for bias factor evolution bj ~ (1+z)^ab;
%                compute internally for ab < 0; set ab = 1 for high redshift (large aeff_LLS)
% om_m        Omega_m for all matter
% h           H0/ (100 km/s/Mpc)
% zi          Ionization redshift
% zf          Final redshift
% k           Single comoving wavenumber (h/ Mpc)
% Pk          Value of power spectrum P(k) at k at z = 0 (comoving Mpc^3/ h^3)
% bfrac       Fraction of baryons in diffuse component of IGM
% bet         Exponent for power law HI distribution of Lyman Limits Systems
% cfrac       Fraction of intergalactic attenuation due to Lyman Limit Systems
% iQmod       QSO evolution model
% bj          Source bias (set to <=0 to override bQ and bG and use -bj instead)
% bnH         Hydrogen density bias
% baA         Radiative recombination rate bias
% bLLS        Bias factor for LLS contribution to attenuation
% tau_Q       Lifetime of QSO sources (Myr)
% lgQbolmin   Log10 minimum QSO bolometric luminosity
% lgQbolmax   Log10 maximum QSO bolometric luminosity
% bG          Galaxy bias factor
% tau_G       Lifetime of galaxy sources (Myr)							%
%
% RETURNS
%  e24            Comoving emissivity coefficient (1e24 erg/s/Hz/Mpc^3)
%  zred           Array of redshifts from zf to zi
%  dGammaCorr     Power spectrum of fluctuation in photoionization rate Gamma from zf to zi
%  dGammaCorr_nsn Power spectrum of fluctuation in photoionization rate Gamma from zf to zi without shotnoise
%  dGammaCorr_sn  Power spectrum of fluctuation in photoionization rate Gamma from zf to zi for shotnoise term alone
%  Gamma          Mean metagalactic photoionization rate
%  aeff_d         Array of effective diffuse attenuation coefficient, matching zred array
%  aeff_LLS       Array of effective LLS attenuation coefficient, matching zred array
%  S              Source q*bj - (2*bnH + baA) array, matching zred array
%
% COMPATIBILITY: Matlab, Octave
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  04 06 18 Creation date. (Adapted from LyASolvedGammakCorrSS.m.)
%  05 07 18 Add lgQbolmin and lgQbolmax arguments
%  15 08 18 Modify for phi .ne. 1
%  24 08 18 Add factor convcomov to convert 1/ neff to comoving (Mpc/ h)^3
%  04 09 18 Add option for Lorentzian power spectrum approximation
%  21 09 18 Fix faulty denom_sn; allow for evolution in comoving neff through alpha_n
%  19 11 21 Fix asymptotic values to steady-state case for high k (for ab = 1) and low k
%  02 12 21 Modified to correct galaxy contribution to match required e24
%%function [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrAsymptoticApprox(Gammai,aj,aS,ab,om_m,h,zi,zf,k,Pk,bfrac,bet,cfrac,iQmod,bj,bnH,baA,bLLS,tau_Q,lgQbolmin,lgQbolmax,bG,tau_G);
function [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrAsymptoticApprox(Gammai,aj,aS,ab,om_m,h,zi,zf,k,Pk,bfrac,bet,cfrac,iQmod,bj,bnH,baA,bLLS,tau_Q,M1450min,M1450max,bG,tau_G);
mpc = 3.0856e24;
%alpha_n = 5; %comoving source density evolves like (1+z)^(-alpha_n); 5 ok at z=2.
h_LF = 0.7; %h = 0.7 used for luminosity functions
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
  disp('no EmissMG.mat available');
%  disp('generating new EmissMG.mat');
%  zred_MG = linspace(zf,zi,lzred_MG);
%  zp1_MG = 1 + zred_MG;
%  Gamma_i = Gammai*ones(1,lzred_MG);
%  [e24,Gamma_MG,xHI,aeffd_MG,aeffLLS_MG,otsfac_MG,err] = LyASolveEmissMG(Gamma_i,aj,aS,om_m,om_v,h,bfrac,bet,cfrac,zred_MG);
% save('EmissMG.mat','e24','Gamma_MG','xHI','aeffd_MG','aeffLLS_MG','otsfac_MG','zred_MG','err');
end
aeff_MG = aeffd_MG + aeffLLS_MG;
zred = zred_MG;
lzred = length(zred);
zp1 = 1 + zred;
zp1_3 = zp1.*zp1.*zp1;
%factor to convert shot noise from 1/ neff to comoving (Mpc/h)^3:
convcomov = h_LF*h_LF*h_LF;
%convcomov = 1; %USE THIS FOR NOW
aexp_flip = 1./ fliplr(zp1);
gg_flip = cdenGrowth(aexp_flip);
gg = fliplr(gg_flip);
deltk = gg*deltk_0/ gg_0;
bQ = 0.278*zp1.*zp1 + 0.57; %from BOSS
%bQ = 8./zp1.^0.5; %TEST CASE
%bQ = zp1; %TEST CASE
%bQ = 0.25.*zp1.*zp1; %TEST CASE
%bQ = 4; %TEST CASE
Gamma = Gamma_MG;
aeff = aeff_MG;
aeff_d = aeffd_MG;
aeff_LLS = aeffLLS_MG;
aeff_wtot = aeff_d - bGamLLS*aeff_LLS;
otsfac = otsfac_MG;
Hdc = h*(om_m*zp1.*zp1.*zp1+om_v).^0.5/(2.99792458e3*mpc);
aTot = aeff + Hdc.*otsfac;
%[neff_Q,avgL_Q,avgL2_Q] = LyAGetNeffQSO(zred,lgQbolmin,lgQbolmax,iQmod);
[neff_Q,avgL_Q,avgL2_Q] = LyAGetNeffQSO_KWH19(zred,M1450min,M1450max,iQmod);
eL24 = e24./ zp1.^aS;
if(bG > 0)
  [neff_G,avgL_G,avgL2_G] = LyAGetNeffGal(zred,lmin);
  fixfescG = (1e24*eL24 - avgL_Q)./ avgL_G;
  avgL_G = fixfescG.*avgL_G;
  avgL2_G = fixfescG.*fixfescG.*avgL2_G;
else
  fixQ = 1e24*eL24./ avgL_Q;
  avgL_Q = fixQ.*avgL_Q;
  avgL2_Q = fixQ.*fixQ.*avgL2_Q;
end
avgL = eL24*1e24;
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
  shot_G = zeros(lzred);
  avgavgL2_G = zeros(lzred);
end
jdf = q./ (Gamma*mpc);
if((cfrac < 1) || (abs(bGamLLS)>0))
   q = jdf./ aeff_wtot;
   S = bj.*q - (2*bnH + baA)*(aeff_d./ aeff_wtot) - bLLS*(aeff_LLS./ aeff_wtot);
else
   q = jdf./ aeff;
   S = bj.*q - bLLS*(aeff_LLS./ aeff);
end
S = S.*deltk;
if(ab < 0)
  ab = gradient(log(bj))./ gradient(log(zp1));
end
chiH = aeff./ Hdc;
chiH_d = aeff_d./ Hdc;
chiH_LLS = aeff_LLS./ Hdc;
%chiHb = chiH + otsfac; %otsfac is beta_H
chiHb = jdf./ Hdc; %phi*(chi+zeta)
chiHdL = chiH_d - bGamLLS*chiH_LLS;
gamH = chiHb - ab + 0.5;
kap = k*h./ (mpc*Hdc);
%zp1i = zp1(lzred);
%xi = zp1./ zp1i;
% kap -> 0 limit
denom = bGamLLS*chiH + chiHb + 1 - ab; %otsfac is beta_H
pi2 = pi/ 2;
%Use Lorentzian formulation for non-shotnoise signal
kaps = (pi2./ zp1).*denom;
fkap = -chiHdL./ (denom.*(1 + (kap./ kaps).^2).^0.5);
S = S.*fkap;
S2ns = S'*S;
% Add shotnoise term in E-deS limit
%%%%%%%%%%
%%Use Lorentzian formulation for shot noise
%Use extended Lorentzian formulation for shot noise
if(bG > 0)
  neff = (avgL.*avgL)./ (avgL2_Q + avgL2_G);
else
  neff = (avgL.*avgL)./ avgL2_Q;
end
alpha_n = -gradient(log(neff))./ gradient(log(zp1)); %comoving neff evolution exponent
HMyr = 3.15e20*h*(om_m*zp1.*zp1.*zp1+om_v).^0.5/ mpc;
tau_S = tau_Q; %assume tau_Q dominates
denom_sn = (1 - bet)*chiH + chiHb + 0.75 - 0.5*alpha_n;
% Use sliding scale-length
kaps_sn = 2.*denom_sn./ (HMyr*tau_S) - ((chiHb - 0.5)/(pi2*pi2) + (1 - bet)*chiH).^2;
kaps_sn = (pi2./ zp1).*kaps_sn.^0.5;
kaps_sn = kaps_sn./ (1 + pi2*((chiHb - 0.5)/(pi2*pi2) + (1 - bet)*chiH)./ (kap.*zp1));
AA = chiHb.*chiHb;
AA = AA./ (2.*denom_sn./ (HMyr*tau_S) - ((chiHb - 0.5)/(pi2*pi2)+(1 - bet)*chiH).^2);
fkap_sn = (AA./ (1 + (kap./ kaps_sn).^2)).^0.5;
S_sn = (convcomov.^0.5).*fkap_sn./ neff.^0.5;
shot = S_sn'*S_sn;
S2 = S2ns + shot;
dGammaCorr = S2;
dGammaCorr_nsn = S2ns;
dGammaCorr_sn = shot;
