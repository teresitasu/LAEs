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
convcomov = h_LF*h_LF*h_LF*zp1_3;
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
%%%%
%% High kap limit
%fackap0 = chiHdL.*(1 - xi.^denom)./ denom;
%ui= 2.*zp1.*kap.*(1 - xi.^0.5);
%fackap1 = pi2*chiHdL./ (zp1.*kap);
%fackap21 = pi2*chiHdL;
%fackap22 =  -(chiHb - ab + 0.5)/ pi2;
%fackap23 = -(xi.^gamH./ (pi*(1 - xi.^0.5))).*cos(ui);
%fackap2 = (fackap21 + fackap22 + fackap23)./ (zp1.*kap);
%fackap = fackap1.*(1 + fackap2);
%%%%
%fkap = 1./ (1./ fackap0 + 1./ fackap); %produces large glitch at transition
%%%%%%%%%%
%%Use break formulation
%maskfk0 = find(kap<2);
%maskfkp = find(kap>=2);
%fkap(maskfk0) = fackap0(maskfk0);
%fkap(maskfkp) = fackap(maskfkp);
%%%%%%%%%%
%Use Lorentzian formulation for non-shotnoise signal
kaps = (pi2./ zp1).*denom;
fkap = -chiHdL./ (denom.*(1 + (kap./ kaps).^2).^0.5);
%%%%%%%%%%
S = S.*fkap;
S2ns = S'*S;
% Add shotnoise term in E-deS limit
%%%%%%%%%%
%%Use Lorentzian formulation for shot noise
%Use extended Lorentzian formulation for shot noise
if(bG > 0)
  neff = (avgL.*avgL.*zp1_3)./ (avgL2_Q + avgL2_G);
else
  neff = (avgL.*avgL.*zp1_3)./ avgL2_Q;
end
alpha_n = -gradient(log(neff))./ gradient(log(zp1));
HMyr = 3.15e20*h*(om_m*zp1.*zp1.*zp1+om_v).^0.5/ mpc;
tau_S = tau_Q; %assume tau_Q dominates
denom_sn = (1 - bet)*chiH + chiHb + 0.75 - 0.5*alpha_n;
%kaps_sn1 = (pi*pi./ (4*zp1)).*(denom_sn./ gamH)./(HMyr*tau_S*(bet - 1).*chiH);
%kaps_sn2 = (pi./ zp1).*(denom_sn./ (2*HMyr*tau_S)).^0.5;
%fkap_sn = (0.5.^0.5)*chiHb./ (1 + (kap./ kaps_sn2).^2).^0.5;
%fkap_sn = (0.5.^0.5)*chiHb./ (1 - kap./ kaps_sn1 + (kap./ kaps_sn2).^2).^0.5;
% add empirical fudge depending on 2.5*kaps
%fkap_sn = (0.5.^0.5)*chiHb./ (1 + log(1+0.4*kap./ kaps) - kap./ kaps_sn1 + (kap./ kaps_sn2).^2).^0.5;
%S_sn = (convcomov.^0.5).*fkap_sn.*(HMyr*tau_S./ (denom_sn.*neff)).^0.5;
% Use sliding scale-length
kaps_sn = 2.*denom_sn./ (HMyr*tau_S) - ((chiHb - 0.5)/(pi2*pi2) + (1 - bet)*chiH).^2;
kaps_sn = (pi2./ zp1).*kaps_sn.^0.5;
kaps_sn = kaps_sn./ (1 + pi2*((chiHb - 0.5)/(pi2*pi2) + (1 - bet)*chiH)./ (kap.*zp1));
AA = chiHb.*chiHb;
AA = AA./ (2.*denom_sn./ (HMyr*tau_S) - ((chiHb - 0.5)/(pi2*pi2)+(1 - bet)*chiH).^2);
fkap_sn = (AA./ (1 + (kap./ kaps_sn).^2)).^0.5;
S_sn = (convcomov.^0.5).*fkap_sn./ neff.^0.5;
shot = S_sn'*S_sn;
%%%%%%%%%%
%%%%%%%%%%
% SN TERM NOT YET WORKING: NEEDS larger zred array including lzm factor?
%Delta_t0 = 6523/ (h*om_m^0.5);
%for iz = 1:lzred
%  Delta_t = Delta_t0*abs(1./ zp1.^1.5 - 1./ zp1(iz)^1.5);
%  mask_t = find(Delta_t<tau_Q);
%  avgavgL2_Q(mask_t) = 0.5*(zp1_3(iz)*avgL2_Q(iz) + zp1_3(mask_t).*avgL2_Q(mask_t));
%  shot_Q(iz,mask_t) = avgavgL2_Q(mask_t)./ (zp1_3(iz)*avgL(iz)*zp1_3(mask_t).*avgL(mask_t));
%  shot_Q(iz,mask_t) = shot_Q(iz,mask_t).*(1 - Delta_t(mask_t)/ tau_Q).*q(mask_t);
%  shot_Q(iz,mask_t) = convcomov(iz)*shot_Q(iz,mask_t)*q(iz);
%  clear mask_t;
%  if(bG > 0)
%    mask_t = find(Delta_t<tau_G);
%    avgavgL2_G(mask_t) = 0.5*(zp1_3(iz)*avgL2_G(iz) + zp1_3(mask_t).*avgL2_G(mask_t));
%    shot_G(iz,mask_t) = avgavgL2_G(mask_t)./ (zp1_3(iz)*avgL(iz)*zp1_3(mask_t).*avgL(mask_t));
%    shot_G(iz,mask_t) =shot_G(iz,mask_t).*(1 - Delta_t(mask_t)/ tau_G).*q(mask_t);
%    shot_G(iz,mask_t) = convcomov(iz)*shot_G(iz,mask_t)*q(iz);
%  end
%end
%if(bG> 0)
%  shot = shot_Q + shot_G;
%else
%  shot = shot_Q;
%end
%% Add shotnoise terms in infinite tau_Q and tau_G limit
%avgavgL2_Q = zp1_3.*avgL2_Q;
%shot_Q = avgavgL2_Q./ (zp1_3.*avgL.*zp1_3.*avgL);
%shot_Q = convcomov.*shot_Q.*q.*q;
%if(bG > 0)
%    avgavgL2_G =zp1_3.*avgL2_G;
%    shot_G = avgavgL2_G./ (zp1_3.*avgL.*zp1_3.*avgL);
%    shot_G = convcomov.*shot_G.*q.*q;
%end
%if(bG> 0)
%    shot = shot_Q + shot_G;
%else
%    shot = shot_Q;
%end
%diagshot = diag(shot)';
%S2 = S2ns + diagshot;
%%%%%%%%%%
S2 = S2ns + shot;
dGammaCorr = S2;
dGammaCorr_nsn = S2ns;
%dGammaCorr_sn = diagshot;
dGammaCorr_sn = shot;
