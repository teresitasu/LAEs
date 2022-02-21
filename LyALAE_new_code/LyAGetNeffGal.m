%
%********************************************
% [neff,avgL,avgL2] = LyAGetNeffGal(zred,lmin)
%********************************************
%********************************************
%
%
% Computes effective comoving number density of galaxies from
% galaxy luminosity function of Bouwens et al. (2015) (ApJ 803:34)
%
% ARGUMENTS
% zred     List of redshifts
% lmin     Minimum Lum/Lums for integration
%
% RETURNS
%  neff     effective number density of galaxies as function or redshift (cMpc^-3)
%  avgL    Mean emissivity
%  avgL2  Mean L^2-emissivity
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
%  22 12 17 Creation date.
%
function [neff,avgL,avgL2] = LyAGetNeffGal(zred,lmin)
e = exp(1);
lenz = length(zred);
zp1 = 1 + zred;
LumLs = 1.80e28;
alpha = -1.87 + 0.1*(6 - zred);
alpha = max(alpha,-1.9);
g2 = gamma(alpha+2).*gammainc(lmin,alpha+2,'upper');
g3 = gamma(alpha+3).*gammainc(lmin,alpha+3,'upper');
coeff1 = 2.5*log10(e)*g2.*LumLs;
coeff2 = 2.5*log10(e)*g3.*LumLs.*LumLs;
phis = 0.44e-3*10.^(0.28*(6 - zred));
fesc = 1.8e-4*zp1.^3.4;
avgL = coeff1.*fesc.*phis;
avgL2 = coeff2.*fesc.*fesc.*phis;
neff = (avgL.*avgL).*(1./ avgL2);
neff = neff';
