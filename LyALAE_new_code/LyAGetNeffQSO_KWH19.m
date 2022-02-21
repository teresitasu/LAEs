%
%***************************************************************************
% [neff,avgL,avgL2] = LyAGetNeffQSO_KWH19(zred,M1450min,M1450max,imod)
%***************************************************************************
%***************************************************************************
%
%
% Computes effective comoving number density of QSOs from
% QSO luminosity function
%
% ARGUMENTS
% zred           List of redshifts
% M1450min   Minimum 1450 absolute magnitude
% M1450max  Maximum 1450 absolute magnitude
% imod          1, 2 or 3
%
% RETURNS
%  neff     effective number density of QSOs as function or redshift (cMpc^-3)
%  avgL    Mean comoving emissivity (erg/s/Hz/cMpc^3)
%  avgL2  Mean comoving L^2-emissivity [(erg/s/Hz)^2/cMpc^3]
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
%  11 09 19 Creation date. (After LyAGetNeffQSO.m)
%
function [neff,avgL,avgL2] = LyAGetNeffQSO_KWH19(zred,M1450min,M1450max,imod)
lenz = length(zred);
lenM = 1001;
M1450 = linspace(M1450min,M1450max,lenM);
dM1450 = gradient(M1450);
Phi = LyAGetQSOLF_KWH19(M1450,zred,imod);
% Convert to luminosity at Lyman edge (KWH19 Eqs.(20) and (21))
LumL = ((912/1450)^0.61)*10.^(-0.4*(M1450 - 51.60));
LumLmat = repmat(LumL',1,lenz);
avgL = (LumLmat.*Phi)'*dM1450';
avgL2 = (LumLmat.*LumLmat.*Phi)'*dM1450';
neff = (avgL.*avgL).*(1./ avgL2);
avgL = avgL';
avgL2 = avgL2';
neff = neff';
