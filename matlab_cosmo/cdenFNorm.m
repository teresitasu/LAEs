%
%******************************************************
%*  cdenFNorm(lgfk,pnorm,Gamma,an,Tk,xf,ismooth,ips)  *
%******************************************************
%******************************************************
% Computes normalization integrand for power spectrum.
%
% ARGUMENTS
%  lgfk       Array of (logarithmically-spaced) wavenumber.
%  pnorm      Normalization of matter power spectrum.
%  Gamma      Curvature of matter power spectrum.
%  an         Spectral index of matter power spectrum.
%  Tk         Transfer function for matter power spectrum.
%  xf         Filtering scale (inverse units of dlfkp).
%  ismooth    Method for smoothing fluctuations.
%  ips        Method of computing matter power spectrum.
%
% COMPATIBILITY: Matlab(?), Octave
%
% REQUIREMENTS:
%         cdenPowspInit.m, cdenPowsp.m
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  19 04 16 Creation date. (After FNORM subfunction in lss/src/Delta2.f.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function FNORM = cdenFNorm(lgfk,pnorm,Gamma,an,Tk,xf,ismooth,ips);
fk = exp(lgfk);
%fk = lgfk; %TEST for linear integration
y = xf*fk;
if(ismooth == 0)
  W = 1.d0;
end
if(ismooth == 1)
  W = 3.*(sin(y) - y.*cos(y))./ y./ y./ y;
end
if(ismooth == 2)
  W = 0.;
  maskl = find(y<10);
  W(maskl) = 1./ exp(y(maskl).*y(maskl));
end
if(ismooth == 3)
  W = 1./ (1 + y.*y);
end
PS = cdenPowsp(fk,pnorm,Gamma,an,Tk,ips);
FNORM = fk.*fk.*fk.*PS.*W.*W/ (2*pi^2);
%FNORM = fk.*fk.*PS.*W.*W/ (2*pi^2); %TEST for linear integration
