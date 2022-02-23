%
%********************************************************
% [g_LAE] = LyALAEComputeG_LAE(tau_eff,x0,sighrat,alpha)
%********************************************************
%********************************************************
%
%
% Computes g_LAE overlap integral for LyALAEPkl(Comps).m
%
% ARGUMENTS
% tau_eff  Effective LyA optical depth (may be array)
% x0       LAE profile peak location (units 2^0.5*sigma_emline) (single value)
% sighrat  Ratio of sigma_halo/ sigma_emline
%
% RETURNS
%  g_LAE   Overlap integral
%
% COMPATIBILITY: Matlab, [Octave]
%
% REQUIREMENTS
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  17 02 22 Creation date. (After LyALAEComputeAvgGG2_LAE.m)
%
function g_LAE = LyALAEComputeG_LAE(tau_eff,x0,sighrat,alpha);
lent = length(tau_eff);
g_LAE = zeros(1,lent);
ff = @(x) min(1e20,erfc(x)./ (1 + erf(x)));
for i = 1:lent
  fncd = @(x) exp(-(x-x0).^2/ sighrat^2).*((1+ff(x))./(1+ff(x)*exp(-tau_eff(i)))).^(alpha+1);
  fncn = @(x) fncd(x).*ff(x)*exp(-tau_eff(i))./(1+ff(x)*exp(-tau_eff(i)));
  gnorm_LAE = integral(fncd,-5+x0,5+x0);
  g_LAE(i) = integral(fncn,-5+x0,5+x0)/ gnorm_LAE;
end
