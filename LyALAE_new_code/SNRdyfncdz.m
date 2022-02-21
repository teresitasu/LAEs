%
%*****************************************************
%*  dydz = SNRdyfncdz(z)                             *
%*****************************************************
%*****************************************************
%
% Integrand for SNRDLum to compute luminosity distance.
%
function dydz = SNRdyfncdz(z);
global omega;
Efnc = (omega(1)*(1 + z).*(1 + z).*(1 + z) + omega(2)).^0.5;
dydz = 1./ Efnc;
