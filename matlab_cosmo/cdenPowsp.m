%
%********************
%*  PS = cdenPowsp(fk,pnorm,Gamma,an,Tk,ips)  *
%********************
%********************
% Computes matter power spectrum.
% Assumes k in Tk file is in conventional comoving h/ Mpc.
% Conventionally pnorm is for P(k) in units (comoving Mpc/h)^3,
% so PS is returned in units (comoving Mpc/h)^3.
%
% ARGUMENTS
%  fkp        Array of wavenumbers.
%  pnorm      Normalization of power spectrum.
%  Gamma      Curvature of power spectrum.
%  an         Tilt of power spectrum.
%  Tk         Transfer function.
%  ips        Method for computing power spectrum.
%
% COMPATIBILITY: Matlab(?), Octave
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  19 04 16 Creation date. (After lss/src/powsp.f.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function PS = cdenPowsp(fk,pnorm,Gamma,an,Tk,ips)
if(ips == 1)
% use B&E
  a = 6.4/ Gamma;
  b = 3.0/ Gamma;
  c = 1.7/ Gamma;
  PS = 1. + (fk.*(a + c*c*fk) + (b*fk).^1.5).^1.13;
  PS = pnorm*fk.^an./ PS.^(2./ 1.13);
end
if(ips == 2)
% use BBKS
  q = fk/ Gamma;
  PS = 1 + q.*(3.89+q.*(16.1*16.1+q.*(5.46*5.46*5.46+q*6.71*6.71*6.71*6.71)));
  PS = log(1. + 2.34.*q)./ PS.^0.25./ (2.34*q);
  PS = pnorm.*fk.^an.*PS.*PS;
end
if(ips == 3)
  tki = interp1(Tk.fk,Tk.tk,fk,'spline','extrap');
  PS = pnorm*(fk.^an).*tki.*tki;
  maskl = find(fk < Tk.fkmin);
  masku = find(fk > Tk.fkmax);
  q = fk(maskl)/ Gamma;
  PS(maskl) = 1 + q.*(3.89+q.*(16.1*16.1+q.*(5.46*5.46*5.46+q*6.71*6.71*6.71*6.71)));
  PS(maskl) = log(1. + 2.34.*q)./ PS(maskl).^0.25./ (2.34*q);
  PS(maskl) = pnorm.*fk(maskl).^an.*PS(maskl).*PS(maskl);
  clear q;
  q = fk(masku)/ Gamma;
  PS(masku) = 1 + q.*(3.89+q.*(16.1*16.1+q.*(5.46*5.46*5.46+q*6.71*6.71*6.71*6.71)));
  PS(masku) = log(1. + 2.34.*q)./ PS(masku).^0.25./ (2.34*q);
  PS(masku) = pnorm.*fk(masku).^an.*PS(masku).*PS(masku);
end

