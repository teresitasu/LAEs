%
%********************
%*  [pnorm,Gamma,Tk] = cdenPowspInitFile(om_m,om_v,om_bh2,h,an,sigma8,ips,Tk_file,nrows)  *
%********************
%********************
% Initialise matter power spectrum computation using input file list of Tk files.
% Assumes k in Tk file is in conventional comoving h/ Mpc.
% Returned pnorm is for P(k) in units (comoving Mpc/h)^3.
%
% ARGUMENTS
%  om_m       Total mass parameter.
%  om_v       Vacuum energy parameter.
%  om_bh2     Baryonic mass parameter x h^2.
%  h          Hubble constant H0/ (100 km/s/Mpc).
%  an         Tilt of matter power spectrum.
%  sigma8     Sigma on 8 Mpc/h scale to normalize power spectrum to.
%  ips        Method for computing power spectrum.
%             ips = 1: B&E; ips = 2: BBKS; ips = 3: Tk file
%  Tk_file    T(k) file for ips = 3
%  nrows      Number of rows in T(k) file
%
% OUTPUT
%  pnorm      Normalization of power spectrum  ([cMpc/h]^3).
%  Gamma      Curvature of power spectrum.
%  Tk         Transfer function.
%
% COMPATIBILITY: Matlab(?), Octave
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  28 11 17 Creation date. (From cdenPowspInit.m.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [pnorm,Gamma,Tk] = cdenPowspInitFile(om_m,om_v,om_bh2,h,an,sigma8,ips,Tk_file,nrows)
%ng = 8388608;
ng = 1e6;
boxsize = 1000;
tol = 1e-6;
fkl = 2*pi/ boxsize;
fku = 0.5*ng*fkl;
lgfkl = log(fkl);
lgfku = log(fku);
%lgfkl = fkl %TEST for linear integration
%lgfku = fku %TEST for linear integration
%Uses Bunn & White (1997) COBE normalization.
tilt = an - 1;
om_b = om_bh2/ (h*h);
if(om_v > 0)
  dh=1.94e-5*(om_m^(-0.785-0.05*log(om_m)))*exp(-0.95*tilt-0.169*tilt*tilt);
else
  dh=1.95e-5*(om_m^(-0.35-0.19*log(om_m)-0.17*tilt))*exp(-tilt-0.14*tilt*tilt);
end
pnorm = 2*pi*pi*dh*dh*(2997.9)^(3+an);
Gamma = om_m*h/ (exp(om_b*(1 + 1./ om_m)));
if(ips == 3)
    %  Tk_file = input('Enter name of Tk file: ','s');
    %  nr = input('Enter number of rows: ');
  fid = fopen(Tk_file,'r');
  Tk_table = fscanf(fid,'%g %g',[2,nrows]);
  Tk_table = Tk_table';  
  Tk.fk = Tk_table(:,1);
  Tk.fkmin = min(Tk.fk);
  Tk.fkmax = max(Tk.fk);
  Tk.tk = Tk_table(:,2);
else
  Tk = 0;
end
disp('Bunn & White COBE pnorm: ');
disp(pnorm);
% Now normalize to desired sigma8.
xf = 8;
ismooth = 1; %spherical top hat
sigma8_BW = quad(@(u)cdenFNorm(u,pnorm,Gamma,an,Tk,xf,ismooth,ips),lgfkl,lgfku,tol);
%% TEST for linear integration
%fk=linspace(lgfkl,lgfku,1000000);
%fk(1)
%fk(1000000)
%dk = (lgfku-lgfkl)/1000000;
%u = cdenFNorm(fk,pnorm,Gamma,an,Tk,xf,ismooth,ips);
%sigma8_BW = dk*trapz(u);
%%
sigma8_BW = sigma8_BW^0.5;
pnorm = pnorm*(sigma8/ sigma8_BW)*(sigma8/ sigma8_BW);
