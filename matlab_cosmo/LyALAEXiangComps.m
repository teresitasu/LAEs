%
%***************************************************************
%* [theta,xiaLL,xiaLG,xiaGG_nsn,xiaGG] = LyALAEXiangComps(sigr);
%***************************************************************
%***************************************************************
%
% Returns angular correlation functions from
% LyALAEXillComps_fft output, for density, density-Gamma and Gamma
% fluctuations separately.
%
% ARGUMENTS
%  sigr        Comoving width of top-hat shell in r-space (cMpc/h)
%
% RETURNS
%  theta        Angular separation (arcsecs)
%  xiaLL        LyA emitter system-system auto-correlation contribution
%  xiaLG        LyA emitter system-Gamma cross-correlation contribution
%  xiaGG_nsn    Gamma-Gamma auto-correlations without shot noise contribution
%  xialGG       Gamma-Gamma auto-correlations with shot noise contribution
%
% COMPATIBILITY: Octave
%
% REQUIREMENTS:
%	         LyAGetDAng.m, cdenCosparamInit.m
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  04 10 21 Creation date.
%
function [theta,xiaLL,xiaLG,xiaGG_nsn,xiaGG] = LyALAEXiangComps(sigr);
if(exist('LyALAEXilfftComps.mat')==2)
  disp('using existing LyALAEXilfftComps.mat file');
  load('LyALAEXilfftComps.mat');
else
  disp('no LyALAEXilfftComps.mat file');
  return;
end
lenz = length(zred_out);
lenr = length(r);
theta = zeros(lenz,lenr);
xiaLL = zeros(lenz,lenr);
xiaLG = zeros(lenz,lenr);
xiaGG = zeros(lenz,lenr);
xiaGG_nsn = zeros(lenz,lenr);
xia = zeros(lenz,lenr);
xia_nsn = zeros(lenz,lenr);
r2 = r.*r; %r is comoving in cMpc/h
rmax = r(lenr);
[DAng,hubb] = LyAGetDAng(zred_out); %DAng in Mpc (proper)
DAng = hubb*(1 + zred_out).*DAng; %express DAng in cMpc/h
for iz = 1:lenz
  thetar = r./ DAng(iz);
  thetas = thetar/ 4.8481e-6; %theta in arcsec
  maskthp = find(thetas>0);
  maskthmax = find(thetas<10000);
  maskth = intersect(maskthp,maskthmax);
  lenmth = length(maskth);
  theta(iz,1:lenmth) = thetas(maskth);
  % integrate along line of sight using 16-pt gauss-legendre quadrature
  nglwt = 16;
  glxpt = [0.0052995   0.0277125   0.0671844   0.1222978   0.1910619 ...
           0.2709916   0.3591982 0.4524937   0.5475063   0.6408018 ...
           0.7290084   0.8089381   0.8777022   0.9328156 0.9722875   0.9947005];
  glwt = [0.013576   0.031127   0.047579   0.062314   0.074798 ...
          0.084578   0.091302   0.094725 0.094725   0.091302   0.084578   0.074798   0.062314   0.047579   0.031127   0.013576];
  % integrate along line of sight using 8-pt gauss-legendre quadrature
  %  nglwt = 8;
  %  glxpt = [0.019855, 0.101667, 0.237234, 0.408283, 0.591717, 0.762766, ...
  %           0.898333, 0.980145];
  %  glwt = [0.050614, 0.111191, 0.156853, 0.181342, 0.181342, 0.156853, ...
  %          0.111191, 0.050614];
  for igl = 1:nglwt
      ugl = sigr*glxpt(igl);
      uth = (ugl.*ugl + DAng(iz)*thetar(maskth)).^0.5;
      iuth = floor(lenr*uth/ rmax);
      mu = 1./ (1 + (ugl./ (DAng(iz)*thetar(maskth))).^2).^0.5;
      L2 = (3*mu.*mu - 1)/ 2;
      L4 = (5*mu.*mu.*(7*mu.*mu - 6) + 3)/ 8;
      xiaLL(iz,1:lenmth) = xiaLL(iz,1:lenmth) + glwt(igl)*(r2xi0LL(iz,iuth) + L2.*r2xi2LL(iz,iuth) + ...
          L4.*r2xi4LL(iz,iuth))./ r2(iuth);
      %      xiaLL(iz,1:lenmth) = xiaLL(iz,1:lenmth)./ r2(iuth);
      xiaLG(iz,1:lenmth) = xiaLG(iz,1:lenmth) + glwt(igl)*(r2xi0LG(iz,iuth) + L2.*r2xi2LG(iz,iuth) + ...
          L4.*r2xi4LG(iz,iuth))./ r2(iuth);
      %      xiaLG(iz,1:lenmth) = xiaLG(iz,1:lenmth)./ r2(iuth);
      xiaGG(iz,1:lenmth) = xiaGG(iz,1:lenmth) + glwt(igl)*(r2xi0GG(iz,iuth) + L2.*r2xi2GG(iz,iuth) + ...
          L4.*r2xi4GG(iz,iuth))./ r2(iuth);
      %      xiaGG(iz,1:lenmth) = xiaGG(iz,1:lenmth)./ r2(iuth);
      xiaGG_nsn(iz,1:lenmth) = xiaGG_nsn(iz,1:lenmth) + glwt(igl)* ...
          (r2xi0GG_nsn(iz,iuth) + L2.*r2xi2GG_nsn(iz,iuth) + L4.*r2xi4GG_nsn(iz,iuth))./ r2(iuth);
      %      xiaGG_nsn(iz,1:lenmth) = xiaGG_nsn(iz,1:lenmth)./ r2(iuth);
  end
  clear thetar;
  clear thetas;
  clear maskthp;
  clear maskthmax;
  clear maskth;
  clear uth;
  clear mu;
  clear L2;
  clear L4;
end
