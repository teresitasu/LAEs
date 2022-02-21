%
%********************************************************
%* [r,r2xilLL,r2xilLG,r2xilGG_nsn,r2xilGG] = LyALAEXilComps_fft(zred_out,lorder);
%********************************************************
%********************************************************
%
% Returns 3D spatial correlation function for given legendre order from
% LyALAEPklComps output based on density, density-Gamma and Gamma
% fluctuations separately..
%
% ARGUMENTS
% zred_out  Output redshifts used for dGammakCorr arrays (low to high)
% lorder    Legendre order
%
% RETURNS
%  rad            Separation (Mpc/ h)
%  r2xilLL        LyA emitter system-system auto-correlation contribution
%  r2xilLG        LyA emitter system-Gamma cross-correlation contribution
%  r2xilGG_nsn    Gamma-Gamma auto-correlations without shot noise contribution
%  r2xilGG        Gamma-Gamma auto-correlations with shot noise contribution
%
% COMPATIBILITY: Octave
%
% REQUIREMENTS:
%	         cdenCosparamInit.m called previously
%
%
% AUTHOR: Avery Meiksin
%
% HISTORY:
%  01 10 21 Creation date. (Adapted from LyATransmissionXilComps_fft.m.)
%
function [rad,r2xilLL,r2xilLG,r2xilGG_nsn,r2xilGG] = LyALAEXilComps_fft(zred_out,lorder);
lok = 0;
if(lorder==0)
    lok = 1;
end
if(lorder==2)
    lok = 1;
end
if(lorder==4)
    lok = 1;
end
if(lok==0)
    disp('need lorder = 0, 2 or 4');
    return;
end
lenz = length(zred_out);
if(exist('LyALAEPklComps.mat')==2)
  disp('using existing LyALAEPklComps.mat file');
  load('LyALAEPklComps.mat');
else
  disp('no LyALAEPklComps.mat file');
  return;
end
lenk = length(fk);
%rmax = 1024.; %comoving length in Mpc/ h
rmax = 2048.; %comoving length in Mpc/ h; well-converged to 300 cMpc/h (010618)
%rmax = 4096.; %comoving length in Mpc/ h
%lenr = 1024;
%lenr = 4096;
lenr = 8192; %well-converged to 300 cMpc/h (010618)
%lenr = 16384;
%lenr = 32768;
drad = rmax/ lenr;
rad = linspace(0,rmax-drad,lenr);
maskrp = find(rad>0);
dfk = 2*pi/ rmax;
fkk = linspace(0,(lenr-1)*dfk,lenr);
Yk = zeros(5,lenr);
xilLL = zeros(lenz,lenr);
xilLG = zeros(lenz,lenr);
xilGG = zeros(lenz,lenr);
Yk_nsn = zeros(5,lenr);
xilGG_nsn = zeros(lenz,lenr);
for iz = 1:lenz
    if(lorder==0)
        Pk_LL_z(1,:) = PkLAE0(iz,:);
        Pk_LG_z(1,:) = PkLAEG0(iz,:);
        Pk_GG_z(1,:) = Pk0_LAE(iz,:) - PkLAE0(iz,:) - PkLAEG0(iz,:);
        Pk_GG_nsn_z(1,:) = Pk0_LAE_nsn(iz,:) - PkLAE0(iz,:) - PkLAEG0(iz,:);
    end
    if(lorder==2)
        Pk_LL_z(1,:) = PkLAE2(iz,:);
        Pk_LG_z(1,:) = PkLAEG2(iz,:);
        Pk_GG_z(1,:) = Pk2_LAE(iz,:) - PkLAE2(iz,:) - PkLAEG2(iz,:);
        Pk_GG_nsn_z(1,:) = Pk2_LAE_nsn(iz,:) - PkLAE2(iz,:) - PkLAEG2(iz,:);
    end
    if(lorder==4)
        Pk_LL_z(1,:) = PkLAE4(iz,:);
        Pk_LG_z(1,:) = PkLAEG4(iz,:);
        Pk_GG_z(1,:) = Pk4_LAE(iz,:) - PkLAE4(iz,:) - PkLAEG4(iz,:);
        Pk_GG_nsn_z(1,:) = Pk4_LAE_nsn(iz,:) - PkLAE4(iz,:) - PkLAEG4(iz,:);
    end
    YkLL = LyATransXilInt_fft(fk,fkk,Pk_LL_z,lorder);
    YkLG = LyATransXilInt_fft(fk,fkk,Pk_LG_z,lorder);
    YkGG = LyATransXilInt_fft(fk,fkk,Pk_GG_z,lorder);
    YkGG_nsn = LyATransXilInt_fft(fk,fkk,Pk_GG_nsn_z,lorder);
    if(lorder==0)
        Yk1(1,:) = YkLL(1,:);
        Yr = ifft(Yk1)/ drad;
        Yr(maskrp) = Yr(maskrp)./ rad(maskrp);
        xilLL(iz,:) = imag(Yr(:));
        clear Yk1;
        Yk1(1,:) = YkLG(1,:);
        Yr = ifft(Yk1)/ drad;
        Yr(maskrp) = Yr(maskrp)./ rad(maskrp);
        xilLG(iz,:) = imag(Yr(:));
        clear Yk1;
        Yk1(1,:) = YkGG(1,:);
        Yr = ifft(Yk1)/ drad;
        Yr(maskrp) = Yr(maskrp)./ rad(maskrp);
        xilGG(iz,:) = imag(Yr(:));
        clear Yk1;
        Yk1(1,:) = YkGG_nsn(1,:);
        Yr = ifft(Yk1)/ drad;
        Yr(maskrp) = Yr(maskrp)./ rad(maskrp);
        xilGG_nsn(iz,:) = imag(Yr(:));
    end
    if(lorder==2)
        Yk1(1,:) = YkLL(1,:);
        Yk2(1,:) = YkLL(2,:);
        Yk3(1,:) = YkLL(3,:);
        Yr1 = -ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ (rad(maskrp).*rad(maskrp).*rad(maskrp));
        Yr2 = -ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp);
        Yr3 = -ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ (rad(maskrp).*rad(maskrp));
        xilLL(iz,:) = imag(Yr1(:) + Yr2(:)) + real(Yr3(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        Yk1(1,:) = YkLG(1,:);
        Yk2(1,:) = YkLG(2,:);
        Yk3(1,:) = YkLG(3,:);
        Yr1 = -ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ (rad(maskrp).*rad(maskrp).*rad(maskrp));
        Yr2 = -ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp);
        Yr3 = -ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ (rad(maskrp).*rad(maskrp));
        xilLG(iz,:) = imag(Yr1(:) + Yr2(:)) + real(Yr3(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        Yk1(1,:) = YkGG(1,:);
        Yk2(1,:) = YkGG(2,:);
        Yk3(1,:) = YkGG(3,:);
        Yr1 = -ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ (rad(maskrp).*rad(maskrp).*rad(maskrp));
        Yr2 = -ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp);
        Yr3 = -ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ (rad(maskrp).*rad(maskrp));
        xilGG(iz,:) = imag(Yr1(:) + Yr2(:)) + real(Yr3(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        Yk1(1,:) = YkGG_nsn(1,:);
        Yk2(1,:) = YkGG_nsn(2,:);
        Yk3(1,:) = YkGG_nsn(3,:);
        Yr1 = -ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ (rad(maskrp).*rad(maskrp).*rad(maskrp));
        Yr2 = -ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp);
        Yr3 = -ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ (rad(maskrp).*rad(maskrp));
        xilGG_nsn(iz,:) = imag(Yr1(:) + Yr2(:)) + real(Yr3(:));
    end
    if(lorder==4)
        Yk1(1,:) = YkLL(1,:);
        Yk2(1,:) = YkLL(2,:);
        Yk3(1,:) = YkLL(3,:);
        Yk4(1,:) = YkLL(4,:);
        Yk5(1,:) = YkLL(5,:);
        Yr1 = ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ rad(maskrp).^5;
        Yr2 = ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp).^3;
        Yr3 = ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ rad(maskrp);
        Yr4 = ifft(Yk4)/ drad;
        Yr4(maskrp) = Yr4(maskrp)./ rad(maskrp).^4;
        Yr5 = ifft(Yk5)/ drad;
        Yr5(maskrp) = Yr5(maskrp)./ rad(maskrp).^2;
        xilLL(iz,:) = imag(Yr1(:) + Yr2(:) + Yr3(:)) + real(Yr4(:) + Yr5(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        clear Yk4;
        clear Yk5;
        Yk1(1,:) = YkLG(1,:);
        Yk2(1,:) = YkLG(2,:);
        Yk3(1,:) = YkLG(3,:);
        Yk4(1,:) = YkLG(4,:);
        Yk5(1,:) = YkLG(5,:);
        Yr1 = ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ rad(maskrp).^5;
        Yr2 = ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp).^3;
        Yr3 = ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ rad(maskrp);
        Yr4 = ifft(Yk4)/ drad;
        Yr4(maskrp) = Yr4(maskrp)./ rad(maskrp).^4;
        Yr5 = ifft(Yk5)/ drad;
        Yr5(maskrp) = Yr5(maskrp)./ rad(maskrp).^2;
        xilLG(iz,:) = imag(Yr1(:) + Yr2(:) + Yr3(:)) + real(Yr4(:) + Yr5(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        clear Yk4;
        clear Yk5;
        Yk1(1,:) = YkGG(1,:);
        Yk2(1,:) = YkGG(2,:);
        Yk3(1,:) = YkGG(3,:);
        Yk4(1,:) = YkGG(4,:);
        Yk5(1,:) = YkGG(5,:);
        Yr1 = ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ rad(maskrp).^5;
        Yr2 = ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp).^3;
        Yr3 = ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ rad(maskrp);
        Yr4 = ifft(Yk4)/ drad;
        Yr4(maskrp) = Yr4(maskrp)./ rad(maskrp).^4;
        Yr5 = ifft(Yk5)/ drad;
        Yr5(maskrp) = Yr5(maskrp)./ rad(maskrp).^2;
        xilGG(iz,:) = imag(Yr1(:) + Yr2(:) + Yr3(:)) + real(Yr4(:) + Yr5(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        clear Yk4;
        clear Yk5;
        Yk1(1,:) = YkGG_nsn(1,:);
        Yk2(1,:) = YkGG_nsn(2,:);
        Yk3(1,:) = YkGG_nsn(3,:);
        Yk4(1,:) = YkGG_nsn(4,:);
        Yk5(1,:) = YkGG_nsn(5,:);
        Yr1 = ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ rad(maskrp).^5;
        Yr2 = ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp).^3;
        Yr3 = ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ rad(maskrp);
        Yr4 = ifft(Yk4)/ drad;
        Yr4(maskrp) = Yr4(maskrp)./ rad(maskrp).^4;
        Yr5 = ifft(Yk5)/ drad;
        Yr5(maskrp) = Yr5(maskrp)./ rad(maskrp).^2;
        xilGG_nsn(iz,:) = imag(Yr1(:) + Yr2(:) + Yr3(:)) + real(Yr4(:) + Yr5(:));
    end
end
r2=repmat(rad.*rad,lenz,1);
r2xilLL = r2.*xilLL;
r2xilLG = r2.*xilLG;
r2xilGG = r2.*xilGG;
r2xilGG_nsn = r2.*xilGG_nsn;
