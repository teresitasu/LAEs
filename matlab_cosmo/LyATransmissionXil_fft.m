% STILL TO WRITE
%
%********************************************************
%* [r,xil_nsn,xil] = LyATransmissionXil_fft(zred_out,lorder);
%********************************************************
%********************************************************
%
% Returns 3D spatial correlation function for given legendre order from
% LyATransmissionPk output.
%
% ARGUMENTS
% zred_out  Output redshifts used for dGammakCorr arrays (low to high)
% lorder    Legendre order
%
% RETURNS
%  rad       Separation (Mpc/ h)
%  xil_nsn  LyA transmission fluctuation correlation function without shot noise
%  xil        LyA transmission fluctuation correlation function with shot noise
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
%  15 05 19 Creation date. (Adapted from LyATransmissionXil.m.)
%
function [rad,r2xil_nsn,r2xil] = LyATransmissionXil_fft(zred_out,lorder);
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
if(exist('LyATransmissionPkl.mat')==2)
  disp('using existing LyATransmissionPkl.mat file');
  load('LyATransmissionPkl.mat');
  
  data_t = load('LyASolvedGammakCorr.txt');
  fk = data_t(:,1); % column 1 of the data text file is assigned the variable x
  dGammakCorr_nsn = data_t(:,2); % column 2 is assigned the variable y
  dGammakCorr = data_t(:,3);
  
  data_t = load('cdenPowsp_omm_0.3_omv_0.7_ombh2_0.0_h_0.7_an_1.0_s8_0.8_ips_3.out');
  PS = data_t(:,2); % column 1 of the data text file is assigned the variable x

else
  disp('no LyATransmissionPkl.mat file');
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
xil = zeros(lenz,lenr);
Yk_nsn = zeros(5,lenr);
xil_nsn = zeros(lenz,lenr);
for iz = 1:lenz
    if(lorder==0)
        Pk_alpha_z(1,:) = Pk0_alpha(iz,:);
        Pk_alpha_nsn_z(1,:) = Pk0_alpha_nsn(iz,:);
    end
    if(lorder==2)
        Pk_alpha_z(1,:) = Pk2_alpha(iz,:);
        Pk_alpha_nsn_z(1,:) = Pk2_alpha_nsn(iz,:);
    end
    if(lorder==4)
        Pk_alpha_z(1,:) = Pk4_alpha(iz,:);
        Pk_alpha_nsn_z(1,:) = Pk4_alpha_nsn(iz,:);
    end
    Yk = LyATransXilInt_fft(fk,fkk,Pk_alpha_z,lorder);
    Yk_nsn = LyATransXilInt_fft(fk,fkk,Pk_alpha_nsn_z,lorder);
    if(lorder==0)
        Yk1(1,:) = Yk(1,:);
        Yr = ifft(Yk1)/ drad;
        Yr(maskrp) = Yr(maskrp)./ rad(maskrp);
        xil(iz,:) = imag(Yr(:));
        clear Yk1;
        Yk1(1,:) = Yk_nsn(1,:);
        Yr = ifft(Yk1)/ drad;
        Yr(maskrp) = Yr(maskrp)./ rad(maskrp);
        xil_nsn(iz,:) = imag(Yr(:));
    end
    if(lorder==2)
        Yk1(1,:) = Yk(1,:);
        Yk2(1,:) = Yk(2,:);
        Yk3(1,:) = Yk(3,:);
        Yr1 = -ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ (rad(maskrp).*rad(maskrp).*rad(maskrp));
        Yr2 = -ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp);
        Yr3 = -ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ (rad(maskrp).*rad(maskrp));
        xil(iz,:) = imag(Yr1(:) + Yr2(:)) + real(Yr3(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        Yk1(1,:) = Yk_nsn(1,:);
        Yk2(1,:) = Yk_nsn(2,:);
        Yk3(1,:) = Yk_nsn(3,:);
        Yr1 = -ifft(Yk1)/ drad;
        Yr1(maskrp) = Yr1(maskrp)./ (rad(maskrp).*rad(maskrp).*rad(maskrp));
        Yr2 = -ifft(Yk2)/ drad;
        Yr2(maskrp) = Yr2(maskrp)./ rad(maskrp);
        Yr3 = -ifft(Yk3)/ drad;
        Yr3(maskrp) = Yr3(maskrp)./ (rad(maskrp).*rad(maskrp));
        xil_nsn(iz,:) = imag(Yr1(:) + Yr2(:)) + real(Yr3(:));
    end
    if(lorder==4)
        Yk1(1,:) = Yk(1,:);
        Yk2(1,:) = Yk(2,:);
        Yk3(1,:) = Yk(3,:);
        Yk4(1,:) = Yk(4,:);
        Yk5(1,:) = Yk(5,:);
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
        xil(iz,:) = imag(Yr1(:) + Yr2(:) + Yr3(:)) + real(Yr4(:) + Yr5(:));
        clear Yk1;
        clear Yk2;
        clear Yk3;
        clear Yk4;
        clear Yk5;
        Yk1(1,:) = Yk_nsn(1,:);
        Yk2(1,:) = Yk_nsn(2,:);
        Yk3(1,:) = Yk_nsn(3,:);
        Yk4(1,:) = Yk_nsn(4,:);
        Yk5(1,:) = Yk_nsn(5,:);
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
        xil_nsn(iz,:) = imag(Yr1(:) + Yr2(:) + Yr3(:)) + real(Yr4(:) + Yr5(:));
    end
end
r2=repmat(rad.*rad,lenz,1);
r2xil = r2.*xil;
r2xil_nsn = r2.*xil_nsn;
