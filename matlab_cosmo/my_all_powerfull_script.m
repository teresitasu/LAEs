%% Runs everything from scratch
%cdenCosparamInit;
%
%my_zred_out = 6.0;
%my_b_delta = 0.707;
%my_b_Gamma = -0.173;
%my_lorder = 0;
%my_b_LAE = 3.0;
%my_tau_eff = 6.0;
%
%%[fk,Pk,PkH,PkHG,Pk0_alpha_nsn,Pk0_alpha] = LyATransmissionPkl(my_zred_out,my_b_delta,my_b_Gamma,0);
%%[fk,Pk,PkH,PkHG,Pk2_alpha_nsn,Pk2_alpha] = LyATransmissionPkl(my_zred_out,my_b_delta,my_b_Gamma,2);
%%[fk,Pk,PkH,PkHG,Pk4_alpha_nsn,Pk4_alpha] = LyATransmissionPkl(my_zred_out,my_b_delta,my_b_Gamma,4);
%%save('LyATransmissionPkl','Pk0_alpha_nsn','Pk0_alpha','Pk2_alpha_nsn','Pk2_alpha','Pk4_alpha_nsn','Pk4_alpha');
%
%
%%[fk,Pk,Pk0LAE,Pk0LAEG,Pk0_LAE_nsn,Pk0_LAE] =     LyALAEPkl(my_zred_out,my_b_LAE,my_b_delta,my_b_Gamma,my_tau_eff,0);
%%[fk,Pk,Pk2LAE,Pk2LAEG,Pk2_LAE_nsn,Pk2_LAE] = LyALAEPkl(my_zred_out,my_b_LAE,my_b_delta,my_b_Gamma,my_tau_eff,2);
%%[fk,Pk,Pk4LAE,Pk4LAEG,Pk4_LAE_nsn,Pk4_LAE] = LyALAEPkl(my_zred_out,my_b_LAE,my_b_delta,my_b_Gamma,my_tau_eff,4);
%%save('LyALAEPkl.mat','fk','Pk','Pk0_LAE_nsn','Pk0_LAE','Pk2_LAE_nsn','Pk2_LAE','Pk4_LAE_nsn','Pk4_LAE');
%[fk,Pk,Pk0LAE,Pk0LAEG,Pk0_LAE_nsn,Pk0_LAE] = LyALAEPklComps_test(my_zred_out,my_b_LAE,my_b_delta,my_b_Gamma,my_tau_eff,my_lorder);
%
%save('LyALAEPkl.mat','fk','Pk','Pk0_LAE_nsn','Pk0_LAE');
%
%fprintf('Pk_LAE %d \n',size(Pk0_LAE)');
%size(Pk0_LAE)
%data_t = load('LyASolvedGammakCorrSS_z6.out');
%fk = data_t(:,1); % column 1 of the data text file is assigned the variable x
%
%%Pk_LAE,Pk_LAE_nsn, Pk
%hold on 
%plot(fk,Pk0_LAE,'g');
%loglog(fk,Pk0_LAE)
%%plot(fk,Pk0_LAE_nsn,'b');
%%plot(fk,Pk0LAEG,'r');
%%plot(fk,Pk0LAE,'m');
%hold off
%%
%zred_out = my_zred_out;
%[rad,r2xil_nsn,r2xi0LL] = LyALAEXil_fft(my_zred_out,my_lorder);
%%[r,xil_nsn,xil] = LyATransmissionXil_fft(my_zred_out,my_lorder);
%[r,r2xilLL,r2xilLG,r2xilGG_nsn,r2xilGG] = LyALAEXilComps_fft(my_zred_out,my_lorder);
%save('LyALAEXilfftComps.mat', 'rad', 'r2xil_nsn', 'r2xi0LL','zred_out');
%%[theta,xiaLL,xiaLG,xiaGG_nsn,xiaGG] = LyALAEXiangComps(6);
%
%fprintf('TEST %d \n');
%%fprintf('Pk_LAE %f \n',Pk0_LAE(0)');
%%fprintf('Pk_LAE_NSN %f \n',Pk0_LAE_nsn(0)');
%%fprintf('Pk %d \n',Pk0');
%
%%hold on 
%%plot(rad,r2xil_nsn,'b');
%%plot(rad,r2xil,'r');
%%hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs everything from scratch
cdenCosparamInit;

my_zred_out = 6.0;
my_b_delta = 0.707;
my_b_Gamma = -0.17;
my_lorder = 0;
my_b_LAE = 3;
my_tau_eff = 6.0;

%[fk,Pk,PkLAEl,PkLAEGl,Pk_LAE_nsn,Pk_LAE]
[fk,Pk,PkLAE0,PkLAEG0,Pk0_LAE_nsn,Pk0_LAE] = LyALAEPklComps(my_zred_out,my_b_LAE,my_b_delta,my_b_Gamma,my_tau_eff,0);
[fk,Pk,PkLAE2,PkLAEG2,Pk2_LAE_nsn,Pk2_LAE] = LyALAEPklComps(my_zred_out,my_b_LAE,my_b_delta,my_b_Gamma,my_tau_eff,2);
[fk,Pk,PkLAE4,PkLAEG4,Pk4_LAE_nsn,Pk4_LAE] = LyALAEPklComps(my_zred_out,my_b_LAE,my_b_delta,my_b_Gamma,my_tau_eff,4);
save('LyALAEPklComps.mat', 'fk','Pk','PkLAE0','PkLAEG0','Pk0_LAE_nsn','Pk0_LAE','PkLAE2','PkLAEG2','Pk2_LAE_nsn','Pk2_LAE','PkLAE4','PkLAEG4','Pk4_LAE_nsn','Pk4_LAE');

fk(1),PkLAE0(1),PkLAEG0(1),Pk0_LAE_nsn(1),Pk0_LAE(1),Pk(1)
fk(701),PkLAE0(701),PkLAEG0(701),Pk0_LAE_nsn(701),Pk0_LAE(701),Pk(701)

fk(1),PkLAE2(1),PkLAEG2(1),Pk2_LAE_nsn(1),Pk2_LAE(1),Pk(1)
fk(701),PkLAE2(701),PkLAEG2(701),Pk2_LAE_nsn(701),Pk2_LAE(701),Pk(701)

fk(1),PkLAE4(1),PkLAEG4(1),Pk4_LAE_nsn(1),Pk4_LAE(1),Pk(1)
fk(701),PkLAE4(701),PkLAEG4(701),Pk4_LAE_nsn(701),Pk4_LAE(701),Pk(701)

[r,r2xi0LL,r2xi0LG,r2xi0GG_nsn,r2xi0GG] = LyALAEXilComps_fft(my_zred_out,0);
[r,r2xi2LL,r2xi2LG,r2xi2GG_nsn,r2xi2GG] = LyALAEXilComps_fft(my_zred_out,2);
[r,r2xi4LL,r2xi4LG,r2xi4GG_nsn,r2xi4GG] = LyALAEXilComps_fft(my_zred_out,4);

ss = 4; 
r(ss),r2xi0LL(ss),r2xi0LG(ss),r2xi0GG_nsn(ss),r2xi0GG(ss),r2xi2LL(ss),r2xi2LG(ss),r2xi2GG_nsn(ss),r2xi2GG(ss),r2xi4LL(ss),r2xi4LG(ss),r2xi4GG_nsn(ss),r2xi4GG(ss)

zred_out = my_zred_out;
save('LyALAEXilfftComps.mat','zred_out','r', 'fk','Pk','PkLAE0','PkLAEG0','Pk0_LAE_nsn','Pk0_LAE','r2xi0LL','r2xi0LG','r2xi0GG_nsn','r2xi0GG','r2xi2LL','r2xi2LG','r2xi2GG_nsn','r2xi2GG','r2xi4LL','r2xi4LG','r2xi4GG_nsn','r2xi4GG');

%[theta,xiaLL,xiaLG,xiaGG_nsn,xiaGG] = LyALAEXiangCompsnew(14.2);
[theta,xiaLL,xiaLG,xiaGG_nsn,xiaGG] = LyALAEXiangComps_25012022(14.2);

xia_nsn = xiaLL + xiaLG + xiaGG_nsn;
xia = xiaLL + xiaLG + xiaGG;

ss = 1;
theta(ss),xiaLL(ss),xiaLG(ss),xiaGG_nsn(ss),xiaGG(ss),xia_nsn(ss),xia(ss);

data_xi = load('LyALAEXiangAALComps_W14B21_bG3_tG100_tQ1_z6.0.out');
av_theta = transpose(data_xi(:,1)); % co lumn 1 of the data text file is assigned the variable x
av_xia_nsn = transpose(data_xi(:,6)); % column 2 is assigned the variable y
av_xia = transpose(data_xi(:,7));

loglog(theta,xia_nsn,'b');
hold all 
loglog(theta,xia,'r');
ylim([1d-3, 10]);
xlim([10, 1d3]);

loglog(av_theta,av_xia_nsn,'b--');
loglog(av_theta,av_xia,'r--');
ylim([1d-3, 10]);
xlim([10, 1d3]);
hold off
