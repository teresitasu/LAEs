more off;
[om_m,om_v,om_bh2,h,an,sigma8] = cdenCosparamInit;
aj = 1.8;
aS = 0.8;
%abphi = 2.0;
zi = 8;
zf = 4.8;
lzm = 32;
bfrac = 0.;
bet = 1.2;
cfrac = 1.;
iQmod = 3;
%bj = 1;
bnH = 1;
baA = -0.5;
bLLS = 1;
tau_Q = 10;
bG = 3;
b_LAE = 3;
a_LAE = -2.5;
tau_G = 100;
%b_delta = -0.17;
%b_Gamma = 0.13;
b_delta = [0.608,0.663,0.654,0.600,0.567,0.542,0.516,0.707];
b_Gamma = [-0.193,-0.313,-0.309,-0.225,-0.208,-0.200,-0.192,-0.173];
tau_eff = [1.58,1.951,2.216,2.8,3.8,4.25,4.7,6.0];
f_LAE = [0.2568,0.2111,0.1820,0.1283,0.0661,0.0479,0.0344,0.0125]; %z0 = 0, sigh = sig
f2_LAE = [0.1262,0.0957,0.0778,0.0481,0.0199,0.0131,0.0085,0.0024]; %z0 = 0, sigh = sig
%f_LAE = [0.0528,0.0393,0.0316,0.0192,0.0078,0.0051,0.0034,0.0010]; % z0 = 1, sigh = sig
%f2_LAE = [0.0117,0.0074,0.0053,0.0025,0.0006,0.0003,0.0002,0.0000]; % z0 = 1, sigh = sig
%f_LAE = [0.3445,0.3119,0.2899,0.2455,0.1826,0.1594,0.1391,0.0937]; %z0 = 0, sigh = 2*sig
%f2_LAE = [0.2540,0.2268,0.2090,0.1743,0.1277,0.1111,0.0968,0.0654]; %z0 = 0, sigh = 2*sig
%f_LAE = [0.1506,0.1302,0.1170,0.0920,0.0603,0.0497,0.0411,0.0238]; % z0 = 1, sigh = 2*sig
%f2_LAE = [0.0923,0.0783,0.0695,0.0536,0.0346,0.0285,0.0236,0.0139]; % z0 = 1, sigh = 2*sig
% Convert to strings
saj = num2str(aj,'%.1f');
saS = num2str(aS,'%.1f');
%sabphi = num2str(abphi,'%.1f');
szi = num2str(zi,'%.1f');
szf = num2str(zf,'%.1f');
slzm = num2str(lzm,'%i');
sbfrac = num2str(bfrac,'%.1f');
sbet = num2str(bet,'%.1f');
scfrac = num2str(cfrac,'%.1f');
siQmod= num2str(iQmod,'%.1f');
%sbj = num2str(bj,'%.1f');
sbnH = num2str(bnH,'%.1f');
sbaA = num2str(abs(baA),'%.1f');
sbLLS = num2str(bLLS,'%.1f');
stau_Q = num2str(tau_Q,'%.1f');
sbG = num2str(bG,'%.1f');
sbLAE = num2str(b_LAE,'%.1f');
smaLAE = num2str(-a_LAE,'%.1f');
stau_G = num2str(tau_G,'%.1f');
%sb_d = num2str(b_delta,'%.1f');
%sb_G = num2str(b_Gamma,'%.1f');
suffix = '_aj_';
suffix = strcat(suffix,saj);
suffix = strcat(suffix,'_aS_');
suffix = strcat(suffix,saS);
%suffix = strcat(suffix,'_abphi_');
%suffix = strcat(suffix,sabphi);
suffix = strcat(suffix,'_zi_');
suffix = strcat(suffix,szi);
suffix = strcat(suffix,'_zf_');
suffix = strcat(suffix,szf);
suffix = strcat(suffix,'_lzm_');
suffix = strcat(suffix,slzm);
suffix = strcat(suffix,'_bf_');
suffix = strcat(suffix,sbfrac);
suffix = strcat(suffix,'_bet_');
suffix = strcat(suffix,sbet);
suffix = strcat(suffix,'_tf_');
suffix = strcat(suffix,scfrac);
suffix = strcat(suffix,'_iQ_');
suffix = strcat(suffix,siQmod);
%suffix = strcat(suffix,'_bj_');
%suffix = strcat(suffix,sbj);
suffix = strcat(suffix,'_bnH_');
suffix = strcat(suffix,sbnH);
if(baA > 0)
  suffix = strcat(suffix,'_baA_');
 else
   suffix = strcat(suffix,'_baA_m');
end
suffix = strcat(suffix,sbaA);
suffix = strcat(suffix,'_bLLS_');
suffix = strcat(suffix,sbLLS);
suffix = strcat(suffix,'_tQ_');
suffix = strcat(suffix,stau_Q);
suffix = strcat(suffix,'_bG_');
suffix = strcat(suffix,sbG);
suffix = strcat(suffix,'_bLAE_');
suffix = strcat(suffix,sbLAE);
suffix = strcat(suffix,'_maLAE_');
suffix = strcat(suffix,smaLAE);
suffix = strcat(suffix,'_tG_');
suffix = strcat(suffix,stau_G);
%suffix = strcat(suffix,'_bd_');
%suffix = strcat(suffix,sb_d);
%suffix = strcat(suffix,'_bGam_');
%suffix = strcat(suffix,sb_G);
%suffix = strcat(suffix,'_kext_0.10_lgk_m4t3');  %match to k range used
%suffix = strcat(suffix,'_kext_0.30_lgk_m4t3');  %match to k range used
% Use AAE to denote added asymptotic extrapolation for high k is used
% Use AAL to denote used asymptotic lorentzian approximation for all k
prefix = 'LyALAEPklAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z001'; %z0 = 0, sigh = sig
%prefix = 'LyALAEPklAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z011'; %z0 = 1, sigh = sig (was _z01)
%prefix = 'LyALAEPklAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z002'; %z0 = 0, sigh = 2*sig
%prefix = 'LyALAEPklAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z012'; %z0 = 1, sigh = 2*sig
outfile = strcat(prefix,suffix);
finalmatoutfile = strcat(outfile,'.mat');
outfileg = outfile;
%zred_out = [4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5];
zred_out = [4.8,5.0,5.2,5.4,5.6,5.7,5.8,6.0];
lzred_out = length(zred_out);
[fk,Pk,PkLAE0,PkLAEG0,Pk0_LAE_nsn,Pk0_LAE] = LyALAEPklComps(zred_out,b_LAE,b_delta,b_Gamma,a_LAE,f_LAE,f2_LAE,tau_eff,0);
[fk,Pk,PkLAE2,PkLAEG2,Pk2_LAE_nsn,Pk2_LAE] = LyALAEPklComps(zred_out,b_LAE,b_delta,b_Gamma,a_LAE,f_LAE,f2_LAE,tau_eff,2);
[fk,Pk,PkLAE4,PkLAEG4,Pk4_LAE_nsn,Pk4_LAE] = LyALAEPklComps(zred_out,b_LAE,b_delta,b_Gamma,a_LAE,f_LAE,f2_LAE,tau_eff,4);
for iz = 1:lzred_out
  Aout = [fk(1,:)',PkLAE0(iz,:)',PkLAEG0(iz,:)',Pk0_LAE_nsn(iz,:)',Pk0_LAE(iz,:)',Pk(iz,:)'];
  Bout = [fk(1,:)',PkLAE2(iz,:)',PkLAEG2(iz,:)',Pk2_LAE_nsn(iz,:)',Pk2_LAE(iz,:)',Pk(iz,:)'];
  Cout = [fk(1,:)',PkLAE4(iz,:)',PkLAEG4(iz,:)',Pk4_LAE_nsn(iz,:)',Pk4_LAE(iz,:)',Pk(iz,:)'];
  szout = num2str(zred_out(iz),'%.1f');
  outfilez = strcat(outfileg,'_zout_');
  outfilez = strcat(outfilez,szout);
  outfilez = strcat(outfilez,'.out');
  fido = fopen(outfilez,'w');
  fprintf(fido,'%8.6f %16.6e %16.6e %16.6e %16.6e %16.6e\n',Aout');
  fprintf(fido,'%8.6f %16.6e %16.6e %16.6e %16.6e %16.6e\n',Bout');
  fprintf(fido,'%8.6f %16.6e %16.6e %16.6e %16.6e %16.6e\n',Cout');
  fclose(fido);
end
save(finalmatoutfile, 'fk', 'Pk', 'PkLAE0', 'PkLAEG0', 'Pk0_LAE_nsn', 'Pk0_LAE', 'PkLAE2', 'PkLAEG2', 'Pk2_LAE_nsn', 'Pk2_LAE', 'PkLAE4', 'PkLAEG4', 'Pk4_LAE_nsn', 'Pk4_LAE', 'zred_out', 'bG', 'b_LAE', 'b_delta', 'b_Gamma', 'tau_eff');
