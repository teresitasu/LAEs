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
siQmod = num2str(iQmod,'%.1f');
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
%suffix = strcat(suffix,'_lgk_m4t0');  % match to k range used
%suffix = strcat(suffix,'_kext_0.10_lgk_m4t3_rmax_2048_lenr_8192');  % match to k range used
%suffix = strcat(suffix,'_kext_0.30_lgk_m4t3_rmax_2048_lenr_8192');  % match to k range used
% AAE denotes uses added asymptotic expansion for dGammak for high k
% AAL denotes uses asymptotic lorentzian approximation for dGammak for all k
% Note that currently (070518), the shot noise doesn't work for the extention in k
%prefix = 'LyALAEXilfftCompsAAE_W14B21mfp_hiresk2';
prefix = 'LyALAEXilfftAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z001'; %z0 = 0, sigh = sig
%prefix = 'LyALAEXilfftAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z011'; %z0 = 1, sigh = sig (was _z01)
%prefix = 'LyALAEXilfftAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z002'; %z0 = 0, sigh = 2*sig
%prefix = 'LyALAEXilfftAALComps_W14B21mfp_QSOLFKWH19_hiresk2_z012'; %z0 = 1, sigh = 2*sig
outfile = strcat(prefix,suffix);
finalmatoutfile = strcat(outfile,'.mat');
outfileg = outfile;
outfile = strcat(outfile,'.out');
%zred_out = [4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5];
zred_out = [4.8,5.0,5.2,5.4,5.6,5.7,5.8,6.0];
lzred_out = length(zred_out);
[r,r2xi0LL,r2xi0LG,r2xi0GG_nsn,r2xi0GG] = LyALAEXilComps_fft(zred_out,0);
[r,r2xi2LL,r2xi2LG,r2xi2GG_nsn,r2xi2GG] = LyALAEXilComps_fft(zred_out,2);
[r,r2xi4LL,r2xi4LG,r2xi4GG_nsn,r2xi4GG] = LyALAEXilComps_fft(zred_out,4);
[lenz,lenk] = size(r2xi0LL);
for iz = 1:lzred_out
  Aout = [r(1,:)',r2xi0LL(iz,:)',r2xi0LG(iz,:)',r2xi0GG_nsn(iz,:)',r2xi0GG(iz,:)',r2xi2LL(iz,:)',r2xi2LG(iz,:)',r2xi2GG_nsn(iz,:)',r2xi2GG(iz,:)',r2xi4LL(iz,:)',r2xi4LG(iz,:)',r2xi4GG_nsn(iz,:)',r2xi4GG(iz,:)'];
  szout = num2str(zred_out(iz),'%.1f');
  outfilez = strcat(outfileg,'_zout_');
  outfilez = strcat(outfilez,szout);
  outfilez = strcat(outfilez,'.out');
  fido = fopen(outfilez,'w');
  fprintf(fido,'%8.6f %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n',Aout');
  fclose(fido);
end
save(finalmatoutfile, 'zred_out', 'r', 'r2xi0LL','r2xi0LG','r2xi0GG_nsn', 'r2xi0GG', 'r2xi2LL','r2xi2LG','r2xi2GG_nsn', 'r2xi2GG', 'r2xi4LL','r2xi4LG','r2xi4GG_nsn', 'r2xi4GG');

