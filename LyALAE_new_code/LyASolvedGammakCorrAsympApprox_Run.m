more off;
[om_m,om_v,om_bh2,h,an,sigma8] = cdenCosparamInit;
ips = 3;
Tk_file = 'tkO3111L6889h6766N2242_eh_PLANCK18.dat';
nrows = 650;
Gammai=0.249e-12;
%Gammai=0.852e-12;
aj = 1.8;
aS = 0.8;
%ab = 2; %example
ab = 1; % turn off delta_DM evolution for high z (high chiLLS)
%ab = -1; %let LyASolvedGammakCorrAsymptoticApprox.m compute
zi = 8;
zf = 4.8;
bfrac = 0.;
bet = 1.2;
cfrac = 1.;
iQmod = 3;
bj = 1;
bnH = 1;
baA = -0.5;
bLLS = 1;
tau_Q = 10;
%lgQbolmin = 10;
%lgQbolmax = 15;
M1450min = -31;
M1450max = -21;
bG = 3;
tau_G = 100;
% Convert to strings
saj = num2str(aj,'%.1f');
saS = num2str(aS,'%.1f');
if(ab>0)
 sab = num2str(ab,'%.1f');
end
szi = num2str(zi,'%.1f');
szf = num2str(zf,'%.1f');
sbfrac = num2str(bfrac,'%.1f');
sbet = num2str(bet,'%.1f');
scfrac = num2str(cfrac,'%.1f');
siQmod = num2str(iQmod,'%.1f');
sbj = num2str(bj,'%.1f');
sbnH = num2str(bnH,'%.1f');
sbaA = num2str(abs(baA),'%.1f');
sbLLS = num2str(bLLS,'%.1f');
stau_Q = num2str(tau_Q,'%.1f');
sbG = num2str(bG,'%.1f');
stau_G = num2str(tau_G,'%.1f');
suffix = '_aj_';
suffix = strcat(suffix,saj);
suffix = strcat(suffix,'_aS_');
suffix = strcat(suffix,saS);
if(ab>0)
  suffix = strcat(suffix,'_ab_');
  suffix = strcat(suffix,sab);
end
suffix = strcat(suffix,'_zi_');
suffix = strcat(suffix,szi);
suffix = strcat(suffix,'_zf_');
suffix = strcat(suffix,szf);
suffix = strcat(suffix,'_bf_');
suffix = strcat(suffix,sbfrac);
suffix = strcat(suffix,'_bet_');
suffix = strcat(suffix,sbet);
suffix = strcat(suffix,'_tf_');
suffix = strcat(suffix,scfrac);
suffix = strcat(suffix,'_iQ_');
suffix = strcat(suffix,siQmod);
suffix = strcat(suffix,'_bj_');
suffix = strcat(suffix,sbj);
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
suffix = strcat(suffix,'_tG_');
suffix = strcat(suffix,stau_G);
prefix = 'LyASolvedGammakCorrAAL_W14B21mfp_QSOLFKWH19_hiresk2_lgk_m4t3';
outfile = strcat(prefix,suffix);
finalmatoutfile = strcat(outfile,'.mat');
outfileg = outfile;
outfile = strcat(outfile,'.out');
%zred_out = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0];
zred_out = [4.8,5.0,5.2,5.4,5.6,5.7,5.8,6.0];
lzred_out = length(zred_out);
%fk = [0.0001, 0.001, 0.01, 0.1, 1.];
%fk = [0.0001, 0.001];
%fk = 10.^linspace(-4,0,81); %default
%fk = 10.^linspace(-4,3,351); %hiresk
%fk = 10.^linspace(-4,3,751); %hiresk2
fk = 10.^linspace(-4,3,701); %hiresk2

[pnorm,Gamma,Tk] = cdenPowspInitFile(om_m,om_v,om_bh2,h,an,sigma8, ...
                                     ips,Tk_file,nrows);
PS = cdenPowsp(fk,pnorm,Gamma,an,Tk,ips);
lenfk = length(fk);
dGammakCorr = zeros(lzred_out,lenfk);
dGammakCorr_nsn = zeros(lzred_out,lenfk);
dGammakCorr_sn = zeros(lzred_out,lenfk);
[neff_Q,avgL_Q,avgL2_Q] = LyAGetNeffQSO_KWH19(zred_out,M1450min,M1450max,iQmod);

for ik = 1:lenfk
    k = fk(ik)
    Pk = PS(ik);
%    [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrAsymptoticApprox(Gammai,aj,aS,ab,om_m,h,zi,zf,k,Pk,bfrac,bet,cfrac,iQmod,bj,bnH,baA,bLLS,tau_Q,lgQbolmin,lgQbolmax,bG,tau_G);
    [e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S] = LyASolvedGammakCorrAsymptoticApprox(Gammai,aj,aS,ab,om_m,h,zi,zf,k,Pk,bfrac,bet,cfrac,iQmod,bj,bnH,baA,bLLS,tau_Q,M1450min,M1450max,bG,tau_G);
    disp(dGammaCorr(1,1));
    disp(dGammaCorr_sn(1,1));
%    sk = num2str(k,'%.4f');
%    matoutfile = strcat(prefix,suffix);
%    matoutfile = strcat(matoutfile,'_k_');
%    matoutfile = strcat(matoutfile,sk);
%    matoutfile = strcat(matoutfile,'.mat');
%    if(ik==1)
%      save(matoutfile, 'k', 'Pk', 'e24', 'tau_Q', 'bG', 'tau_G', 'zred', 'dGammaCorr', 'dGammaCorr_nsn', 'dGammaCorr_sn','Gamma', 'aeff_d', 'aeff_LLS', 'S');
%    end
%    if(ik==21)
%      save(matoutfile, 'k', 'Pk', 'e24', 'tau_Q', 'bG', 'tau_G', 'zred', 'dGammaCorr', 'dGammaCorr_nsn', 'dGammaCorr_sn','Gamma', 'aeff_d', 'aeff_LLS', 'S');
%    end
%    if(ik==41)
%      save(matoutfile, 'k', 'Pk', 'e24', 'tau_Q', 'bG', 'tau_G', 'zred', 'dGammaCorr', 'dGammaCorr_nsn', 'dGammaCorr_sn','Gamma', 'aeff_d', 'aeff_LLS', 'S');
%    end
%    if(ik==61)
%      save(matoutfile, 'k', 'Pk', 'e24', 'tau_Q', 'bG', 'tau_G', 'zred', 'dGammaCorr', 'dGammaCorr_nsn', 'dGammaCorr_sn','Gamma', 'aeff_d', 'aeff_LLS', 'S');
%    end
%    if(ik==81)
%      save(matoutfile, 'k', 'Pk', 'e24', 'tau_Q', 'bG', 'tau_G', 'zred', 'dGammaCorr', 'dGammaCorr_nsn', 'dGammaCorr_sn','Gamma', 'aeff_d', 'aeff_LLS', 'S');
%    end
    stop
    for iz = 1:lzred_out
      masku = find(zred>=zred_out(iz));
      izo = masku(1);
      dGammakCorr(iz,ik) = dGammaCorr(izo,izo);
      dGammakCorr_nsn(iz,ik) = dGammaCorr_nsn(izo,izo);
      dGammakCorr_sn(iz,ik) = dGammaCorr_sn(izo,izo);
      Aout = [fk(1,:)',dGammakCorr(iz,:)',dGammakCorr_nsn(iz,:)',dGammakCorr_sn(iz,:)',PS(1,:)'];
%      sizo = num2str(izo,'%i');
%     outfilez = strcat(outfileg,'_izo_');
%     outfilez = strcat(outfilez,sizo);
      szout = num2str(zred_out(iz),'%.1f');
      outfilez = strcat(outfileg,'_zout_');
      outfilez = strcat(outfilez,szout);
      outfilez = strcat(outfilez,'.out');
      fido = fopen(outfilez,'w');
      fprintf(fido,'%8.6f %16.6e %16.6e %16.6e %16.6e\n',Aout');
      fclose(fido);
    end
end
save(finalmatoutfile, 'fk', 'PS', 'e24', 'tau_Q', 'bG', 'tau_G', 'zred', 'dGammakCorr', 'dGammakCorr_nsn', 'dGammakCorr_sn','Gamma', 'aeff_d', 'aeff_LLS', 'S');
