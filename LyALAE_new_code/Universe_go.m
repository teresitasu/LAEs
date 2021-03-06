LyASolvedGammakCorrAsympApprox_Run;
%LyALAEPklComps_Run;




load('EmissMG.mat')

Aout = [Gamma_MG(1,:)',xHI(1,:)',aeffd_MG(1,:)',aeffLLS_MG(1,:)',otsfac_MG(1,:)',zred_MG(1,:)'];


fido = fopen('EmissMG.txt','w');
fprintf(fido,'%16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n',Aout');
fclose(fido);

[e24,Gamma_MG,xHI,aeffd_MG,aeffLLS_MG,otsfac_MG,err] = LyASolveEmissMG(Gamma_i,aj,aS,om_m,om_v,h,bfrac,bet,cfrac,zred_MG);
% save('EmissMG.mat','e24','Gamma_MG','xHI','aeffd_MG','aeffLLS_MG','otsfac_MG','zred_MG','err');
