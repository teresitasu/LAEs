LAE output files at z = 6

r and k are in comoving units in the output files
k is on a logarithmic grid of 701 points from lgk = -4 to 3
r is on a linear grid of 8192 points from 0 to 2148-0.25

For generating PGammaGamma(k), b_chi,delta = 1 and b_chi,Gamma = -0.2 were used.
For generating Pkl values, b_delta = 0.707, b_Gamma = -0.17 and tau_eff = 6.0 were used.

1. LyASolvedGammaCorrAAL_W14B21_bG3_tG100_tQ1_z6.0.out
Columns: (1) k (h/Mpc) (2) PGammaGamma (total) (3) PGammaGamma_nsn (non-shotnoise term) (4) PGammaGamma_sn (shotnoise term) (5) PS (dm power spectrum at z = 0)

2. LyALAEPklAALComps_W14B21_bG3_tG100_tQ1_z6.0.out (output from LyALAEPklComps.m)
Output is in 3 batches of 701 rows for each batch (one row for each k value), output for each output redshift.
The three batches are for l = 0, 2 and 4. So for each of the output redshifts,
Columns: (1) k (h/Mpc) (2) PkLAE0 (3) PkLAEG0 (4) Pk0_LAE_nsn (5) Pk0_LAE (6) Pk (701 rows)
         (1) k (h/Mpc) (2) PkLAE2 (3) PkLAEG2 (4) Pk2_LAE_nsn (5) Pk2_LAE (6) Pk  (701 rows)
         (1) k (h/Mpc) (2) PkLAE4 (3) PkLAEG4 (4) Pk4_LAE_nsn (5) Pk4_LAE (6) Pk  (701 rows)
The final column (Pk) is the dark matter power spectrum at z = 6

3. LyALAEXilfftAALComps_W14B21_bG3_tG100_tQ1_z6.0.out (output from LyALAEXilComps_fft.m)
Columns: (1) r (Mpc/h) (2) r2xi0LL (3) r2xi0LG (4) r2xi0GG_nsn (5) r2xi0GG (6) r2xi2LL, (7) r2xi2LG  (8) r2xi2GG_nsn  (9) r2xi2GG  (10) r2xi4LL  (11) r2xi4LG (12) r2xi4GG_nsn  (13) r2xi4GG

4. LyALAEXiangAALComps_W14B21_bG3_tG100_tQ1_z6.0.out (output from LyALAEXiangComps.m)
Columns: (1) theta (degrees) (2) xiaLL (3) xiaLG (4) xiaGG_nsn (5) xiaGG (6) xia_nsn (7) xia
Here,
      xia_nsn = xiaLL + xiaLG + xiaGG_nsn
and
      xia = xiaLL + xiaLG + xiaGG

5. LyAMakeTableGammakCoeffs_W14B21_lzMG_2561.out
Columns: (1) redshift (2) phi*(chi + zeta) (3) chi