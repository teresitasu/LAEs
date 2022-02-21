import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
import pandas as pd


# Define Cosmology
om_m = 0.3111
om_v = 1.0 - om_m
om_k = 0.0
om_bh2 = 0.02242
h = 0.6766
an = 0.9665
sigma8 = 0.8102

omega = [om_m, om_v, 0., om_bh2]

#Columns: (1) k (h/Mpc) (2) PGammaGamma (total) (3) PGammaGamma_nsn (non-shotnoise term) (4) PGammaGamma_sn (shotnoise term) (5) PS (dm power spectrum at z = 0)
fk,dGammakCorr,dGammakCorr_nsn,PGG,PS = np.loadtxt('./LyASolvedGammakCorrAAL_W14B21mfp_QSOLFKWH19_hiresk2_lgk_m4t3_aj_1.8_aS_0.8_ab_1.0_zi_8.0_zf_4.8_bf_0.0_bet_1.2_tf_1.0_iQ_3.0_bj_1.0_bnH_1.0_baA_m0.5_bLLS_1.0_tQ_10.0_bG_3.0_tG_100.0_zout_6.0.out',unpack = True)

## Load previously calculated data
#data_t = np.loadtxt('./cdenPowsp_omm_0_3_omv_0_7_ombh2_0_0_h_0_7_an_1_0_s8_0_8_ips_3.out')
#PS = data_t[:,1]


def cdenGrowth(aexp):
    '''
     INPUT
      aexp            Cosmological expansion factor aexp = 1/ (1 + z).
    
     OUTPUT
      D               Growth factor, normalised arbitrarily to initial a_i = 0.001.
      dlnDdt          d(lnD)/dt for growth factor D, normalised
                      arbitrarily to initial a_i = 0.001.

     Based on Peebles LSSU Eq.(10.12).
     If aexp is an array, it must be ordered from small to large values and
     with a spacing sufficient for good accuracy using the trapezoidal rule.
    '''
    aexp0 = 0.001
    gint = 0
    lena = np.size(aexp)
    

    if lena > 1:
        gint = np.zeros(lena)
        aexp = np.array(aexp)
        for i in range(lena):
            g1 = integrate.quad(cdenGrowthFunc,aexp0,aexp[i])
            gint[i] = g1[0] + gint[i]
    else:
        g1 = integrate.quad(cdenGrowthFunc,aexp0,aexp)
        gint = g1[0] + gint
    # Debug
    #print(gint,(omega[0]/aexp**3 + omega[1] + omega[2]/aexp**2)**0.5)
    D = gint*(omega[0]/aexp**3 + omega[1] + omega[2]/aexp**2)**0.5

    dlnDdt = -1.5*omega[0]/aexp**3 - omega[2]/aexp**2
    dlnDdt = dlnDdt / (omega[0]/aexp**3 + omega[1] + omega[2]/aexp**2)**0.5
    dlnDdt = dlnDdt + 1.0/ (gint*(omega[0]/aexp + omega[1]*aexp**2 + omega[2]))

    return D , dlnDdt


def cdenGrowthFunc(aexp):
    dfdap = 1.0/ (omega[0]/ aexp + omega[1]*aexp*aexp + omega[2])**1.5
    return dfdap


def LyATransXilInt_fft(fk,fkk,Pk_alpha_z,lorder=0,debug = False):
    from scipy.interpolate import PchipInterpolator as pchip

    #Pk_alpha_z_arr = np.interp(fkk,fk,Pk_alpha_z)
    pk_interp = pchip(fk,Pk_alpha_z, axis=0, extrapolate=True)
    Pk_alpha_z_arr = pk_interp.__call__(fkk)
    n = len(fkk)
    pk_interp = pchip(fk,Pk_alpha_z, axis=0, extrapolate=True)
    maskp = fkk > 0
    Yk = np.zeros([5,n])
    lok = 0
    if lorder == 0:
        Yk[0][maskp] = Pk_alpha_z_arr[maskp]
        lok = 1
    if lorder == 2:
        Yk[0][maskp] = 3*Pk_alpha_z_arr[maskp]/ fkk[maskp]**2
        Yk[1][maskp] = -Pk_alpha_z_arr[maskp]
        Yk[2][maskp] = -3*Pk_alpha_z_arr[maskp]/ fkk[maskp]
        lok = 1
    if lorder == 4:
        Yk[0][maskp] = 105*Pk_alpha_z_arr[maskp]/ fkk[maskp]**4
        Yk[1][maskp] = -45*Pk_alpha_z_arr[maskp]/ fkk[maskp]**2
        Yk[2][maskp] = Pk_alpha_z_arr[maskp]
        Yk[3][maskp] = -105*Pk_alpha_z_arr[maskp]/ fkk[maskp]**3
        Yk[4][maskp] = 10*Pk_alpha_z_arr[maskp]    / fkk[maskp]
        lok = 1

    if lok == 0:
        print('LyATransXilInt: lorder must be 0, 2 or 4');
        return None
    
    for i in range(5):
        Yk[i][maskp] = fkk[maskp] * Yk[i][maskp]/ (2.0 * np.pi)

    for j in np.arange(np.int(n/2)+1,n):
        if lorder == 0:
            
            Yk[0][j] = -Yk[0][n-j]
            if debug:
                if (j == np.int(n/2)+1) or (j == n-1):
                    print("n {} {} {} {:.5e}".format(n,j,n-j,Yk[0][j]))
        if lorder == 2:
            Yk[0][j] = -Yk[0][n-j]
            Yk[1][j] = -Yk[1][n-j]
            Yk[2][j] =  Yk[2][n-j]
        
        if lorder == 4:
            Yk[0][j] = -Yk[0][n-j]
            Yk[1][j] = -Yk[1][n-j]
            Yk[2][j] = -Yk[2][n-j]
            Yk[3][j] =  Yk[3][n-j]
            Yk[4][j] =  Yk[4][n-j]

    return Yk

def LyALAEPkl(zred_out,b_LAE=3.0,b_delta=0.707,
              b_Gamma=-0.173,tau_eff=6,
              a_LAE=1,f_LAE=1,f2_LAE=1, lorder=0,debug=False):
    '''
     Returns 3D power spectrum of LyA emitters, including correlations in
     metagalactic photoionization rate. Allows for QSO bias from Laurent et al. (2017)
     and galaxy bias.

     ARGUMENTS
     zred_out    Output redshifts used for dGammakCorr arrays (low to high)
     b_LAE       LyA emitter bias factor
     b_delta     Gas density bias factor of tau_eff
     b_Gamma     Ionization rate bias factor of tau_eff
     tau_eff     Median tau_eff of IGM
     lorder      Legendre order l

     RETURNS
     fk          Wavenumber (h/ Mpc)
     Pk          Dark matter power spectrum (zred,fk)
     PkLAEl      Legendre component of LAE power spectrum from density alone
     PkLAEGl     Legendre component of LAE flux power spectrum from density and
                            density-Gamma cross term
     Pk_LAE_nsn  LyA emitter power spectrum without shot noise
     Pk_LAE      LyA emitter power spectrum with shot noise

     COMPATIBILITY: Matlab, Octave

     REQUIREMENTS:
                cdenCosparamInit.m called previously


     AUTHOR: Avery Meiksin

     HISTORY:
      01 10 21 Creation date. (Adapted from LyATransmissionPklComps.m.)
    '''
    if lorder == 0:
        lok = 1.0
    elif lorder == 2:
        lok = 1.0
    elif lorder == 4:
        lok = 1.0
    else:
        print("must use lorder = 0, 2 or 4")
        return None

    lenk  = len(fk)
    dGammak_nsn = dGammakCorr_nsn**0.5
    deltk_0 = PS**0.5
    gg_0,gv_0 = cdenGrowth(1.)
    
    
    zp1 = 1.0 + zred_out
    zp1_3 = zp1*zp1*zp1
    Hfac = (om_m*zp1_3 + om_v + om_k*zp1*zp1)**0.5
    aexp_flip = 1./ zp1
    
    gg,gv = cdenGrowth(aexp_flip)

    if debug: print("gg {:.4f}, gv {:.4f}".format(gg,gv))
    if debug: print("Hfac {:.4f}".format(Hfac))
    #gg = fliplr(gg_flip)
    #gv = fliplr(gv_flip)
    deltk = gg * deltk_0/ gg_0
    Pk = deltk * deltk
    ddGammak_nsn = deltk * dGammak_nsn
    
    beta = gv / Hfac
    beta2 = beta * beta
    betafac = np.ones(lenk) * beta
    beta2fac = np.ones(lenk) * beta2
    PkLAEa = np.ones(lenk) * b_LAE * deltk
    PkLAEb = np.ones(lenk) * (a_LAE+1.0)* f_LAE*tau_eff*b_delta* deltk    
    PkLAEb2= np.ones(lenk) * (a_LAE+1.0)*f2_LAE*tau_eff*b_delta* deltk
    PkLAE = PkLAEa**2 + 2.0 * PkLAEa*PkLAEb + PkLAEb**2
    if debug: print("taus {:.4f}".format((b_LAE - tau_eff*b_delta)* \
                    (tau_eff*b_Gamma)))
    PkLAEG = 2.0 *(a_LAE +1.0) *np.ones(lenk)*((b_LAE*f_LAE +\
                    (a_LAE+1.0)*f2_LAE*tau_eff*b_delta)* \
                    (tau_eff*b_Gamma)) * ddGammak_nsn
    if debug: print("PkLAEG {:.4f}".format(PkLAEG[9]))

    if lorder == 0:
        PkLAEl = (1.0 + betafac/ 1.5 + beta2fac/ 5.) * PkLAE
        PkLAEGl = (1.0 + betafac/ 3.) * PkLAEG
        
    if lorder == 2:
        PkLAEl = 4.0 * betafac * (1./ 3 + betafac/ 7.) * PkLAE
        PkLAEGl = (2.0 * betafac/ 3.) * PkLAEG
        b_Gamma = 0.0
        
    if lorder == 4:
        PkLAEl = 8.0 * beta2fac * PkLAE / 35.
        PkLAEGl = 0 * PkLAEl
        b_Gamma = 0.0

    Pk_LAE_nsn = PkLAEl + PkLAEGl + np.ones(lenk)*(a_LAE+1.0)**2 * \
                 f2_LAE*(tau_eff*tau_eff)*(b_Gamma*b_Gamma) * \
                 dGammakCorr_nsn
    Pk_LAE = PkLAEl + PkLAEGl + np.ones(lenk)*(a_LAE+1.0)**2 * \
                 f2_LAE*(tau_eff*tau_eff)*(b_Gamma*b_Gamma) * \
                 dGammakCorr

    return Pk,PkLAEl,PkLAEGl,Pk_LAE,Pk_LAE_nsn


def LyALAEXil_fft(zred_out,lorder=0,b_LAE=3.0,b_delta=0.707,
              b_Gamma=-0.173,tau_eff=6,
              a_LAE=1,f_LAE=1,f2_LAE=1,debug=False):
    '''
    % Returns 3D spatial correlation function for given legendre order from
     LyATransmissionPk output.
    
     INPUT
      zred_out    Output redshifts used for dGammakCorr arrays (low to high)
      lorder      Legendre order
    
     OUTPUT
      rad        Separation (Mpc/ h)
      r2xil_nsn  LyA emitter correlation function (times r^2) without shot noise
      r2xil      LyA emitter correlation function (times r^2) with shot noise
    
    '''
    lenk = len(fk)
    rmax = 2048 #; %comoving length in Mpc/ h; well-converged to 300 cMpc/h (010618)
    lenr = 8192 #; %well-converged to 300 cMpc/h (010618)
    drad = rmax/ lenr
    rad = np.linspace(0,rmax-drad,lenr)
    maskrp = rad > 0
    dfk = 2.0*np.pi/ rmax
    fkk = np.linspace(0,(lenr-1)*dfk,lenr)
    Yk = np.zeros([5,lenr])
    xil = np.zeros(lenr)

    xilLL = np.zeros(lenr)
    xilLG = np.zeros(lenr)
    xilGG = np.zeros(lenr)
    Yk_nsn = np.zeros([5,lenr]) 
    xilGG_nsn = np.zeros(lenr)

    # Pk,PkLAEl,PkLAEGl,Pk_LAE,Pk_LAE_nsn

    Pk_z,PkLAEl_z,PkLAEGl_z,Pk_LAE_z,Pk_LAE_nsn_z = LyALAEPkl(zred_out,b_LAE=b_LAE,
                        b_delta=b_delta,b_Gamma=b_Gamma,tau_eff=tau_eff,
                        a_LAE=a_LAE,f_LAE=f_LAE,f2_LAE=f2_LAE,lorder = lorder)
    
    aa = 9

    PK_LL_z = PkLAEl_z
    PK_LG_z = PkLAEGl_z
    PK_GG_z = Pk_LAE_z - PkLAEl_z - PkLAEGl_z
    PK_GG_nsn_z = Pk_LAE_nsn_z -  PkLAEl_z - PkLAEGl_z
    
    if debug: print('PK_GG_z {:10.3e} {:10.3e} {:10.3e}'.format(Pk_LAE_z[aa], PkLAEl_z[aa], PkLAEGl_z[aa]))
    if debug: print("pkl {:10.3e} {:10.3e} {:10.3e} {:10.3e} ".format(PK_LL_z[aa]
                        ,PK_LG_z[aa],PK_GG_z[aa],PK_GG_nsn_z[aa]))
    

    YkLL = np.array(LyATransXilInt_fft(fk,fkk,PK_LL_z,lorder))
    YkLG = np.array(LyATransXilInt_fft(fk,fkk,PK_LG_z,lorder))
    YkGG = np.array(LyATransXilInt_fft(fk,fkk,PK_GG_z,lorder))
    YkGG_nsn =  np.array(LyATransXilInt_fft(fk,fkk,PK_GG_nsn_z,lorder))
    
    if debug: print("Ykll {:10.3e} {:10.3e} {:10.3e} {:10.3e} ".format(YkLL[0][aa]
                        ,YkLG[0][aa],YkGG[0][aa],YkGG_nsn[0][aa]))
    if lorder == 0:
        Yk1 = YkLL[0]
        Yr = np.fft.ifft(Yk1)/ drad
        Yr[maskrp] = Yr[maskrp]/ rad[maskrp]
        xilLL = np.imag(Yr)
        #---------------
        Yk1 = YkLG[0]
        Yr = np.fft.ifft(Yk1)/ drad
        Yr[maskrp] = Yr[maskrp]/ rad[maskrp]
        xilLG = np.imag(Yr)
        #---------------
        Yk1 = YkGG[0]
        Yr = np.fft.ifft(Yk1)/ drad
        Yr[maskrp] = Yr[maskrp]/ rad[maskrp]
        xilGG = np.imag(Yr)
        #---------------
        Yk1 = YkGG_nsn[0]
        Yr = np.fft.ifft(Yk1)/ drad
        Yr[maskrp] = Yr[maskrp]/ rad[maskrp]
        xilGG_nsn = np.imag(Yr)
        
        if debug: print("xill {:10.3e} {:10.3e} {:10.3e} {:10.3e} ".format(xilLL[aa]
                        ,xilLG[aa],xilGG[aa],xilGG_nsn[aa]))
    if lorder == 2:
        Yk1 = YkLL[0]
        Yk2 = YkLL[1]
        Yk3 = YkLL[2]
        Yr1 = -np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ (rad[maskrp]**3)
        Yr2 = -np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ rad[maskrp]
        Yr3 = -np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ (rad[maskrp]**2)
        xilLL = np.imag(Yr1 + Yr2) + np.real(Yr3)
        # -----------
        #print(Yk1[0],YkLG[0])
        Yk1 = YkLG[0]
        Yk2 = YkLG[1]
        Yk3 = YkLG[2]
        Yr1 = -np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ (rad[maskrp]**3)
        Yr2 = -np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ (rad[maskrp])
        Yr3 = -np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ (rad[maskrp]**2)
        xilLG = np.imag(Yr1 + Yr2) + np.real(Yr3)
        # -----------
        Yk1 = YkGG[0]
        Yk2 = YkGG[1]
        Yk3 = YkGG[2]
        Yr1 = -np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ (rad[maskrp]**3)
        Yr2 = -np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ (rad[maskrp])
        Yr3 = -np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ (rad[maskrp]**2)
        xilGG = np.imag(Yr1 + Yr2) + np.real(Yr3)
        # -----------
        Yk1 = YkGG_nsn[0]
        Yk2 = YkGG_nsn[1]
        Yk3 = YkGG_nsn[2]
        Yr1 = -np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ (rad[maskrp]**3)
        Yr2 = -np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ (rad[maskrp])
        Yr3 = -np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ (rad[maskrp]**2)
        xilGG_nsn = np.imag(Yr1 + Yr2) + np.real(Yr3)

        if debug: print("xill {:10.3e} {:10.3e} {:10.3e} {:10.3e} ".format(xilLL[aa]
                        ,xilLG[aa],xilGG[aa],xilGG_nsn[aa]))
    if lorder == 4:
        Yk1 = YkLL[0]
        Yk2 = YkLL[1]
        Yk3 = YkLL[2]
        Yk4 = YkLL[3]
        Yk5 = YkLL[4]
        Yr1 = np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ rad[maskrp]**5
        Yr2 = np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ rad[maskrp]**3
        Yr3 = np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ rad[maskrp]
        Yr4 = np.fft.ifft(Yk4)/ drad
        Yr4[maskrp] = Yr4[maskrp]/ rad[maskrp]**4
        Yr5 = np.fft.ifft(Yk5)/ drad
        Yr5[maskrp] = Yr5[maskrp]/ rad[maskrp]**2
        xilLL = np.imag(Yr1 + Yr2 + Yr3) + np.real(Yr4 + Yr5)
        # -----------------
        Yk1 = YkLG[0]
        Yk2 = YkLG[1]
        Yk3 = YkLG[2]
        Yk4 = YkLG[3]
        Yk5 = YkLG[4]
        Yr1 = np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ rad[maskrp]**5
        Yr2 = np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ rad[maskrp]**3
        Yr3 = np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ rad[maskrp]
        Yr4 = np.fft.ifft(Yk4)/ drad
        Yr4[maskrp] = Yr4[maskrp]/ rad[maskrp]**4
        Yr5 = np.fft.ifft(Yk5)/ drad
        Yr5[maskrp] = Yr5[maskrp]/ rad[maskrp]**2
        xilLG = np.imag(Yr1 + Yr2 + Yr3) + np.real(Yr4 + Yr5)
        # -----------------
        Yk1 = YkGG[0]
        Yk2 = YkGG[1]
        Yk3 = YkGG[2]
        Yk4 = YkGG[3]
        Yk5 = YkGG[4]
        Yr1 = np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ rad[maskrp]**5
        Yr2 = np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ rad[maskrp]**3
        Yr3 = np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ rad[maskrp]
        Yr4 = np.fft.ifft(Yk4)/ drad
        Yr4[maskrp] = Yr4[maskrp]/ rad[maskrp]**4
        Yr5 = np.fft.ifft(Yk5)/ drad
        Yr5[maskrp] = Yr5[maskrp]/ rad[maskrp]**2
        xilGG = np.imag(Yr1 + Yr2 + Yr3) + np.real(Yr4 + Yr5)
        # -----------------
        Yk1 = YkGG_nsn[0]
        Yk2 = YkGG_nsn[1]
        Yk3 = YkGG_nsn[2]
        Yk4 = YkGG_nsn[3]
        Yk5 = YkGG_nsn[4]
        Yr1 = np.fft.ifft(Yk1)/ drad
        Yr1[maskrp] = Yr1[maskrp]/ rad[maskrp]**5
        Yr2 = np.fft.ifft(Yk2)/ drad
        Yr2[maskrp] = Yr2[maskrp]/ rad[maskrp]**3
        Yr3 = np.fft.ifft(Yk3)/ drad
        Yr3[maskrp] = Yr3[maskrp]/ rad[maskrp]
        Yr4 = np.fft.ifft(Yk4)/ drad
        Yr4[maskrp] = Yr4[maskrp]/ rad[maskrp]**4
        Yr5 = np.fft.ifft(Yk5)/ drad
        Yr5[maskrp] = Yr5[maskrp]/ rad[maskrp]**2
        xilGG_nsn = np.imag(Yr1 + Yr2 + Yr3) + np.real(Yr4 + Yr5)
        
        if debug: print("xill {:10.3e} {:10.3e} {:10.3e} {:10.3e} ".format(xilLL[aa]
                        ,xilLG[aa],xilGG[aa],xilGG_nsn[aa]))
    r2 = (rad*rad)
    r2xilLL = r2 * xilLL
    r2xilLG = r2 * xilLG
    r2xilGG = r2 * xilGG
    r2xilGG_nsn = r2 * xilGG_nsn


    return rad, r2xilLL, r2xilLG, r2xilGG, r2xilGG_nsn



def LyALAEXiangComps(sigr,zred_out,b_LAE=3.0,b_delta=0.707,
              b_Gamma=-0.173,tau_eff=6,
              a_LAE=1,f_LAE=1,f2_LAE=1,debug=False):
    '''
     Returns angular correlation functions from
     LyALAEXillComps_fft output, for density, density-Gamma and Gamma
     fluctuations separately.
    
     ARGUMENTS
      sigr        Comoving width of top-hat shell in r-space (cMpc/h)
      zred_out    Output redshifts used for dGammakCorr arrays (low to high)

     RETURNS
      theta        Angular separation (arcsecs)
      xiaLL        LyA emitter system-system auto-correlation contribution
      xiaLG        LyA emitter system-Gamma cross-correlation contribution
      xiaGG_nsn    Gamma-Gamma auto-correlations without shot noise contribution
      xialGG       Gamma-Gamma auto-correlations with shot noise contribution
    
     COMPATIBILITY: Octave
    
     REQUIREMENTS:
                LyAGetDAng.m, cdenCosparamInit.m
    
    
     AUTHOR: Avery Meiksin
    
     HISTORY:
      04 10 21 Creation date.
    '''



    r,r2xi0LL,r2xi0LG,r2xi0GG,r2xi0GG_nsn = LyALAEXil_fft(zred_out,0,b_LAE=b_LAE,
                        b_delta=b_delta,b_Gamma=b_Gamma,tau_eff=tau_eff,
                        a_LAE=a_LAE,f_LAE=f_LAE,f2_LAE=f2_LAE,debug=debug)
    r,r2xi2LL,r2xi2LG,r2xi2GG,r2xi2GG_nsn = LyALAEXil_fft(zred_out,2,b_LAE=b_LAE,
                        b_delta=b_delta,b_Gamma=b_Gamma,tau_eff=tau_eff,
                        a_LAE=a_LAE,f_LAE=f_LAE,f2_LAE=f2_LAE,debug=debug)
    r,r2xi4LL,r2xi4LG,r2xi4GG,r2xi4GG_nsn = LyALAEXil_fft(zred_out,4,b_LAE=b_LAE,
                        b_delta=b_delta,b_Gamma=b_Gamma,tau_eff=tau_eff,
                        a_LAE=a_LAE,f_LAE=f_LAE,f2_LAE=f2_LAE,debug=debug)

    lenr = len(r)

    xiaLL = np.zeros(lenr)
    xiaLG = np.zeros(lenr)
    xiaGG = np.zeros(lenr)
    xiaGG_nsn = np.zeros(lenr)
    theta = np.zeros(lenr)

    r2 = r * r  # r is comoving in cMpc/h
    if r2[0] ==0.0 : r2[0] = r[1]*r[1]/4.0
    rmax = r[lenr-1]
    #print(r[0],rmax)

    Dang, hubb = LyAGetDang(zred_out) # Dang in Mpc (proper)
    Dang = hubb * (1.0 + zred_out) * Dang # express Dang in cMpc/h
    #print(Dang)
    thetar = r/Dang
    thetas = thetar/4.8481e-6 #theta in arcsec
    maskthp = thetas > 0
    maskthmax = thetas < 10000.0

    maskth = (maskthp * maskthmax)
    lenmth = np.sum(maskth)
    #print(lenmth)
    theta[:lenmth] = thetas[maskth]
    
    # integrate along line of sight using 16-pt gauss-legendre 
    # quadrature
    nglwt = 16
    glxpt = [0.0052995, 0.0277125, 0.0671844, 0.1222978, 0.1910619, 
             0.2709916, 0.3591982, 0.4524937, 0.5475063, 0.6408018,
             0.7290084, 0.8089381, 0.8777022, 0.9328156, 0.9722875,
             0.9947005]
    glwt = [0.013576, 0.031127, 0.047579, 0.062314, 0.074798,
            0.084578, 0.091302, 0.094725, 0.094725, 0.091302,
            0.084578, 0.074798, 0.062314, 0.047579, 0.031127, 0.013576]

    for igl in range(nglwt):
        ugl = sigr * glxpt[igl]
        uth = (ugl * ugl + (Dang*thetar[maskth])**2)**0.5
        ## The iuth line is NOT working.
        #print(lenr,uth,rmax)
        iuth = (lenr*uth/rmax)
        #print('lenr',lenr,rmax,uth[1])
        #print('iuth',iuth[1],np.round(iuth[1]))
        iuth = iuth.astype(int) -1
        #print(len(iuth),lenmth)
        mu = 1.0/(1.0 + (Dang * thetar[maskth]/ugl)**2)**0.5

        L2 = (3 * mu**2 -1.0)/2.0
        L4 = (5 * mu**2 * (7.0 * mu**2 - 6.0)+3.0)/8.0

        aa = 49
        if (igl ==2) & debug: print('L4',L2[aa],L4[aa])
        
        xiaLL[:lenmth] = xiaLL[:lenmth] + glwt[igl] * (r2xi0LL[iuth] + \
                        L2 * r2xi2LL[iuth] + L4 * r2xi4LL[iuth])/r2[iuth]

        xiaLG[:lenmth] = xiaLG[:lenmth] + glwt[igl] * (r2xi0LG[iuth] + \
                        L2 * r2xi2LG[iuth] + L4 * r2xi4LG[iuth])/r2[iuth]

        xiaGG[:lenmth] = xiaGG[:lenmth] + glwt[igl] * (r2xi0GG[iuth] + \
                        L2 * r2xi2GG[iuth] + L4 * r2xi4GG[iuth])/r2[iuth]

        xiaGG_nsn[:lenmth] = xiaGG_nsn[:lenmth] + glwt[igl] * (r2xi0GG_nsn[iuth] + \
                        L2 * r2xi2GG_nsn[iuth] + L4 * r2xi4GG_nsn[iuth])/r2[iuth]
    
    return theta, xiaLL, xiaLG,  xiaGG, xiaGG_nsn


def LyAGetDang(zred):
    '''
     ARGUMENTS
      zred         Redshift.
    
    
     RETURNS
      DAng         Proper angular distance (Mpc).
      h            H0/(100 km/s/Mpc)
    
     COMPATIBILITY: Matlab
    
     REQUIREMENTS: 
                cdenCosparamInit.m called previously
    
     AUTHOR: Avery Meiksin
    
     HISTORY:
      05 10 21 Creation date. (After SNRGetDLum.m)


    '''
    z_min = 0.0
    tol = 1.0e-5
    D,kk = integrate.quad(SNRdyfncdz,z_min,zred,epsrel = tol)
    D = 2.99792458e3/h * D
    Dang = D / (1.0 + zred)

    return Dang, h

def SNRdyfncdz(zred):
    #  Integrand for SNRDLum to compute luminosity distance.
    Efnc = (omega[0]*(1.0 + zred)**3 + omega[1])**0.5
    dydz = 1./ Efnc
    return dydz
