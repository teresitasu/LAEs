import numpy as np
#from scipy import interpolate
from scipy.interpolate import PchipInterpolator as pchip
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import pandas as pd
from scipy.special import gamma,gammaincc,gammainc,gammainccinv,gammaincinv




# Define Cosmology
om_m = 0.3111
om_v = 1.0 - om_m
om_k = 0.0
om_bh2 = 0.02242
h = 0.6766
an = 0.9665
sigma8 = 0.8102

omega = [om_m, om_v, 0., om_bh2]

fk = np.logspace(-4.0,3.0,701)

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

def LyAGetNeffGal(zred,lmin):
    '''
     Computes effective comoving number density of galaxies from
     galaxy luminosity function of Bouwens et al. (2015) (ApJ 803:34)

     ARGUMENTS
     zred     List of redshifts
     lmin     Minimum Lum/Lums for integration

     RETURNS
      neff     effective number density of galaxies as function or redshift (cMpc^-3)
      avgL    Mean emissivity
      avgL2  Mean L^2-emissivity
      

     COMPATIBILITY: Matlab, Octave

     REQUIREMENTS


     AUTHOR: Avery Meiksin

     HISTORY:
      22 12 17 Creation date.
    '''


    e = np.exp(1.0)
    zp1 = 1.0 + zred
    LumLs = 1.80e28
    alpha = -1.87 + 0.1 * (6.0 - zred)
    #print(alpha)
    ss = alpha <= -1.9
    alpha[ss] = -1.9
    #alpha = np.max([alpha,-1.9],axis=-)
    g2 = gamma(alpha+2.0) * gammaincc(alpha+2.0,lmin)
    g3 = gamma(alpha+3.0) * gammaincc(alpha+3.0,lmin)
    #print(g2)
    coeff1 = 2.5 * np.log10(e)*g2*LumLs
    coeff2 = 2.5 * np.log10(e)*g3*LumLs**2
    phis = 0.44e-3 * 10**(0.28*(6.0 - zred))
    fesc = 1.8e-4 * zp1**3.4
    avgL = coeff1 *fesc*phis
    avgL2 = coeff2 *fesc * fesc * phis
    neff = (avgL *avgL) * (1.0/ avgL2)
    #neff = neff

    return neff, avgL, avgL2

def LyAGetNeffQSO_KWH19(zred,M1450min,M1450max,imod):
    '''
     Computes effective comoving number density of QSOs from
     QSO luminosity function

     ARGUMENTS
     zred           List of redshifts
     M1450min   Minimum 1450 absolute magnitude
     M1450max  Maximum 1450 absolute magnitude
     imod          1, 2 or 3

     RETURNS
      neff     effective number density of QSOs as function or redshift (cMpc^-3)
      avgL    Mean comoving emissivity (erg/s/Hz/cMpc^3)
      avgL2  Mean comoving L^2-emissivity [(erg/s/Hz)^2/cMpc^3]
      

     COMPATIBILITY: Matlab, Octave

     REQUIREMENTS


     AUTHOR: Avery Meiksin

     HISTORY:
      11 09 19 Creation date. (After LyAGetNeffQSO.m)
    '''
    lenM = 1001
    lenz = len(zred)
    M1450 = np.linspace(M1450min,M1450max,lenM)
    ### ! Missing here next two lines!
    dM1450 = np.gradient(M1450)
    Phi = LyAGetQSOLF_KWH19(M1450,zred,imod)
    # Convert to luminosity at Lyman edge (KWH19 Eqs.(20) and (21)
    LumL = ((912./1450.)**0.61)*10**(-0.4*(M1450 - 51.60))
    #print(Phi.shape,LumL.shape)
    LumLmat = np.ones((lenz,lenM))* (LumL)
    #print(LumL)
    #print(np.transpose(LumLmat) * (LumL))
    #print((np.transpose(LumLmat) * Phi).shape,dM1450.shape)
    avgL = np.dot(LumLmat * np.transpose(Phi), dM1450)
    avgL2 = np.dot((LumLmat * LumLmat) * np.transpose(Phi), dM1450)
    neff = (avgL * avgL)*(1./ avgL2)
    #avgL = avgL
    #avgL2 = avgL2'
    #neff = neff'
    return neff, avgL, avgL2

def LyAGetQSOLF_KWH19(M1450,zred,imod):
    '''
     Returns QSO luminosity function
     Uses Kulkarni et al. (2019) MNRAS 488:1035

     ARGUMENTS
     M1450  Absolute magnitude at 1450A
     zred     List of redshifts (CURRENTLY FOR A SINGLE REDSHIFT)
     imod    1, 2 or 3

     RETURNS
      Phi       2D array containing Phi(M1450,z) (row,col)
                 (Phi(M1450,z)dM1450 = comoving number density in Mpc^-3 of QSOs at z
               with M1450 absolute magnitude between M1450 and M1450 + dM1450)
      

     COMPATIBILITY: Matlab, Octave

     REQUIREMENTS


     AUTHOR: Avery Meiksin

     HISTORY:
      11 09 19 Creation date. (From LyAGetQSOLF.m)
    '''
    #M1450 = M1450
    lenM = len(M1450)
    lenz = len(zred)
    M1450mat = np.ones((lenz,lenM))*M1450
    Phi = np.zeros((lenM,lenz))
    zp1 = 1 + zred
    # Define Chebyshev poloynomials of first kind, orders 0 to 3
    T0 = 1.0
    T1 = zp1
    T2 = 2*zp1*zp1 - 1.0
    T3 = 4*zp1*zp1*zp1 - 3.0*zp1
    if(imod<1) or (imod>3):
        print('LyAGetQSOLF_KWH19: invalid imod')
        return None
    if (imod == 1):
        c0 = np.array([-7.798,1.128,-0.120])
        c1 = np.array([-17.163,-5.512,0.593,-0.024])
        c2 = np.array([-3.223,-0.258])
        c3 = np.array([-2.312,0.559,3.773,141.884,-0.171])
        phis = 10**(np.dot(c0,[T0,T1,T2]))
        Ms = np.dot(c1,[T0,T1,T2,T3])
        Msmat = np.ones(lenM)*Ms
        alpha = np.dot(c2,[T0,T1])
        zeta = np.log10(zp1/ (1.0 + c3[2]));
        beta = c3[0] + c3[1]/ (10**(c3[3]*zeta) + 10**(c3[4]*zeta))
        alphap1 = 1 + alpha
        betap1 = 1 + beta
        x = 10**(0.4*(np.transpose(M1450mat) - Msmat))
        phismat = np.ones((lenM,lenz))*phis
        alphap1mat = np.ones((lenM,lenz))*alphap1
        betap1mat = np.ones((lenM,lenz))*betap1
        Phi = phismat/(x**alphap1mat + x**betap1mat)

    if (imod==2):
        c0 = np.array([-7.432,0.953,-0.112])
        c1 = np.array([-15.412,-6.869,0.778,-0.032])
        c2 = np.array([-2.959,-0.351])
        c3 = np.array([-2.264,0.530,2.379,12.527,-0.229])
        phis = 10**(np.dot(c0,[T0,T1,T2]))
        Ms = np.dot(c1,[T0,T1,T2,T3])
        Msmat = np.ones(lenM)*Ms
        alpha = np.dot(c2,[T0,T1])
        zeta = np.log10(zp1/ (1 + c3[2]))
        beta = c3[0] + c3[1]/ (10**(c3[3]*zeta) + 10**(c3[4]*zeta))
        alphap1 = 1 + alpha
        betap1 = 1 + beta
        x = 10**(0.4*(np.transpose(M1450mat) - Msmat))
        phismat = np.ones((lenM,lenz))*phis
        alphap1mat = np.ones((lenM,lenz))*alphap1
        betap1mat = np.ones((lenM,lenz))*betap1
        Phi = phismat/(x**alphap1mat + x**betap1mat)
    if (imod==3):
        c0 = np.array([-6.942,0.629,-0.086])
        c1 = np.array([-15.038,-7.046,0.772,-0.030])
        c2 = np.array([-2.888,-0.383])
        c3 = np.array([-1.602,-0.082])
        phis = 10**(np.dot(c0,np.array([T0,T1,T2])))
        Ms = np.dot(c1,np.array([T0,T1,T2,T3]))
        #print(Ms)
        Msmat = np.ones((lenM,lenz))*Ms
        alpha = np.dot(c2,np.array([T0,T1]))
        beta = c3[0] + c3[1]*zp1
        alphap1 = 1.0 + alpha
        betap1 = 1.0 + beta
        #print(M1450mat.shape,Msmat.shape)
        x = 10**(0.4*(np.transpose(M1450mat) - Msmat))
        phismat = np.ones((lenM,lenz))*phis
        alphap1mat = np.ones((lenM,lenz))*alphap1
        betap1mat = np.ones((lenM,lenz))*betap1
        #print(phismat.shape,alphap1mat.shape,betap1mat.shape)
        Phi = phismat/(x**alphap1mat + x**betap1mat)
    return Phi

def LyASolvedGammakCorrAsymptoticApprox(k,Pk,Gammai=None,aj=None,aS=None,
                            ab=None,zi=None,
                            zf=None,bfrac=None,
                            bet=None,cfrac=None,iQmod=None,
                            bj=None,bnH=None,baA=None,
                            bLLS=None,tau_Q=None,M1450min=None,
                            M1450max=None,bG=None,tau_G=None,
                            e24=None):
    '''
     Returns deltaGammaCorr based on asympytotic limits alone: power spectrum of
     photoionization rate Gamma fluctuations, proportional to dark matter density
     power spectrum. Adopts conventional units of the dark matter power spectrum
     as comoving volume units (Mpc/h)^3. Allows for QSO bias from
     Laurent et al. (2017).

     ARGUMENTS
     Gammai      Initial guess on ionization rate (for LyASolveEmissMG.m)
     aj          Source comoving emissivity spectral power-law exponent
     aS          Source comoving emissivity evolution power-law exponent
     ab          Value of a_b for bias factor evolution bj ~ (1+z)^ab;
                    compute internally for ab < 0; set ab = 1 for high redshift (large aeff_LLS)
     om_m        Omega_m for all matter
     h           H0/ (100 km/s/Mpc)
     zi          Ionization redshift
     zf          Final redshift
     k           Single comoving wavenumber (h/ Mpc)
     Pk          Value of power spectrum P(k) at k at z = 0 (comoving Mpc^3/ h^3)
     bfrac       Fraction of baryons in diffuse component of IGM
     bet         Exponent for power law HI distribution of Lyman Limits Systems
     cfrac       Fraction of intergalactic attenuation due to Lyman Limit Systems
     iQmod       QSO evolution model
     bj          Source bias (set to <=0 to override bQ and bG and use -bj instead)
     bnH         Hydrogen density bias
     baA         Radiative recombination rate bias
     bLLS        Bias factor for LLS contribution to attenuation
     tau_Q       Lifetime of QSO sources (Myr)
     lgQbolmin   Log10 minimum QSO bolometric luminosity
     lgQbolmax   Log10 maximum QSO bolometric luminosity
     bG          Galaxy bias factor
     tau_G       Lifetime of galaxy sources (Myr)                          %

     RETURNS
      e24            Comoving emissivity coefficient (1e24 erg/s/Hz/Mpc^3)
      zred           Array of redshifts from zf to zi
      dGammaCorr     Power spectrum of fluctuation in photoionization rate Gamma from zf to zi
      dGammaCorr_nsn Power spectrum of fluctuation in photoionization rate Gamma from zf to zi without shotnoise
      dGammaCorr_sn  Power spectrum of fluctuation in photoionization rate Gamma from zf to zi for shotnoise term alone
      Gamma          Mean metagalactic photoionization rate
      aeff_d         Array of effective diffuse attenuation coefficient, matching zred array
      aeff_LLS       Array of effective LLS attenuation coefficient, matching zred array
      S              Source q*bj - (2*bnH + baA) array, matching zred array

     COMPATIBILITY: Matlab, Octave


     AUTHOR: Avery Meiksin

     HISTORY:
      04 06 18 Creation date. (Adapted from LyASolvedGammakCorrSS.m.)
      05 07 18 Add lgQbolmin and lgQbolmax arguments
      15 08 18 Modify for phi .ne. 1
      24 08 18 Add factor convcomov to convert 1/ neff to comoving (Mpc/ h)^3
      04 09 18 Add option for Lorentzian power spectrum approximation
      21 09 18 Fix faulty denom_sn; allow for evolution in comoving neff through alpha_n
      19 11 21 Fix asymptotic values to steady-state case for high k (for ab = 1) and low k
      02 12 21 Modified to correct galaxy contribution to match required e24
    '''



    #[neff_Q,avgL_Q,avgL2_Q] = LyAGetNeffQSO_KWH19(zred_out,M1450min,M1450max,iQmod);
    if Gammai == None:
        Gammai = 0.249e-12
    if aj == None:
        aj = 1.8
    if aS == None:
        aS = 0.8
    if ab == None:
        # turn off delta_DM evolution for high z (high chiLLS)
        ab = 1 
    if zi == None:
        zi = 8.0
    if zf == None:
        zf = 4.8
    if bfrac == None:
        bfrac = 0.0
    if bet == None:
        bet = 1.2
    if cfrac == None:
        cfrac = 1.0
    if iQmod == None:
        iQmod = 3
    if bj == None:
        bj = 1.0
    if bnH == None:
        bnH = 1.0
    if baA == None:
        baA = -0.5
    if bLLS == None:
        bLLS = 1.0
    if tau_Q == None:
        tau_Q = 10.0
    if M1450min == None:
        M1450min = -31.0
    if M1450max == None:
        M1450max = -21.0
    if bG == None:
        bG = 3.0
    if tau_G == None:
        tau_G = 100.0
    if e24 == None:
        e24 = 29.329250318837154

    mpc = 3.0856e24
    h_LF = 0.7
    lmin = 0.01
    bGamLLS = 1.0 - bet
    om_v = 1.0 - om_m
    deltk_0 = Pk**0.5
    gg_0,ll = cdenGrowth(1.0)
    #print(gg_0)
    ##### MISSING LyASolveEmissMG
    #flags = True
    Gamma_MG,xHI,aeffd_MG,aeffLLS_MG,otsfac_MG,zred_MG = np.loadtxt('EmissMG.txt',unpack=True)
    '''
    if flags:
        print('Missing')
    else:
        zred_MG = np.linspace(zf,zi,lzred_MG)
        zp1_MG = 1.0 + zred_MG
        Gamma_i = Gammai * np.ones(1.0,lzred_MG)
        e24,Gamma_MG,xHI,aeffd_MG,aeffLLS_MG,otsfac_MG,err = LyASolveEmissMG(Gamma_i,aj,aS,om_m,om_v,h,bfrac,bet,cfrac,zred_MG)
    '''
    ##### MISSING LyASolveEmissMG

    aeff_MG = aeffd_MG + aeffLLS_MG
    zred = zred_MG
    lzred = len(zred)

    zp1 = 1.0 + zred
    zp1_3 = zp1**3
    #factor to convert shot noise from 1/ neff to comoving (Mpc/h)^3
    convcomov = h_LF**3

    aexp_flip = 1.0/ zp1
    gg,ll = cdenGrowth(aexp_flip[::-1])
    deltk = np.array(gg)[::-1]*deltk_0/ gg_0
    bQ = 0.278*zp1*zp1 + 0.57 #from BOSS
    Gamma = Gamma_MG
    aeff = aeff_MG
    aeff_d = aeffd_MG
    aeff_LLS = aeffLLS_MG
    aeff_wtot = aeff_d - (bGamLLS * aeff_LLS)
    otsfac = otsfac_MG
    Hdc = h*(om_m*zp1**3 +om_v)**0.5/(2.99792458e3*mpc)
    aTot = aeff + Hdc * otsfac

    neff_Q,avgL_Q,avgL2_Q = LyAGetNeffQSO_KWH19(zred,M1450min,M1450max,iQmod)
    eL24 = e24/ zp1**aS
    #print('bG',bG)
    if (bG > 0):
        neff_G,avgL_G,avgL2_G = LyAGetNeffGal(zred,lmin)
        fixfescG = (1e24*eL24 - avgL_Q)/ avgL_G
        avgL_G = fixfescG*avgL_G
        avgL2_G = fixfescG*fixfescG*avgL2_G
    else:
        fixQ = 1e24*eL24/ avgL_Q
        avgL_Q = fixQ*avgL_Q
        avgL2_Q = fixQ*fixQ*avgL2_Q

    avgL = eL24*1e24
    if (bj <= 0):
        bj = -bj
    else:
        if (bG > 0):
            bj = (bQ*avgL_Q + bG*avgL_G)/ avgL
        else:
            bj = bQ
            avgL = eL24*1e24
    q = 0.999e-16 * e24 * zp1**(3 - aS)/ (3 + aj)
    shot = np.zeros(lzred)
    shot_Q = np.zeros(lzred)
    avgavgL2_Q = np.zeros(lzred)
    if (bG > 0):
        shot_G = np.zeros(lzred)
        avgavgL2_G = np.zeros(lzred)

    jdf = q/ (Gamma*mpc)
    if ((cfrac < 1.0) or (abs(bGamLLS) > 0)):
        q = jdf/ aeff_wtot
        S = bj*q - (2.0*bnH + baA)*(aeff_d/ aeff_wtot) - \
                     bLLS*(aeff_LLS/ aeff_wtot)
    else:
        q = jdf/ aeff
        S = bj*q - bLLS*(aeff_LLS/ aeff)
    S = S * deltk

    if (ab < 0.0):
        ab = np.gradient(np.log(bj))/ np.gradient(np.log(zp1))
    
    chiH = aeff/ Hdc
    chiH_d = aeff_d/ Hdc
    chiH_LLS = aeff_LLS/ Hdc
    #chiHb = chiH + otsfac #otsfac is beta_H
    chiHb = jdf/ Hdc #phi*(chi+zeta)
    chiHdL = chiH_d - bGamLLS*chiH_LLS
    gamH = chiHb - ab + 0.5
    kap = (k * h )/ (mpc*Hdc)
    #zp1i = zp1(lzred)
    #xi = zp1./ zp1i;
    # kap -> 0 limit
    denom = bGamLLS*chiH + chiHb + 1.0 - ab #otsfac is beta_H
    pi2 = np.pi / 2.0
    #Use Lorentzian formulation for non-shotnoise signal
    kaps = (pi2/ zp1) * denom
    fkap = -chiHdL/ (denom*(1.0 + (kap/ kaps)**2)**0.5)
    S = S*fkap
    S2ns = S**2
    # Add shotnoise term in E-deS limit
    ############
    ##Use Lorentzian formulation for shot noise
    #Use extended Lorentzian formulation for shot noise        
    if (bG > 0):
      neff = (avgL*avgL)/ (avgL2_Q + avgL2_G)
    else:
      neff = (avgL*avgL)/ avgL2_Q
    
    # comoving neff evolution exponent
    alpha_n = -np.gradient(np.log(neff))/ np.gradient(np.log(zp1))
    HMyr = 3.15e20*h*(om_m*zp1**3+om_v)**0.5/ mpc
    tau_S = tau_Q ### assume tau_Q dominates
    denom_sn = (1.0 - bet)*chiH + chiHb + 0.75 - 0.5*alpha_n
    # Use sliding scale-length
    kaps_sn = 2.0 * denom_sn/ (HMyr*tau_S) - \
            ((chiHb - 0.5)/(pi2*pi2) + (1.0 - bet)*chiH)**2
    kaps_sn = (pi2/ zp1)* kaps_sn**0.5
    kaps_sn = kaps_sn/ (1.0 + pi2*((chiHb - 0.5)/(pi2*pi2) +\
             (1.0 - bet)*chiH)/ (kap*zp1))
    AA = chiHb*chiHb
    AA = AA/ (2.0 * denom_sn/ (HMyr*tau_S) - \
             ((chiHb - 0.5)/(pi2*pi2)+(1.0 - bet)*chiH)**2)
    fkap_sn = (AA/ (1.0 + (kap/ kaps_sn)**2))**0.5
    S_sn = (convcomov**0.5)*fkap_sn/ neff**0.5
    shot = S_sn**2
    S2 = S2ns + shot
    dGammaCorr = S2
    dGammaCorr_nsn = S2ns
    dGammaCorr_sn = shot


    return e24,zred,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S


def cdenPowspInitFile(ips=3,Tk_file=None):

    if Tk_file == None:
        Tk_file = './tkO3111L6889h6766N2242_eh_PLANCK18.dat'
    ng = 1e6
    boxsize = 1000
    tol = 1e-6
    fkl = 2*np.pi/ boxsize
    fku = 0.5*ng*fkl
    lgfkl = np.log(fkl)
    lgfku = np.log(fku)
    #Uses Bunn & White (1997) COBE normalization.
    tilt = an - 1.0
    om_b = om_bh2/ (h*h)
    if (om_v > 0.0):
        dh=1.94e-5*(om_m**(-0.785-0.05*np.log(om_m)))*np.exp(-0.95*tilt-0.169*tilt*tilt)
    else:
        dh=1.95e-5*(om_m**(-0.35-0.19*np.log(om_m)-0.17*tilt))*np.exp(-tilt-0.14*tilt*tilt)

    pnorm = 2.0*np.pi*np.pi*dh*dh*(2997.9)**(3.0 + an)
    Gamma = om_m*h/ (np.exp(om_b*(1.0 + 1.0/ om_m)))
    if (ips == 3):
        Tk = np.loadtxt(Tk_file)
        #fk = Tk_table[:,0]
        #fkmin = np.min(Tk.fk)
        #fkmax = np.max(Tk.fk)
        #tk = Tk_table[:,1]
    else:
        Tk = 0
    print('Bunn & White COBE pnorm: {:.5e}'.format(pnorm))
    # Now normalize to desired sigma8.
    xf = 8.0
    ismooth = 1 #spherical top hat
    ## pnorm,Gamma,an,Tk,xf,ismooth,ips
    sigma8_BW,garbage = integrate.quad(cdenFNorm,lgfkl,lgfku,
            args=(pnorm,Gamma,Tk,xf,ismooth,ips),epsrel=tol)
    #sigma8_BW = quad(@(u)cdenFNorm(u,pnorm,Gamma,an,Tk,xf,ismooth,ips),lgfkl,lgfku,tol)
    #% TEST for linear integration
    #fk=linspace(lgfkl,lgfku,1000000);
    #fk(1)
    #fk(1000000)
    #dk = (lgfku-lgfkl)/1000000;
    #u = cdenFNorm(fk,pnorm,Gamma,an,Tk,xf,ismooth,ips);
    #sigma8_BW = dk*trapz(u);
    #%
    sigma8_BW = sigma8_BW**0.5
    pnorm = pnorm*(sigma8/ sigma8_BW)*(sigma8/ sigma8_BW)
    
    return pnorm, Gamma, Tk

def cdenFNorm(lgfk,pnorm,Gamma,Tk,xf,ismooth,ips):
    '''
    ******************************************************
     Computes normalization integrand for power spectrum.
    
     ARGUMENTS
      lgfk       Array of (logarithmically-spaced) wavenumber.
      pnorm      Normalization of matter power spectrum.
      Gamma      Curvature of matter power spectrum.
      an         Spectral index of matter power spectrum.
      Tk         Transfer function for matter power spectrum.
      xf         Filtering scale (inverse units of dlfkp).
      ismooth    Method for smoothing fluctuations.
      ips        Method of computing matter power spectrum.
    
     COMPATIBILITY: Matlab(?), Octave
    
     REQUIREMENTS:
             cdenPowspInit.m, cdenPowsp.m
    
     AUTHOR: Avery Meiksin
    
     HISTORY:
      19 04 16 Creation date. (After FNORM subfunction in lss/src/Delta2.f.)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    '''
    fk = np.exp(lgfk)
    #fk = lgfk; %TEST for linear integration
    y = xf*fk
    if (ismooth == 0):
        W = 1.0
    if (ismooth == 1):
        W = 3.0*(np.sin(y) - y*np.cos(y))/ y**3
    if (ismooth == 2):
        W = 0.0
        maskl = y < 10
        W[maskl] = 1./ np.exp(y[maskl]**2)

    if (ismooth == 3):
        W = 1./ (1.0 + y**2)
    PS = cdenPowsp(fk,pnorm,Gamma,Tk,ips)
    return fk**3 * PS * W**2/ (2.0 * np.pi**2)


def cdenPowsp(fk,pnorm,Gamma,Tk,ips):
    '''
    #********************
    # Computes matter power spectrum.
    # Assumes k in Tk file is in conventional comoving h/ Mpc.
    # Conventionally pnorm is for P(k) in units (comoving Mpc/h)^3,
    # so PS is returned in units (comoving Mpc/h)^3.
    #
    # ARGUMENTS
    #  fkp        Array of wavenumbers.
    #  pnorm      Normalization of power spectrum.
    #  Gamma      Curvature of power spectrum.
    #  an         Tilt of power spectrum.
    #  Tk         Transfer function.
    #  ips        Method for computing power spectrum.
    #
    # COMPATIBILITY: Matlab(?), Octave
    #
    # AUTHOR: Avery Meiksin
    #
    # HISTORY:
    #  19 04 16 Creation date. (After lss/src/powsp.f.)
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    '''
    if (ips == 1):
        # use B&E
        a = 6.4/ Gamma
        b = 3.0/ Gamma
        c = 1.7/ Gamma
        PS = 1.0 + (fk*(a + c*c*fk) + (b*fk)**1.5)**1.13
        PS = pnorm*fk**an/ PS**(2./ 1.13)

    if (ips == 2):
        # use BBKS
        q = fk/ Gamma
        PS = 1.0 + q*(3.89+q*(16.1*2 + q* (5.46**3 + q * 6.71**4)))
        PS = np.log(1. + 2.34*q)/ PS**0.25/ (2.34*q)
        PS = pnorm * (fk**an) * PS**2

    if(ips == 3):
        ggg = interp1d(Tk[:,0],Tk[:,1],kind='cubic',fill_value="extrapolate")
        tki = ggg(fk)
        PS = pnorm*(fk**an)*tki*tki
        PS = np.array(PS)
        maskl = fk < np.min(Tk[:,0])
        masku = fk > np.max(Tk[:,0])
        q = fk[maskl]/ Gamma
        PS[maskl] = 1.0 + q*(3.89+q*(16.1*16.1+q*(5.46*5.46*5.46+q*6.71*6.71*6.71*6.71)))
        PS[maskl] = np.log(1. + 2.34*q)/ PS[maskl]**0.25/ (2.34*q)
        PS[maskl] = pnorm*fk[maskl]**an*PS[maskl]*PS[maskl]
    
        q = fk[masku]/ Gamma
        PS[masku] = 1.0 + q*(3.89+q*(16.1*16.1+q*(5.46*5.46*5.46+q*6.71*6.71*6.71*6.71)))
        PS[masku] = np.log(1. + 2.34*q)/ PS[masku]**0.25/ (2.34*q)
        PS[masku] = pnorm*fk[masku]**an*PS[masku]*PS[masku]
    
    return PS


def LyASolvedGammakCorr(zred,Gammai=None,aj=None,aS=None,
                    ab=None,zi=None,
                    zf=None,bfrac=None,
                    bet=None,cfrac=None,iQmod=None,
                    bj=None,bnH=None,baA=None,
                    bLLS=None,tau_Q=None,M1450min=None,
                    M1450max=None,bG=None,tau_G=None,
                    e24=None):

    if Gammai == None:
        Gammai = 0.249e-12
    if aj == None:
        aj = 1.8
    if aS == None:
        aS = 0.8
    if ab == None:
        # turn off delta_DM evolution for high z (high chiLLS)
        ab = 1 
    if zi == None:
        zi = 8.0
    if zf == None:
        zf = 4.8
    if bfrac == None:
        bfrac = 0.0
    if bet == None:
        bet = 1.2
    if cfrac == None:
        cfrac = 1.0
    if iQmod == None:
        iQmod = 3
    if bj == None:
        bj = 1.0
    if bnH == None:
        bnH = 1.0
    if baA == None:
        baA = -0.5
    if bLLS == None:
        bLLS = 1.0
    if tau_Q == None:
        tau_Q = 10.0
    if M1450min == None:
        M1450min = -31.0
    if M1450max == None:
        M1450max = -21.0
    if bG == None:
        bG = 3.0
    if tau_G == None:
        tau_G = 100.0
    if e24 == None:
        e24 = 29.329250318837154

    dkCorr = []
    dkCorr_nsn = []
    dkCorr_sn = []
    ctr = 0
    for k,Pk in zip(fk,PS):
        e24,zredout,dGammaCorr,dGammaCorr_nsn,dGammaCorr_sn,Gamma,aeff_d,aeff_LLS,S = LyASolvedGammakCorrAsymptoticApprox(k,Pk,
            Gammai=Gammai,aj=aj,aS=aS,
                    ab=ab,zi=zi,
                    zf=zf,bfrac=bfrac,
                    bet=bet,cfrac=cfrac,iQmod=iQmod,
                    bj=bj,bnH=bnH,baA=baA,
                    bLLS=bLLS,tau_Q=tau_Q,M1450min=M1450min,
                    M1450max=M1450max,bG=bG,tau_G=tau_G,
                    e24=e24)
        masku = zredout >= zred

        dkCorr.append(dGammaCorr[masku][0])
        dkCorr_nsn.append(dGammaCorr_nsn[masku][0])
        dkCorr_sn.append(dGammaCorr_sn[masku][0])
        if ctr%10 == 0: print(ctr,fk.size,k)
        ctr +=1

    return fk, np.array(dkCorr),np.array(dkCorr_nsn),np.array(dkCorr_sn),PS


