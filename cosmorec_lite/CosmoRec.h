//===========================================================================================================
//
// Author: Jens Chluba
// first implementation: April 2010
// last modification   : July  2014
// purpose: to access main CosmoRec code
//
//===========================================================================================================
// 23.07.2014: added possibility to change A2s1s two-photon decay rate
// 22.07.2014: new He-diffusion runmodes (runpars[1]== 4, 5, 6)
// 17.07.2014: included HeI-diffusion correction (most of the coding done 02-05.2011)
// 16.07.2014: default runmode for batch version with Diff_iteration_max=1 and nShellsHeI=2
// 17.06.2013: small cosmetic changes
// 18.06.2012: added C-version to communicate Hubble pointer
// 08.06.2012: added version to communicate Hubble
// 23.03.2012: added runmode with helium radiative transfer module;
// 16.02.2012: added module to load precomputed recombination history (useful for runs with Camb)

#ifndef COSMOREC_H
#define COSMOREC_H

#ifdef __cplusplus

//===========================================================================================================
// CosmoRec for run from console
//===========================================================================================================
int CosmoRec(int narg, char *args[]);

//===========================================================================================================
// to call CosmoRec with a list of parameters in batch mode (e.g. when calling it from CosmoMC).
//
// runmode == 0: CosmoRec run with diffusion
//            1: CosmoRec run without diffusion
//            2: Recfast++ run (equivalent of the original Recfast version)
//            3: Recfast++ run with correction function of Chluba & Thomas, 2010
//
// On entry, the array z_arr should contain the redshifts at which Xe and Te are required. nz
// determines the number of redshift points. Xe_arr and Te_arr will contain the solution on exit.
//
// Furthermore, runpars[0] defines the dark matter annihilation efficiency in eV/s.
// runpars[1] switches the accuracy of the recombination model:
//
// runpars[1]==-1: closest equivalent of 'HyRec' case (Haimoud & Hirata, 2010)
// runpars[1]== 0: default setting
// runpars[1]== 1: 2g for n<=4 & Raman for n<=3
// runpars[1]== 2: 2g for n<=8 & Raman for n<=7
// runpars[1]== 3: 2g for n<=8 & Raman for n<=7 + Helium feedback up to n=5
// runpars[1]== 4: default setting              + Helium radiative transfer
// runpars[1]== 5: 2g for n<=4 & Raman for n<=3 + Helium radiative transfer up to n=3
// runpars[1]== 6: 2g for n<=4 & Raman for n<=3 + Helium radiative transfer up to n=5 (full setting)
//
// The value of runpars[1] is only important for runmode 0 & 1.
//
// CosmoRec can be set to output the ionization history & electron temperature, solutions for
// the populations of the resolved levels, and some information about the Cosmology (usually
// directed into the directory "./outputs/"). This can cause problems when using CosmoRec
// with a parallel CosmoMC, however, it possible to control the outputs through runpars[2]:
//
// runpars[2]==0: don't write out anything (default).
// runpars[2]==1: write out only the recombination history.
// runpars[2]==2: write out the recombination history, and the cosmology.
// runpars[2]==3: write out the recombination history, populations, and the cosmology.
//
// runpars[3]>=0: use given number to set const_HI_A2s_1s. For ==0 default value is used.
//
//===========================================================================================================
int CosmoRec(const int runmode, const double runpars[5],
             const double omegac, const double omegab,
             const double omegak, const double Nnu,
             const double h0, const double tcmb, const double yhe,
             const int nz, double *z_arr, double *Xe_arr, double *Te_arr,
             const int label);

 #endif /* __cplusplus */

//===========================================================================================================
// Wrap the C++ Fortran routine to be allow calling from Fortran. Arguments are as above.
// Added 06.03.2011 (Richard Shaw)
//===========================================================================================================
// 18.06.2012: Added version to communicate Hubble using pointer
// 08.06.2012: Added version to communicate Hubble
//===========================================================================================================
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

    void cosmorec_calc_cpp_(const int * runmode, const double * runpars,
                            const double * omega_c, const double * omega_b, const double * omega_k,
                            const double * num_nu, const double * h0,
                            const double * t_cmb, const double * y_he,
                            double * za_in, double * xe_out, double * tb_out,
                            const int * len, const int* label);

    void cosmorec_calc_h_cpp_(const int * runmode, const double * runpars,
                              const double * omega_c, const double * omega_b, const double * omega_k,
                              const double * num_nu, const double * h0,
                              const double * t_cmb, const double * y_he,
                              const double * z_Hz, const double * Hz, const int * nz,
                              double * za_in, double * xe_out, double * tb_out,
                              const int * len, const int* label);

    void cosmorec_calc_hptr_cpp_(const int * runmode, const double * runpars,
                                 const double * omega_c, const double * omega_b, const double * omega_k,
                                 const double * num_nu, const double * h0,
                                 const double * t_cmb, const double * y_he,
                                 double (* Hptr)(const double * z),
                                 double * za_in, double * xe_out, double * tb_out,
                                 const int * len, const int* label);

    void evaluate_TM(double z, double Xe, double fHe, double rho, double Tg, double Hz, double * drho_dt);

#ifdef __cplusplus
}
#endif /* __cplusplus */

//===========================================================================================================
// to cleanup after finishing with CosmoRec
//===========================================================================================================
void cleanup_CosmoRec();

//===========================================================================================================
#endif
