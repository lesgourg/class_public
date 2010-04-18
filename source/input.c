/** @file input.c 
 * Julien Lesgourgues, 18.04.2010    
 */

#include "input.h" 

int input_init(
		 struct background *pba,
		 struct thermo *pth,
		 struct perturbs *ppt,
		 struct bessels * pbs,
		 struct transfers *ptr,
		 struct primordial *ppm,
		 struct spectra *psp,
		 struct output *pop 
	       ) {

  double h;

  /** - assign values to background cosmological parameters */

  h=0.7;
  /* H0 in Mpc^{-1} = h / 2999.7 */
  pba->H0 = h / 2999.7;
  /* cp.Omega0_g = 2.3812e-5/h/h*pow(2.726/2.7,4.); */
  pba->Omega0_g = 5.05064e-5;
  /* cp.Omega0_nur = 0.*cp.Omega0_g*0.68; */
  pba->Omega0_nur = 3.486993e-5;
  pba->Omega0_b = 0.05;
  pba->Omega0_cdm = 0.23;
  pba->Omega0_lambda = 1.-pba->Omega0_b-pba->Omega0_cdm-pba->Omega0_g-pba->Omega0_nur;
  pba->Omega0_de = 0.;

  /** - assign values to thermodynamics cosmological parameters */

  pth->Tcmb=2.726;
  pth->YHe=0.25;
  pth->reio_parametrization=reio_camb;
  pth->reio_z_or_tau=reio_tau;
  pth->z_reio=10.;
  pth->tau_reio=0.1;
  
  /** - define which perturbations and sources should be computed */

  ppt->has_scalars=_TRUE_;
  ppt->has_vectors=_FALSE_;
  ppt->has_tensors=_FALSE_;
  ppt->has_ad=_TRUE_;
  ppt->has_bi=_FALSE_;
  ppt->has_cdi=_FALSE_;
  ppt->has_nid=_FALSE_;
  ppt->has_niv=_FALSE_;
  ppt->has_source_t=_TRUE_;
  ppt->has_source_p=_TRUE_;
  ppt->has_source_l=_FALSE_;

  /** - state whether Cl's (spectra in harmonic space) will be required */

  ptr->has_cls=_TRUE_;

  /** - define the primordial spectrum */

  ppm->primordial_spec_type = smooth_Pk;
  ppm->A_s_ad = 1. ;
  ppm->n_s_ad = 1. ;
  ppm->alpha_s_ad = 0. ;  
  ppm->k_pivot = 0.05; /* in Mpc-1 */

  /** - name of output files */

  pop->cls_ad = "output/cls.dat";

  /** - amount of information sent to standard output (none if all set to zero) */

  pba->background_verbose = 1;
  pth->thermodynamics_verbose = 1;
  ppt->perturbations_verbose = 1;
  pbs->bessels_verbose = 1;
  ppm->primordial_verbose = 1;
  ptr->transfer_verbose = 1;
  psp->spectra_verbose = 1;
  pop->output_verbose = 1;

  return _SUCCESS_;

}
