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
  double Neff;

  /** - assign values to background cosmological parameters */

  /* H0 in Mpc^{-1} = h / 2999.7 */
  h=0.7;
  pba->H0 = h / 2999.7;

  /* effective number of neutrinos (following the usual definition) */
  Neff=3.04;

  /* Omega's */
  pba->Omega0_g = 2.3812e-5/h/h*pow(2.726/2.7,4.);
  pba->Omega0_nur = Neff*7./8.*pow(4./11.,4./3.)*pba->Omega0_g;
  pba->Omega0_b = 0.05;
  pba->Omega0_cdm = 0.23;
  pba->Omega0_lambda = 1.-pba->Omega0_b-pba->Omega0_cdm-pba->Omega0_g-pba->Omega0_nur;
  pba->Omega0_de = 0.;

  /** - assign values to thermodynamics cosmological parameters */

  pth->Tcmb=2.726;
  pth->YHe=0.25;
  pth->reio_parametrization=reio_camb; /* reio_camb means "same form for X_e(z) as in camb" */
  pth->reio_z_or_tau=reio_tau;
  pth->z_reio=10.;   /* used only if above set to reio_z */
  pth->tau_reio=0.1; /* used only if above set to reio_tau */
  
  /** - define which perturbations and sources should be computed */

  ppt->has_scalars=_TRUE_;  
  ppt->has_vectors=_FALSE_;
  ppt->has_tensors=_FALSE_;

  ppt->has_ad=_TRUE_;
  ppt->has_bi=_FALSE_;
  ppt->has_cdi=_FALSE_;
  ppt->has_nid=_FALSE_;
  ppt->has_niv=_FALSE_;

  ppt->has_cl_cmb_temperature = _TRUE_;
  ppt->has_cl_cmb_polarization = _TRUE_;
  ppt->has_cl_cmb_lensing_potential = _FALSE_;
  ppt->has_pk_matter = _FALSE_;

  /** - define the primordial spectrum */

  ppm->primordial_spec_type = smooth_Pk;
  ppm->A_s_ad = 1. ; /* amplitude */
  ppm->n_s_ad = 1. ; /* tilt */
  ppm->alpha_s_ad = 0. ; /* running */ 
  ppm->k_pivot = 0.05; /* pivot wavenumber in Mpc-1 */

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
