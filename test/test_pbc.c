/** @file class.c 
 * Julien Lesgourgues, 20.04.2010    
 */
 
#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct bessels bs;          /* for bessel functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct lensing le;          /* for lensing spectra */
  struct output op;           /* for output files */
  struct spectra_nl nl;       /* for calculation of non-linear spectra */
  ErrorMsg errmsg;

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&le,&op,&nl,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg); 
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (spectra_init(&ba,&pt,&tr,&pm,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (output_init(&ba,&pt,&sp,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }

  /****** done ******/

  int index_mode=0;
  int index_eta;
  int index_ic=0;
  int index_k;
  double delta_rho_bc,rho_bc;
  double delta_i,rho_i;
  double P_bc;
  int last_index_back;
  double * pvecback_long;
  double k,pk;
  FILE * output;
  double z,eta;
  double * tki;

  if(pt.has_matter_transfers == _FALSE_) {
    printf("You need to switch on mTk calculation for this code\n");
    return _FAILURE_;
  }

  output=fopen("output/Pcb.dat","w");
  
  class_alloc(pvecback_long,sizeof(double)*ba.bg_size,errmsg);
  class_alloc(tki,sizeof(double)*sp.ln_k_size*sp.tr_size,errmsg);
  
  z=0.;
  
  if(background_eta_of_z(&ba,z,&eta) == _FAILURE_) {
    printf("\n\nError running background_eta_of_z \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }
  
  if(background_at_eta(&ba,
		       eta,
		       long_info, 
		       normal, 
		       &last_index_back, 
		       pvecback_long) == _FAILURE_) {
    printf("\n\nError running background_at_eta \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }
  
  if (spectra_tk_at_z(&ba,&sp,z,tki)  == _FAILURE_) {
    printf("\n\nError in spectra_tk_at_z \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }
  
  for (index_k=0; index_k<sp.ln_k_size; index_k++) {
    
    k=exp(sp.ln_k[index_k]);
    
    delta_rho_bc=0.;
    rho_bc=0.;
    
    /* T_b(k,eta) */
    
    delta_i = tki[index_k*sp.tr_size+sp.index_tr_b];
    
    rho_i = pvecback_long[ba.index_bg_rho_b];
    
    delta_rho_bc += rho_i * delta_i;
    
    rho_bc += rho_i;
    
    /* T_cdm(k,eta) */
    
    if (ba.has_cdm == _TRUE_) {
      
      delta_i = tki[index_k*sp.tr_size+sp.index_tr_cdm];
      
      rho_i = pvecback_long[ba.index_bg_rho_cdm];
      
      delta_rho_bc += rho_i * delta_i;
      
      rho_bc += rho_i;
      
    }
      
    if (primordial_spectrum_at_k(&pm,index_mode,linear,k,&pk) == _FAILURE_) {
      printf("\n\nError in primordial_spectrum_at_k \n=>%s\n",pm.error_message);
      return _FAILURE_;
    }
    
    P_bc=pk*pow(delta_rho_bc/rho_bc,2)*2.*_PI_*_PI_/k/k/k;
    
    fprintf(output,"%e   %e\n",k,P_bc);
    
  }

  /******************/

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
