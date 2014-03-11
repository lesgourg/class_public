/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

/* this main calls CLASS several times in a loop, with different input
   parameters. It illustrates how the code could be interfaced with a
   parameter extraction code. */

#include "class.h"

int class(
          struct file_content *pfc,
          struct precision * ppr,
          struct background * pba,
          struct thermo * pth,
          struct perturbs * ppt,
          struct primordial * ppm,
          struct nonlinear * pnl,
          struct transfers * ptr,
          struct spectra * psp,
          struct lensing * ple,
          struct output * pop,
          int l_max,
          double ** cl,
          ErrorMsg errmsg) {

  int l;

  class_call(input_init(pfc,ppr,pba,pth,ppt,ptr,ppm,psp,pnl,ple,pop,errmsg),
             errmsg,
             errmsg);

  class_call(background_init(ppr,pba),
             pba->error_message,
             errmsg);

  class_call(thermodynamics_init(ppr,pba,pth),
             pth->error_message,
             errmsg);

  class_call(perturb_init(ppr,pba,pth,ppt),
             ppt->error_message,
             errmsg);

  class_call(primordial_init(ppr,ppt,ppm),
             ppm->error_message,
             errmsg);

  class_call(nonlinear_init(ppr,pba,pth,ppt,ppm,pnl),
             pnl->error_message,
             errmsg);

  class_call(transfer_init(ppr,pba,pth,ppt,pnl,ptr),
             ptr->error_message,
             errmsg);

  class_call(spectra_init(ppr,pba,ppt,ppm,pnl,ptr,psp),
             psp->error_message,
             errmsg);

  class_call(lensing_init(ppr,ppt,psp,pnl,ple),
             ple->error_message,
             errmsg);

  /****** write the Cl values in the input array cl[l]  *******/

  for (l=2; l <= l_max; l++) {

    class_call(output_total_cl_at_l(psp,ple,pop,(double)l,cl[l]),
               psp->error_message,
               errmsg);
  }

  /****** all calculations done, now free the structures ******/

  class_call(lensing_free(ple),
             ple->error_message,
             errmsg);

  class_call(spectra_free(psp),
             psp->error_message,
             errmsg);

  class_call(transfer_free(ptr),
             ptr->error_message,
             errmsg);

  class_call(nonlinear_free(pnl),
             pnl->error_message,
             errmsg);

  class_call(primordial_free(ppm),
             ppm->error_message,
             errmsg);

  class_call(perturb_free(ppt),
             ppt->error_message,
             errmsg);

  class_call(thermodynamics_free(pth),
             pth->error_message,
             errmsg);

  class_call(background_free(pba),
             pba->error_message,
             errmsg);

  return _SUCCESS_;

}

int main() {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  int i;
  int l,l_max;
  int num_ct_max=7;
  double ** cl;
  struct file_content fc;

  /* choose a value of l_max in C_l's */
  l_max=3000;

  /* all parameters for which we don't want to keep default values
     should be passed to the code through a file_content
     structure. Create such a structure with the size you need: 9 in
     this exemple */
  parser_init(&fc,10,"",errmsg);

  /* assign values to these 9 parameters. Some will be fixed, some
     will be varied in the loop. */
  strcpy(fc.name[0],"output");
  strcpy(fc.value[0],"tCl,pCl,lCl");

  strcpy(fc.name[1],"l_max_scalars");
  sprintf(fc.value[1],"%d",l_max);

  strcpy(fc.name[2],"lensing");
  sprintf(fc.value[2],"yes");

  strcpy(fc.name[3],"H0");
  sprintf(fc.value[3],"%e",72.);

  strcpy(fc.name[4],"omega_b");
  sprintf(fc.value[4],"%e",0.024);

  strcpy(fc.name[5],"omega_cdm");
  sprintf(fc.value[5],"%e",0.05);

  strcpy(fc.name[6],"z_reio");
  sprintf(fc.value[6],"%e",10.);

  strcpy(fc.name[7],"A_s");
  sprintf(fc.value[7],"%e",2.3e-9);

  strcpy(fc.name[8],"n_s");
  sprintf(fc.value[8],"%e",1.);

  strcpy(fc.name[9],"T_cmb");
  sprintf(fc.value[9],"%e",2.726);

  /* allocate the array where calculated Cl's will be written (we
     could add another array with P(k), or extract other results from
     the code - here we assume that we are interested in the C_l's
     only */
  cl=malloc((l_max+1)*sizeof(double*));
  for (l=0;l<=l_max;l++)
    cl[l]=malloc(num_ct_max*sizeof(double));

  /* now, loop over values of some of these parameters: in this
     exemple, omega_b */
  for (i=0; i<=10; i++) {

    /* assign one value to omega_b */
    sprintf(fc.value[4],"%e",0.01+i*0.002);
    printf("#run with omega_b = %s\n",fc.value[4]);

    /* calls class and return the C_l's*/
    if (class(&fc,&pr,&ba,&th,&pt,&pm,&nl,&tr,&sp,&le,&op,l_max,cl,errmsg) == _FAILURE_) {
      printf("\n\nError in class \n=>%s\n",errmsg);
      return _FAILURE_;
    }

    /* print the lensed C_l^TT, C_l^EE, C_l^TE's for this run */
    for (l=2;l<=l_max;l++) {
      fprintf(stdout,"%d  %e  %e  %e\n",
	      l,
	      cl[l][sp.index_ct_tt],
	      cl[l][sp.index_ct_ee],
	      cl[l][sp.index_ct_te]);
    }
    fprintf(stdout,"\n");
  }

  for (l=0;l<l_max;l++)
    free(cl[l]);
  free(cl);

  return _SUCCESS_;

}
