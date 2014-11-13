/** @file custom_lensing.c 
 * Julien Lesgourgues, 21.06.2011    
 */
 
/* purpose: compute lensed Cl's, assuming lensing power spectrum from
   another model.

   --------------------------------------------------------
   usage (if you compute the lensing potential with CLASS):
   -------------------------------------------------------- 

   1. prepare two input files, differing only through comsological
   parameters and root name. For instance, one file 'lensing.ini'
   writh root 'lensing_...' and one file 'primordial.ini' with root
   'primordial_...'. The two files should have at least
   'output=tCl,pCl,lCl', and 'lensing = yes'.

   2. run the model for the lensing potential:

   ./class lensing.ini

   3. below, set the name of the file where the lensing potential will
   be read, e.g. sprintf(filename,"output/lensing_cl.dat");
   
   4. run with

   ./custom_lensing primordial.ini

   in the output file 'primordial_cl.dat', you can see the ClTT from
   one model and the Clphiphi from the other model. In the file
   'primordial_cl_lensed.dat', the Cl's of one model have been lensed
   with the potential of the other model.

   5. to cross-check, you can compare with the lensed spectra that you
   would get if not assuming another lensing potential. For this
   purpose, just change the root name of primordial.ini to
   e.g. 'normal_...', and run with

   /class primordial.ini

   There will be differences between 'primordial_cl.dat' and
   'normal_cl.dat' (only in the phiphi column), and between
   'primordial_cl_lensed.dat' and 'normal_cl_lensed.dat' (all
   columns).

   --------------------------------------------------------
   usage (if you compute the lensing potential yourself):
   -------------------------------------------------------- 

   Instead of computing 'lensing_cl.dat' with CLASS, you need to
   compute the lesning potential spectrum yourself, and to write the
   results in a file with at least two columns: l from 2 to l_max (by
   default, l_max=3250), and [l(l+1)/2pi]C_l^phiphi. Below, you will
   need to adapat the name of the file, the number of lines of
   comments, the number of l's to be read, and the number of columns.

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
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  FILE * output;
  char filename[200];
  char junk_string[200];
  int index_l;
  int l,l_read;
  float cl_read;
  double factor;

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
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

  if (spectra_init(&pr,&ba,&pt,&tr,&pm,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(&pr,&ba,&th,&pm,&sp,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  /* read model with appropriate lensing potential */

  sprintf(filename,"output/lensing_cl.dat");

  output = fopen(filename,"r");
  fprintf(stderr,"Read reference in file %s\n",filename);

  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);
  fgets(junk_string,200,output);

  index_l=0;

  for (l=2; l <= sp.l_max_tot; l++) {
    fscanf(output,"%d",&l_read); // l
    if (l_read != l) {
      printf("l_read != l: %d %d\n",l_read,l);
    }

    fscanf(output,"%e",&cl_read); // TT
    fscanf(output,"%e",&cl_read); // TE
    fscanf(output,"%e",&cl_read); // EE
    fscanf(output,"%e",&cl_read); // BB
    fscanf(output,"%e",&cl_read); // phiphi
    
    if (l == (int)sp.l[index_l]) {

      factor = l*(l+1)/2./_PI_;

      sp.cl[sp.index_md_scalars][index_l*sp.ct_size+sp.index_ct_pp]=cl_read/factor;
      index_l++;
      
    }

    fscanf(output,"%e",&cl_read); // Tphi
    fscanf(output,"%e",&cl_read); // Ephi
  }
  
  fclose(output);

  /* spline the spectra (redundent excepted for the new phiphi column) */

  class_call(array_spline_table_lines(sp.l,
				      sp.l_size[sp.index_md_scalars],
				      sp.cl[sp.index_md_scalars],
				      sp.ic_ic_size[sp.index_md_scalars]*sp.ct_size,
				      sp.ddcl[sp.index_md_scalars],
				      _SPLINE_EST_DERIV_,
				      sp.error_message),
	     sp.error_message,
	     sp.error_message);

  /* carry on the calculation */

  if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (output_init(&ba,&pt,&sp,&nl,&le,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }

  /****** all calculations done, now free the structures ******/

  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

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
