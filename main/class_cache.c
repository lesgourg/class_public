/** @file class.c 
 * Julien Lesgourgues, 17.04.2011    
 */
 
#include "class.h"
#include <sys/shm.h>
#include <sys/stat.h>
#include <errno.h>

int main(int argcin, char **argv) {

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

  ErrorMsg errmsg;
  int initialise_cache = _FALSE_;  
  int index_q, index_q_max=16384, argc=argcin, shmid;
  FILE * q_file;

  /** Check last argument to see if we are in cache init mode or cache free mode*/
  if (strncmp(argv[argc-1],"-i",2)==0){
    printf("Initialise cache...%s \n",argv[argc-1]);
    initialise_cache = _TRUE_;
    argc--;
  }
  else if (strncmp(argv[argc-1],"-f",2)==0){
    for(index_q=0; index_q<index_q_max; index_q++){
      shmid = shmget(_SHARED_MEMORY_KEYS_START_+index_q, 0, S_IRUSR | S_IWUSR);
      if (shmid<0)
        continue;
      shmctl (shmid, IPC_RMID, 0);
    }
    printf("All memory segments free'ed\n");
    return _SUCCESS_;
  }
  
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

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (bessel_init(&pr,&bs) == _FAILURE_) {
    printf("\n\nError in bessel_init \n =>%s\n",bs.error_message);
    return _FAILURE_;
  }

  //Initialise cache during call to transfer_init:
  if (initialise_cache == _TRUE_){
    tr.initialise_HIS_cache=_TRUE_;
    bs.get_HIS_from_shared_memory=_FALSE_;

    if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
      printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
      return _FAILURE_;
    }
    q_file=fopen("q_from_cache.dat","w");
    for (index_q=0; index_q<tr.q_size; index_q++){
      fprintf(q_file,"%.16e ",tr.q[index_q]);
    }
    fclose(q_file);
    q_file = fopen("q_list_binary.bin","wb");
    fwrite ( tr.q, sizeof(double),tr.q_size,q_file );
    fclose(q_file);
  }
  else{
    if (transfer_init(&pr,&ba,&th,&pt,&bs,&tr) == _FAILURE_) {
      printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
      return _FAILURE_;
    }

    if (spectra_init(&pr,&ba,&pt,&tr,&pm,&sp) == _FAILURE_) {
      printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
      return _FAILURE_;
    }

    if (nonlinear_init(&pr,&ba,&th,&pt,&bs,&tr,&pm,&sp,&nl) == _FAILURE_) {
      printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
      return _FAILURE_;
    }

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
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (bessel_free(&bs) == _FAILURE_)  {
    printf("\n\nError in bessel_free \n=>%s\n",bs.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
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
