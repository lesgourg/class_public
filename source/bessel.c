/** @file bessel.c Documented Bessel module.
 *
 * Julien Lesgourgues, 26.08.2010
 *
 * This module loads spherical Bessel functions
 * (either read from file or computed from scratch).
 *
 * Hence the following functions can be called from other modules:
 *
 * -# bessel_init() at the beginning (anytime after input_init() and before transfer_init())
 * -# bessel_at_x() at any time for computing a value j_l(x) at any x by interpolation
 * -# bessel_free() at the end
 */

#include "bessel.h"

/**
 * Bessel function for arbitrary argument x.
 *
 * Evaluates the spherical Bessel function x at a given value of x by
 * interpolating in the pre-computed table.  This function can be
 * called from whatever module at whatever time, provided that
 * bessel_init() has been called before, and bessel_free() has not
 * been called yet.
 *
 * @param pbs     Input: pointer to bessels structure
 * @param x       Input: argument x
 * @param index_l Input: index defining l = pbs->l[index_l]
 * @param j       Ouput: \f$j_l(x) \f$
 * @return the error status
 */
int bessel_at_x(
		struct bessels * pbs,
		double x,
		int index_l,
		double * j
		) {


  /** Summary: */

  /** - define local variables */

  int index_x;          /* index in the interpolation table */
  double a;         /* quantities for the splint interpolation formula */

  /** - if index_l is too large to be in the interpolation table, return  an error */

  class_test(index_l > pbs->l_size,
	     pbs->error_message,
	     "index_l=%d>l_size=%d; increase l_max.",index_l,pbs->l_size);

  /** - if x is too small to be in the interpolation table, return 0 */

  if (x < *(pbs->x_min[index_l])) {
    *j=0;
    return _SUCCESS_;
  }
  else {

    /** - if x is too large to be in the interpolation table, return an error (this should never occur since x_max in the table should really be the highest value needed by the code, given the precision parameters) */

    class_test(x > pbs->x_max,
	       pbs->error_message,
	       "x=%e>x_max=%e in bessel structure",x,pbs->x_max);

    /** - otherwise, interpolation is needed: */

    /** (a) find index_x, i.e. the position of x in the table; no complicated algorithm needed, since values are regularly spaced with a known step and known first value */

    index_x = (int)((x-*(pbs->x_min[index_l]))/pbs->x_step);

    /** (b) find result with the splint algorithm (equivalent to the one in numerical recipies, although terms are rearranged differently to minimize number of operations) */
    a = (*(pbs->x_min[index_l])+pbs->x_step*(index_x+1) - x)/pbs->x_step;
    *j= a * pbs->j[index_l][index_x]
      + (1.-a) * ( pbs->j[index_l][index_x+1]
		   - a * ((a+1.) * pbs->ddj[index_l][index_x]
			  +(2.-a) * pbs->ddj[index_l][index_x+1])
		   * pbs->x_step * pbs->x_step / 6.0);
  }

  return _SUCCESS_;

}

/**
 * Get spherical Bessel functions (either read from file or compute
 * from scratch).
 *
 * Each table of spherical Bessel functions \f$ j_l(x) \f$ corresponds
 * to a set of values for:
 *
 * -# pbs->l[index_l]: list of l values l of size pbs->l_size
 * -# pbs->x_step: step dx for sampling Bessel functions \f$ j_l(x) \f$
 * -# pbs->x_max: last value of x (always a multiple of x_step!)
 * -# pbs->j_cut: value of \f$ j_l \f$ below which it is approximated by zero (in the region x << l)
 *
 * This function checks whether there is alread a file "bessels.dat"
 * with the same l's, x_step, x_max, j_cut).
 * If yes, it fills the table of bessel functions (and
 * their second derivatives, needed for spline interpolation) with the
 * values read from the file. If not, it computes all values using
 * bessel_j_for_l(), and stores them both in the bessels
 * stucture pbs, and in a file "bessels.dat" (in view of the next
 * runs).
 *
 * @param ppr Input : pointer to precision strucutre
 * @param pbs Output: initialized bessel structure
 * @return the error status
 */

int bessel_init(
		struct precision * ppr,
		struct bessels * pbs
		) {

  /** Summary: */

  /** - define local variables */

  /* index for l (since first value of l is always 2, l=index_l+2) */
  int index_l;

  /* first numbers to be read in bessels.dat file */
  int l_size_file;
  int * l_file;
  double x_step_file;
  double x_max_file;
  double j_cut_file;
  int has_dj_file;

  int num_j;

  /* bessels.dat file */
  FILE * bessel_file;

  /* This code can be optionally compiled with the openmp option for parallel computation.
     Inside parallel regions, the use of the command "return" is forbidden.
     For error management, instead of "return _FAILURE_", we will set the variable below
     to "abort = _TRUE_". This will lead to a "return _FAILURE_" jus after leaving the
     parallel region. */
  int abort;

#ifdef _OPENMP
  /* instrumentation times */
  double tstart, tstop;
#endif

  if (pbs->use_pbs == _FALSE_) {
    if (pbs->bessels_verbose > 0)
      printf("Bessel functions will be computed on the fly by hyperspherical module. Bessel module skipped.\n");
    return _SUCCESS_;
  }

  if (pbs->l_max == 0) {
    if (pbs->bessels_verbose > 0)
      printf("No harmonic space transfer functions to compute. Bessel module skipped.\n");
    return _SUCCESS_;
  }

  /** - infer l values from precision parameters and from l_max */

  class_call(bessel_get_l_list(ppr,pbs),
	     pbs->error_message,
	     pbs->error_message);

  /** - check x_step, x_max and j_cut from precision parameters and from l_max */

  class_test(pbs->x_step <= 0.,
	     pbs->error_message,
	     "x_step=%e, stop to avoid segmentation fault",pbs->x_step);

  class_test(pbs->x_max <= 0.,
	     pbs->error_message,
	     "x_max=%e, stop to avoid segmentation fault",pbs->x_max);

  pbs->j_cut = ppr->bessel_j_cut;

  /** - do we need to store also j_l'(x) ? */
  // added for new version
  pbs->has_dj = _TRUE_;

  /** - check if file bessels.dat already exists with the same (l's, x_step, x_max, j_cut). If yes, read it. */

  if (pbs->bessel_always_recompute == _FALSE_) {

    bessel_file=fopen(ppr->bessel_file_name,"r");

    if (bessel_file == NULL) {
      if (pbs->bessels_verbose > 1)
	printf("File %s did not exist.\n",ppr->bessel_file_name);
    }
    else {

      class_test(fread(&l_size_file,sizeof(int),1,bessel_file) != 1,
		 pbs->error_message,
		 "Could not read in bessel file");

      class_alloc(l_file,l_size_file * sizeof(int),pbs->error_message);

      for (index_l=0; index_l < l_size_file; index_l++) {
	class_test(fread(&l_file[index_l],sizeof(int),1,bessel_file) != 1,
		   pbs->error_message,
		   "Could not read in bessel file");
      }

      class_test(fread(&x_step_file,sizeof(double),1,bessel_file) != 1,
		 pbs->error_message,
		 "Could not read in bessel file");

      class_test(fread(&x_max_file,sizeof(double),1,bessel_file) != 1,
		 pbs->error_message,
		 "Could not read in bessel file");

      class_test(fread(&j_cut_file,sizeof(double),1,bessel_file) != 1,
		 pbs->error_message,
		 "Could not read in bessel file");

      class_test(fread(&has_dj_file,sizeof(int),1,bessel_file) != 1,
		 pbs->error_message,
		 "Could not read in bessel file");

      index_l=0;

      if (l_size_file == pbs->l_size) {
	while ((pbs->l[index_l] == l_file[index_l]) && (index_l < pbs->l_size-1)) {
	  index_l++;
	}
	if (pbs->l[pbs->l_size-1] == l_file[pbs->l_size-1])
	  index_l++;
      }

      free(l_file);

      if ((index_l == pbs->l_size) &&
	  (x_step_file == pbs->x_step) &&
	  (j_cut_file == pbs->j_cut) &&
	  (x_max_file == pbs->x_max) &&
	  (has_dj_file == pbs->has_dj)) {

	if (pbs->bessels_verbose > 0)
	  printf("Read bessels in file %s\n",ppr->bessel_file_name);

	class_alloc(pbs->x_size,pbs->l_size*sizeof(int*),pbs->error_message);
	class_alloc(pbs->x_min,pbs->l_size*sizeof(double*),pbs->error_message);
	class_alloc(pbs->buffer,pbs->l_size*sizeof(double*),pbs->error_message);
	class_alloc(pbs->j,pbs->l_size*sizeof(double*),pbs->error_message);
	class_alloc(pbs->ddj,pbs->l_size*sizeof(double*),pbs->error_message);
	if (pbs->has_dj == _TRUE_) {
	  class_alloc(pbs->dj,pbs->l_size*sizeof(double*),pbs->error_message);
	  class_alloc(pbs->dddj,pbs->l_size*sizeof(double*),pbs->error_message);
	}

	class_test(fread(pbs->x_size,sizeof(int),pbs->l_size,bessel_file) != pbs->l_size,
		   pbs->error_message,
		   "Could not read in bessel file");

	pbs->x_size_max=0;

	for (index_l=0; index_l < pbs->l_size; index_l++) {

	  if (pbs->x_size[index_l] > pbs->x_size_max)
	    pbs->x_size_max=pbs->x_size[index_l];

          if (pbs->has_dj == _TRUE_) {
	    num_j = 4;
	  }
	  else {
	    num_j = 2;
	  }

	  class_alloc(pbs->buffer[index_l],
		      (1+num_j*pbs->x_size[index_l])*sizeof(double),
		      pbs->error_message);

	  pbs->x_min[index_l] = pbs->buffer[index_l];
	  pbs->j[index_l] = pbs->buffer[index_l]+1;
	  pbs->ddj[index_l] = pbs->j[index_l]+pbs->x_size[index_l];
	  if (pbs->has_dj == _TRUE_) {
	    pbs->dj[index_l] = pbs->ddj[index_l]+pbs->x_size[index_l];
	    pbs->dddj[index_l] = pbs->dj[index_l]+pbs->x_size[index_l];
	  }

	  class_test(fread(pbs->x_min[index_l],sizeof(double),1,bessel_file) != 1,
		     pbs->error_message,
		     "Could not read in bessel file");

	  class_test(fread(pbs->j[index_l],sizeof(double),pbs->x_size[index_l],bessel_file) != pbs->x_size[index_l],
		     pbs->error_message,
		     "Could not read in bessel file");

	  class_test(fread(pbs->ddj[index_l],sizeof(double),pbs->x_size[index_l],bessel_file) != pbs->x_size[index_l],
		     pbs->error_message,
		     "Could not read in bessel file");

	  if (pbs->has_dj == _TRUE_) {

	    class_test(fread(pbs->dj[index_l],sizeof(double),pbs->x_size[index_l],bessel_file) != pbs->x_size[index_l],
		       pbs->error_message,
		       "Could not read in bessel file");

	    class_test(fread(pbs->dddj[index_l],sizeof(double),pbs->x_size[index_l],bessel_file) != pbs->x_size[index_l],
		       pbs->error_message,
		       "Could not read in bessel file");

	  }
	}

	fclose(bessel_file);

	return _SUCCESS_;

      }
      else {
	fclose(bessel_file);
      }
    }
  }

  /** - if not, compute form scratch : */

  if (pbs->bessels_verbose > 0)
    printf("Computing bessels\n");

  class_alloc(pbs->x_size,pbs->l_size*sizeof(int*),pbs->error_message);
  class_alloc(pbs->buffer,pbs->l_size*sizeof(double*),pbs->error_message);
  class_alloc(pbs->x_min,pbs->l_size*sizeof(double*),pbs->error_message);
  class_alloc(pbs->j,pbs->l_size*sizeof(double*),pbs->error_message);
  class_alloc(pbs->ddj,pbs->l_size*sizeof(double*),pbs->error_message);
  if (pbs->has_dj == _TRUE_) {
    class_alloc(pbs->dj,pbs->l_size*sizeof(double*),pbs->error_message);
    class_alloc(pbs->dddj,pbs->l_size*sizeof(double*),pbs->error_message);
  }

  /* initialize error management flag */
  abort = _FALSE_;

  /* beginning of parallel region */

#pragma omp parallel				\
  shared(ppr,pbs,abort)				\
  private(index_l,tstart,tstop)

  {

#ifdef _OPENMP
    tstart = omp_get_wtime();
#endif

#pragma omp for schedule (dynamic)

    /** (a) loop over l and x values, compute \f$ j_l(x) \f$ for each of them */
    for (index_l = 0; index_l < pbs->l_size; index_l++) {

      class_call_parallel(bessel_j_for_l(ppr,pbs,index_l),
			  pbs->error_message,
			  pbs->error_message);

#pragma omp flush(abort)

    } /* end of loop over l */

#ifdef _OPENMP
    tstop = omp_get_wtime();
    if (pbs->bessels_verbose > 1)
      printf("In %s: time spent in parallel region (loop over l's) = %e s for thread %d\n",
	     __func__,tstop-tstart,omp_get_thread_num());
#endif

  } /* end of parallel region */

  if (abort == _TRUE_) return _FAILURE_;

  pbs->x_size_max=0;
  for (index_l=0; index_l < pbs->l_size; index_l++)
    if (pbs->x_size[index_l] > pbs->x_size_max)
      pbs->x_size_max=pbs->x_size[index_l];

  if (pbs->bessel_always_recompute == _FALSE_) {

    if (pbs->bessels_verbose > 0)
      printf(" -> (over)write in file %s\n",ppr->bessel_file_name);

    /** (b) write in file */

    bessel_file = fopen(ppr->bessel_file_name,"w");

    fwrite(&(pbs->l_size),sizeof(int),1,bessel_file);
    fwrite(pbs->l,sizeof(int),pbs->l_size,bessel_file);
    fwrite(&(pbs->x_step),sizeof(double),1,bessel_file);
    fwrite(&(pbs->x_max),sizeof(double),1,bessel_file);
    fwrite(&(pbs->j_cut),sizeof(double),1,bessel_file);
    fwrite(&(pbs->has_dj),sizeof(int),1,bessel_file);

    fwrite(pbs->x_size,sizeof(int),pbs->l_size,bessel_file);

    for (index_l=0; index_l<pbs->l_size; index_l++) {
      fwrite(pbs->x_min[index_l],sizeof(double),1,bessel_file);
      fwrite(pbs->j[index_l],sizeof(double),pbs->x_size[index_l],bessel_file);
      fwrite(pbs->ddj[index_l],sizeof(double),pbs->x_size[index_l],bessel_file);
      if (pbs->has_dj == _TRUE_) {
	fwrite(pbs->dj[index_l],sizeof(double),pbs->x_size[index_l],bessel_file);
	fwrite(pbs->dddj[index_l],sizeof(double),pbs->x_size[index_l],bessel_file);
      }
    }

    fclose(bessel_file);

  }

  return _SUCCESS_;
}

/**
 * Free all memory space allocated by bessel_init().
 *
 * To be called at the end of each run.
 *
 * @param pbs Input : Initialized bessel structure
 * @return the error status
 */

int bessel_free(
		struct bessels * pbs) {

  int index_l;

  if ((pbs->l_max > 0) && (pbs->use_pbs == _TRUE_)) {

    for (index_l = 0; index_l < pbs->l_size; index_l++) {
      free(pbs->buffer[index_l]);
    }
    free(pbs->buffer);
    free(pbs->x_min);
    free(pbs->j);
    free(pbs->ddj);
    if (pbs->has_dj == _TRUE_) {
      free(pbs->dj);
      free(pbs->dddj);
    }
    free(pbs->x_size);
    free(pbs->l);

  }

  return _SUCCESS_;
}

/**
 * Define number and values of mutipoles l. This is crucial since not
 * only the Bessel functions, but also the transfer functions and
 * anisotropy spectra C_l will automatically be sampled at the same
 * values (there would be no logic in having various l lists differing
 * from each other).
 *
 *
 * @param ppr Input : pointer to precision structure
 * @param pbs Input/Output : pointer to besseld structure (result stored here)
 * @return the error status
 */

int bessel_get_l_list(
		      struct precision * ppr,
		      struct bessels * pbs
		      ) {

  /** Summary: */

  /** - define local variables */

  int index_l,increment,current_l;

  /** - start from l = 2 and increase with logarithmic step */

  index_l = 0;
  current_l = 2;
  increment = MAX((int)(current_l * (ppr->l_logstep-1.)),1);

  while (((current_l+increment) < pbs->l_max) &&
	 (increment < ppr->l_linstep)) {

    index_l ++;
    current_l += increment;
    increment = MAX((int)(current_l * (ppr->l_logstep-1.)),1);

  }

  /** - when the logarithmic step becomes larger than some linear step,
      stick to this linear step till l_max */

  increment = ppr->l_linstep;

  while ((current_l+increment) <= pbs->l_max) {

    index_l ++;
    current_l += increment;

  }

  /** - last value set to exactly l_max */

  if (current_l != pbs->l_max) {

    index_l ++;
    current_l = pbs->l_max;

  }

  pbs->l_size = index_l+1;

  /** - so far we just counted the number of values. Now repeat the
      whole thing but fill array with values. */

  class_alloc(pbs->l,pbs->l_size*sizeof(int),pbs->error_message);

  index_l = 0;
  pbs->l[0] = 2;
  increment = MAX((int)(pbs->l[0] * (ppr->l_logstep-1.)),1);

  while (((pbs->l[index_l]+increment) < pbs->l_max) &&
	 (increment < ppr->l_linstep)) {

    index_l ++;
    pbs->l[index_l]=pbs->l[index_l-1]+increment;
    increment = MAX((int)(pbs->l[index_l] * (ppr->l_logstep-1.)),1);

  }

  increment = ppr->l_linstep;

  while ((pbs->l[index_l]+increment) <= pbs->l_max) {

    index_l ++;
    pbs->l[index_l]=pbs->l[index_l-1]+increment;

  }

  if (pbs->l[index_l] != pbs->l_max) {

    index_l ++;
    pbs->l[index_l]= pbs->l_max;

  }

  return _SUCCESS_;

}

/**
 * Get spherical Bessel functions for given value of l.
 *
 * Find the first value x_min(l) at which the function is not
 * negligible (for large l values, Bessel functions are very close to
 * zero nearly until x=l).
 * Then, sample it with step x_step till x_max.
 *
 * @param ppr Input : pointer to precision structure
 * @param pbs Input/Output : pointer to bessel structure (store result here)
 * @return the error status
 */

int bessel_j_for_l(
		   struct precision * ppr,
		   struct bessels * pbs,
		   int index_l
		   ){

  /** Summary: */

  /** - define local variables */

  /* index for x and value x=x_min[index_l]+x_step*index_x */
  int index_x;
  double x;

  /* value of j_l(x), j_{l-1}(x) returned by bessel_j(); plus j_l'(x) */
  double j,jm,jprime;

  /* for computing x_min */
  double x_min_up;
  double x_min_down;
  double x_min;

  int num_j;

  index_x=0;
  j = 0.;

  /** - find x_min[index_l] by bisection */

  x_min_up=(double)pbs->l[index_l]+0.5;
  x_min_down=0.;

  class_call(bessel_j(pbs,
		      pbs->l[index_l], /* l */
		      x_min_up, /* x */
		      &j),  /* j_l(x) */
	     pbs->error_message,
	     pbs->error_message);

  class_test(j < pbs->j_cut,
	     pbs->error_message,
	     "in bisection, wrong initial guess for x_min_up.");

  while ((x_min_up-x_min_down)/x_min_down > ppr->bessel_tol_x_min) {

    class_test((x_min_up-x_min_down) < ppr->smallest_allowed_variation,
	       pbs->error_message,
	       "(x_min_up-x_min_down) =%e < machine precision : maybe kmin=%e is too small",
	       (x_min_up-x_min_down),ppr->bessel_tol_x_min);

    class_call(bessel_j(pbs,
			pbs->l[index_l], /* l */
			0.5 * (x_min_up+x_min_down), /* x */
			&j),  /* j_l(x) */
	       pbs->error_message,
	       pbs->error_message);

    if (j >= pbs->j_cut)
      x_min_up=0.5 * (x_min_up+x_min_down);
    else
      x_min_down=0.5 * (x_min_up+x_min_down);

  }

  x_min = x_min_down;

  //DEBUG/HACK! Force xmin = 1e-5;
  x_min = 1e-5;

  class_call(bessel_j(pbs,
		      pbs->l[index_l], /* l */
		      x_min, /* x */
		      &j),  /* j_l(x) */
	     pbs->error_message,
	     pbs->error_message);

  /** - define number of x values to be stored (one if all values of j_l(x) were negligible for this l) */

  if (x_min >= pbs->x_max)
    pbs->x_size[index_l] = 1;
  else
    pbs->x_size[index_l] = (int)((pbs->x_max-x_min)/pbs->x_step) + 1;

  /** - allocate memory for x_min[index_l], j[index_l], ddj[index_l] in such way that they stand in a contiguous memory location */

  if (pbs->has_dj == _TRUE_) {
    num_j = 4;
  }
  else {
    num_j = 2;
  }

  class_alloc(pbs->buffer[index_l],
	      (1+num_j*pbs->x_size[index_l])*sizeof(double),
	      pbs->error_message);

  pbs->x_min[index_l] = pbs->buffer[index_l];
  pbs->j[index_l] = pbs->buffer[index_l]+1;
  pbs->ddj[index_l] = pbs->j[index_l] + pbs->x_size[index_l];
  if (pbs->has_dj == _TRUE_) {
    pbs->dj[index_l] = pbs->ddj[index_l] + pbs->x_size[index_l];
    pbs->dddj[index_l] = pbs->dj[index_l] + pbs->x_size[index_l];
  }

  /** - case when all values of j_l(x) were negligible for this l*/

  if (x_min >= pbs->x_max) {

    *(pbs->x_min[index_l]) = pbs->x_max;
    pbs->j[index_l][0]=0.;
    pbs->ddj[index_l][0]=0.;
    if (pbs->has_dj == _TRUE_) {
      pbs->dj[index_l][0]=0.;
      pbs->dddj[index_l][0]=0.;
    }
  }

  /** -otherwise, write first non-negligible value and then loop over x */
  else {

    *(pbs->x_min[index_l]) = x_min;

    pbs->j[index_l][0] = j;

    class_call(bessel_j(pbs,
			pbs->l[index_l]-1, /* l-1 */
			x_min, /* x */
			&jm),  /* j_{l-1}(x) */
	       pbs->error_message,
	       pbs->error_message);

    jprime = jm - (pbs->l[index_l]+1)*j/x_min; /* j_l'=j_{l-1}-(l+1)j_l/x */

    pbs->ddj[index_l][0] = - 2./x_min*jprime
      + (pbs->l[index_l]*(pbs->l[index_l]+1)/x_min/x_min-1.)*j; /* j_l'' = -2/x j_l' + (l(l+1)/x/x-1)*j */

    if (pbs->has_dj == _TRUE_) {

      pbs->dj[index_l][0] = jprime;

      pbs->dddj[index_l][0] = - 2./x_min*pbs->ddj[index_l][0]
	+ ((pbs->l[index_l]*(pbs->l[index_l]+1)+2)/x_min/x_min-1.)*jprime
	- 2.*pbs->l[index_l]*(pbs->l[index_l]+1)/x_min/x_min/x_min*j;
    }

    /* loop over other non-negligible values */
    for (index_x=1; index_x < pbs->x_size[index_l]; index_x++) {

      x = *(pbs->x_min[index_l])+index_x*pbs->x_step;

      class_call(bessel_j(pbs,
			  pbs->l[index_l], /* l */
			  x, /* x */
			  &j),  /* j_l(x) */
		 pbs->error_message,
		 pbs->error_message);

      class_call(bessel_j(pbs,
			  pbs->l[index_l]-1, /* l-1 */
			  x, /* x */
			  &jm),  /* j_{l-1}(x) */
		 pbs->error_message,
		 pbs->error_message);

      jprime = jm - (pbs->l[index_l]+1)*j/x; /* j_l'=j_{l-1}-(l+1)j_l/x */

      pbs->j[index_l][index_x] = j;

      pbs->ddj[index_l][index_x] = - 2./x*jprime
	+ (pbs->l[index_l]*(pbs->l[index_l]+1)/x/x-1.)*j; /* j_l'' = -2/x j_l' + (l(l+1)/x/x-1)*j */

      if (pbs->has_dj == _TRUE_) {

	pbs->dj[index_l][index_x] = jprime;

	pbs->dddj[index_l][index_x] = - 2./x*pbs->ddj[index_l][index_x]
	  + ((pbs->l[index_l]*(pbs->l[index_l]+1)+2)/x/x-1.)*jprime
	  - 2.*pbs->l[index_l]*(pbs->l[index_l]+1)/x/x/x*j;
      }
    }
  }

  return _SUCCESS_;
}


/**
 * Compute spherical Bessel function \f$ j_l(x) \f$ for a given l and x.
 *
 * Inspired from Numerical Recipies.
 *
 * @param pbs Input : pointer to bessel structure (used only for error mess ge)
 * @param l   Input: l value
 * @param x   Input: x value
 * @param jl  Output: \f$ j_l(x) \f$ value
 * @return the error status
 */

int bessel_j(
	     struct bessels * pbs,
	     int l,
	     double x,
	     double * jl
	     ) {

  double nu,nu2,beta,beta2;
  double x2,sx,sx2,cx;
  double cotb,cot3b,cot6b,secb,sec2b;
  double trigarg,expterm,fl;
  double l3,cosb;

  class_test(l < 0,
	     pbs->error_message,
	     " ");

  class_test(x < 0,
	     pbs->error_message,
	     " ");

  fl = (double)l;

  x2 = x*x;

  /************* Use closed form for l<7 **********/

  if (l < 7) {

    sx=sin(x);
    cx=cos(x);

    if(l == 0) {
      if (x > 0.1) *jl=sx/x;
      else *jl=1.-x2/6.*(1.-x2/20.);
      return _SUCCESS_;
    }

    if(l == 1) {
      if (x > 0.2) *jl=(sx/x -cx)/x;
      else *jl=x/3.*(1.-x2/10.*(1.-x2/28.));
      return _SUCCESS_;
    }

    if (l == 2) {
      if (x > 0.3) *jl=(-3.*cx/x-sx*(1.-3./x2))/x;
      else *jl=x2/15.*(1.-x2/14.*(1.-x2/36.));
      return _SUCCESS_;
    }

    if (l == 3) {
      if (x > 0.4) *jl=(cx*(1.-15./x2)-sx*(6.-15./x2)/x)/x;
      else *jl=x*x2/105.*(1.-x2/18.*(1.-x2/44.));
      return _SUCCESS_;
    }

    if (l == 4) {
      if (x > 0.6) *jl=(sx*(1.-45./x2+105./x2/x2) +cx*(10.-105./x2)/x)/x;
      else *jl=x2*x2/945.*(1.-x2/22.*(1.-x2/52.));
      return _SUCCESS_;
    }

    if (l == 5) {
      if (x > 1) *jl=(sx*(15.-420./x2+945./x2/x2)/x -cx*(1.0-105./x2+945./x2/x2))/x;
      else *jl=x2*x2*x/10395.*(1.-x2/26.*(1.-x2/60.));
      return _SUCCESS_;
    }

    if (l == 6) {
      if (x > 1) *jl=(sx*(-1.+(210.-(4725.-10395./x2)/x2)/x2)+
		      cx*(-21.+(1260.-10395./x2)/x2)/x)/x;
      else *jl=x2*x2*x2/135135.*(1.-x2/30.*(1.-x2/68.));
      return _SUCCESS_;
    }

  }

  else {

    if (x <= 1.e-40) {
      *jl=0.0;
      return _SUCCESS_;
    }

    nu= fl + 0.5;
    nu2=nu*nu;

    if ((x2/fl) < 0.5) {
      *jl=exp(fl*log(x/nu/2.)+nu*(1-log(2.))-(1.-(1.-3.5/nu2)/nu2/30.)/12./nu)
	/nu*(1.-x2/(4.*nu+4.)*(1.-x2/(8.*nu+16.)*(1.-x2/(12.*nu+36.))));
      return _SUCCESS_;
    }

    if ((fl*fl/x) < 0.5) {

      beta = x - _PI_/2.*(fl+1.);
      *jl = (cos(beta)*(1.-(nu2-0.25)*(nu2-2.25)/8./x2*(1.-(nu2-6.25)*(nu2-12.25)/48./x2))
	     -sin(beta)*(nu2-0.25)/2./x* (1.-(nu2-2.25)*(nu2-6.25)/24./x2*(1.-(nu2-12.25)*(nu2-20.25)/80./x2)) )/x;

      return _SUCCESS_;

    }

    l3 = pow(nu,0.325);

    if (x < nu-1.31*l3) {

      cosb=nu/x;
      sx=sqrt(nu2-x2);
      cotb=nu/sx;
      secb=x/nu;
      beta=log(cosb+sx/x);
      cot3b=cotb*cotb*cotb;
      cot6b=cot3b*cot3b;
      sec2b=secb*secb;
      expterm=((2.+3.*sec2b)*cot3b/24.
	       - ((4.+sec2b)*sec2b*cot6b/16.
		  + ((16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b/5760.
		     + (32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b/128./nu)*cot6b/nu)/nu)/nu;
      *jl=sqrt(cotb*cosb)/(2.*nu)*exp(-nu*beta+nu/cotb-expterm);

      return _SUCCESS_;

    }

    if (x > nu+1.48*l3) {

      cosb=nu/x;
      sx=sqrt(x2-nu2);
      cotb=nu/sx;
      secb=x/nu;
      beta=acos(cosb);
      cot3b=cotb*cotb*cotb;
      cot6b=cot3b*cot3b;
      sec2b=secb*secb;
      trigarg=nu/cotb-nu*beta-_PI_/4.
	-((2.+3.*sec2b)*cot3b/24.
	  +(16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b*cot6b/5760./nu2)/nu;
      expterm=((4.+sec2b)*sec2b*cot6b/16.
	       -(32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b*cot6b/128./nu2)/nu2;
      *jl=sqrt(cotb*cosb)/nu*exp(-expterm)*cos(trigarg);

      return _SUCCESS_;
    }

    /* last possible case */

    beta=x-nu;
    beta2=beta*beta;
    sx=6./x;
    sx2=sx*sx;
    secb=pow(sx,1./3.);
    sec2b=secb*secb;
    *jl=(_GAMMA1_*secb + beta*_GAMMA2_*sec2b
	 -(beta2/18.-1./45.)*beta*sx*secb*_GAMMA1_
	 -((beta2-1.)*beta2/36.+1./420.)*sx*sec2b*_GAMMA2_
	 +(((beta2/1620.-7./3240.)*beta2+1./648.)*beta2-1./8100.)*sx2*secb*_GAMMA1_
	 +(((beta2/4536.-1./810.)*beta2+19./11340.)*beta2-13./28350.)*beta*sx2*sec2b*_GAMMA2_
	 -((((beta2/349920.-1./29160.)*beta2+71./583200.)*beta2-121./874800.)*
	   beta2+7939./224532000.)*beta*sx2*sx*secb*_GAMMA1_)*sqrt(sx)/12./sqrt(_PI_);

    return _SUCCESS_;

  }

  class_test(0==0,
	     pbs->error_message,
	     "value of l=%d or x=%e out of bounds",l,x);

}

