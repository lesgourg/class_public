
/**** Cosmological parameters.
      Include information on starting and ending redshit and timestep  ****/

typedef struct {
   double T0;                   /* CMB temperature today in K*/
   double obh2, omh2, okh2;     /* cosmological parameters */
   double odeh2, w0, wa;        /* dark energy parameters */
   double Y;                    /* primordial helium abundance */
   double Nnueff;               /* effective number of neutrinos */

   /* Secondary parameters, to avoid recalculating every time */
   double nH0;                  /* density of hydrogen today in m^{-3} */
   double fHe;                  /* Helium fraction by number */

   double zstart, zend, dlna;   /* initial and final redshift and step size in log a */
   long nz;                     /* total number of redshift steps */

   /** parameters for energy injection */

   double annihilation; /** parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003) */

   short has_on_the_spot; /** do we want to use the on-the-spot approximation? */

   double decay; /** parameter descibing CDM decay (f/tau, see e.g. 1109.6322)*/

   double annihilation_variation; /** if this parameter is non-zero,
				     the function F(z)=(f <sigma*v> /
				     m_cdm)(z) will be a parabola in
				     log-log scale between zmin and
				     zmax, with a curvature given by
				     annihlation_variation (must ne
				     negative), and with a maximum in
				     zmax; it will be constant outside
				     this range */
   int annihil_coef_num_lines;
   double *annihil_coef_heat;
   double *annihil_coef_ionH;
   double *annihil_coef_ionHe;
   double *annihil_coef_lya;
   double *annihil_coef_lowE;
   double *annihil_coef_xe;
   double *annihil_coef_dd_heat;
   double *annihil_coef_dd_ionH;
   double *annihil_coef_dd_ionHe;
   double *annihil_coef_dd_lya;
   double *annihil_coef_dd_lowE;


   double annihilation_z; /** if annihilation_variation is non-zero,
			     this is the value of z at which the
			     parameter annihilation is defined, i.e.
			     F(annihilation_z)=annihilation */

   double annihilation_zmax; /** if annihilation_variation is non-zero,
				redhsift above which annihilation rate
				is maximal */

   double annihilation_zmin; /** if annihilation_variation is non-zero,
				redhsift below which annihilation rate
				is constant */

   double annihilation_f_halo; /* takes the contribution of DM annihilation in halos into account*/
   double annihilation_z_halo; /*characteristic redshift for DM annihilation in halos*/

} REC_COSMOPARAMS;
