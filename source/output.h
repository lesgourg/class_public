/** @file output.h Documented includes for output module */

#ifndef __OUTPUT__
#define __OUTPUT__

#include "common.h"
#include "lensing.h"

/**
 * Maximum number of values of redshift at which the spectra will be
 * written in output files
 */

#define _Z_PK_NUM_MAX_ 100

/**
 * Structure containing various informations on the output format,
 * all of them initialized by user in input module.
 *
 */

struct output {

  //@{

  char root[_FILENAMESIZE_-32]; /**< root for all file names */

  //@}

  /** @name - number and value(s) of redshift at which P(k,z) and T_i(k,z) should be written */

  //@{

  int z_pk_num; /**< number of redshift at which P(k,z) and T_i(k,z) should be written */
  double z_pk[_Z_PK_NUM_MAX_]; /**< value(s) of redshift at which P(k,z) and T_i(k,z) should be written */

  //@}

  /** @name - extra information on output */

  //@{

  short write_header; /**< flag stating whether we should write a header in output files */

  enum file_format output_format; /**< which format for output files (definitions, order of columns, etc.) */

  short write_background; /**< flag for outputing background evolution in file */
  short write_thermodynamics; /**< flag for outputing thermodynamical evolution in file */
  short write_perturbations; /**< flag for outputing perturbations of selected wavenumber(s) in file(s) */
  short write_primordial; /**< flag for outputing scalar/tensor primordial spectra in files */

  //@}

  /** @name - technical parameters */

  //@{

  short output_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  //@}
};

#endif
