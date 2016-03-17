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

  FileName root; /**< root for all file names */

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

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  int output_total_cl_at_l(
                           struct spectra * psp,
                           struct lensing * ple,
                           struct output * pop,
                           int l,
                           double * cl
                           );

  int output_init(
                  struct background * pba,
                  struct thermo * pth,
                  struct perturbs * ppt,
                  struct primordial * ppm,
                  struct transfers * ptr,
                  struct spectra * psp,
                  struct nonlinear * pnl,
                  struct lensing * ple,
                  struct output * pop
                  );

  int output_cl(
                struct background * pba,
                struct perturbs * ppt,
                struct spectra * psp,
                struct lensing * ple,
                struct output * pop
                );

  int output_pk(
                struct background * pba,
                struct perturbs * ppt,
                struct spectra * psp,
                struct output * pop
                );

  int output_pk_nl(
                   struct background * pba,
                   struct perturbs * ppt,
                   struct spectra * psp,
                   struct output * pop
                   );

  int output_tk(
                struct background * pba,
                struct perturbs * ppt,
                struct spectra * psp,
                struct output * pop
                );

  int output_background(
                        struct background * pba,
                        struct output * pop
                        );

  int output_thermodynamics(
                            struct background * pba,
                            struct thermo * pth,
                            struct output * pop
                            );

  int output_perturbations(
                           struct background * pba,
                           struct perturbs * ppt,
                           struct output * pop
                           );

  int output_primordial(
                        struct perturbs * ppt,
                        struct primordial * ppm,
                        struct output * pop
                        );

  int output_print_data(FILE *out,
                        char titles[_MAXTITLESTRINGLENGTH_],
                        double *dataptr,
                        int tau_size);
  int output_open_cl_file(
                          struct spectra * psp,
                          struct output * pop,
                          FILE ** clfile,
                          FileName filename,
                          char * first_line,
                          int lmax
                          );

  int output_one_line_of_cl(
                            struct background * pba,
                            struct spectra * psp,
                            struct output * pop,
                            FILE * clfile,
                            double l,
                            double * cl,
                            int ct_size
                            );

  int output_open_pk_file(
                          struct background * pba,
                          struct spectra * psp,
                          struct output * pop,
                          FILE ** pkfile,
                          FileName filename,
                          char * first_line,
                          double z
                          );

  int output_one_line_of_pk(
                            FILE * tkfile,
                            double one_k,
                            double one_pk
                            );

  int output_open_pk_nl_file(
                             struct background * pba,
                             struct nonlinear * pnl,
                             struct output * pop,
                             FILE ** pkfile,
                             FileName filename,
                             char * first_line,
                             double z,
                             int k_size
                             );


#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
