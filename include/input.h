/** @file input.h Documented includes for input module */

#ifndef __INPUT__
#define __INPUT__

#include "common.h"
#include "parser.h"

/* Declare the structs this module refers to */
struct precision * ppr;
struct background * pba;
struct thermo * pth;
struct perturbs * ppt;
struct transfers * ptr;
struct primordial * ppm;
struct spectra * psp;
struct nonlinear * pnl;
struct lensing * ple;
struct output * pop;
struct distortions * psd;

#define _N_FILEROOT_ 100 /* Number of files that will be not overwritten for a given root */

/* macro for checking if string begins with certain character */
int string_begins_with(char* thestring,char beginchar);

/* macro for reading parameter values with routines from the parser */
#define class_read_double(name,destination)                                     \
  do {                                                                          \
    double param_temp; int flag_temp;                                           \
    class_call(parser_read_double(pfc,name,&param_temp,&flag_temp,errmsg),      \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == _TRUE_){                                                   \
      destination = param_temp;                                                 \
    }                                                                           \
  } while(0);


#define class_read_int(name,destination)                                        \
  do {                                                                          \
    int int_temp,flag_temp;                                                     \
    class_call(parser_read_int(pfc,name,&int_temp,&flag_temp,errmsg),           \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == _TRUE_){                                                   \
      destination = int_temp;                                                   \
    }                                                                           \
  } while(0);

#define class_read_string(name,destination)                                     \
  do {                                                                          \
    char string_temp[_ARGUMENT_LENGTH_MAX_]; int flag_temp;                     \
    class_call(parser_read_string(pfc,name,&string_temp,&flag_temp,errmsg),     \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == _TRUE_){                                                   \
      strcpy(destination,string_temp);                                          \
    }                                                                           \
  } while(0);

#define class_read_flag(name,destination)                                       \
  do {                                                                          \
    char string_temp[_ARGUMENT_LENGTH_MAX_]; int flag_temp;                     \
    class_call(parser_read_string(pfc,name,&string_temp,&flag_temp,errmsg),     \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == _TRUE_){                                                   \
      if( string_begins_with(string_temp,'y')                                   \
         || string_begins_with(string_temp,'Y') ){                              \
        destination = _TRUE_;                                                   \
      }                                                                         \
      if( string_begins_with(string_temp,'n')                                   \
         || string_begins_with(string_temp,'N') ){                              \
        destination = _FALSE_;                                                  \
      }                                                                         \
      else {                                                                    \
        class_stop(errmsg,"incomprehensible input '%s' for the field '%s'.",    \
                   string_temp, name);                                          \
      }                                                                         \
    }                                                                           \
  } while(0);

#define class_read_flag_or_deprecated(name,oldname,destination)                 \
  do {                                                                          \
    char string_temp[_ARGUMENT_LENGTH_MAX_]; int flag_temp;                     \
    class_call(parser_read_string(pfc,name,&string_temp,&flag_temp,errmsg),     \
               errmsg,                                                          \
               errmsg);                                                         \
    /* Compatibility code BEGIN */                                              \
    if(flag_temp == _FALSE_){                                                   \
      class_call(parser_read_string(pfc,oldname,&string_temp,&flag_temp,errmsg),\
                 errmsg,                                                        \
                 errmsg);                                                       \
    }                                                                           \
    /* Compatibility code END */                                                \
    if (flag_temp == _TRUE_){                                                   \
      if( string_begins_with(string_temp,'y')                                   \
         || string_begins_with(string_temp,'Y') ){                              \
        destination = _TRUE_;                                                   \
      }                                                                         \
      if( string_begins_with(string_temp,'n')                                   \
         || string_begins_with(string_temp,'N') ){                              \
        destination = _FALSE_;                                                  \
      }                                                                         \
      else {                                                                    \
        class_stop(errmsg,"incomprehensible input '%s' for the field '%s'.",    \
                   string_temp, name);                                          \
      }                                                                         \
    }                                                                           \
  } while(0);

#define class_read_double_one_of_two(name1,name2,destination)                   \
  do {                                                                          \
    int flag_temp1,flag_temp2;                                                  \
    double param_temp1,param_temp2;                                             \
    class_call(parser_read_double(pfc,name1,&param_temp1,&flag_temp1,errmsg),   \
               errmsg,                                                          \
               errmsg);                                                         \
    class_call(parser_read_double(pfc,name2,&param_temp2,&flag_temp2,errmsg),   \
               errmsg,                                                          \
               errmsg);                                                         \
    class_test((flag_temp1 == _TRUE_) && (flag_temp2 == _TRUE_),                \
               errmsg,                                                          \
               "You can only enter one of '%s' or '%s'.",                       \
               name1,name2);                                                    \
    if (flag_temp1 == _TRUE_){                                                  \
      destination = param_temp1;                                                \
    }                                                                           \
    if (flag_temp2 == _TRUE_){                                                  \
      destination = param_temp2;                                                \
    }                                                                           \
  } while(0);

#define class_at_least_two_of_three(a,b,c)                                      \
  (((a) == _TRUE_) && ((b) == _TRUE_)) ||                                       \
  (((a) == _TRUE_) && ((c) == _TRUE_)) ||                                       \
  (((b) == _TRUE_) && ((c) == _TRUE_))

#define class_none_of_three(a,b,c)                                              \
  ((a) == _FALSE_) && ((b) == _FALSE_) && ((c) == _FALSE_)

#define class_read_list_of_doubles_or_default(name,destination,val_default,siz) \
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_doubles(pfc,name,                            \
                                           &entries_read_temp,                  \
                                           &(destination),                      \
                                           &flag_temp,                          \
                                           errmsg),                             \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == _TRUE_){                                                   \
      class_test(entries_read_temp != siz, errmsg,                              \
                 "Number of entries of '%s' (%d) does not match expected number (%d).", \
                 name, entries_read_temp, siz);                                 \
    }else{                                                                      \
      class_alloc(destination,siz*sizeof(double),errmsg);                       \
      for(n=0; n<siz; n++){destination[n] = val_default;}                       \
    }                                                                           \
  } while(0);

#define class_read_list_of_integers_or_default(name,destination,val_default,siz)\
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_integers(pfc,name,                           \
                                            &entries_read_temp,                 \
                                            &(destination),                     \
                                            &flag_temp,                         \
                                            errmsg),                            \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == _TRUE_){                                                   \
      class_test(entries_read_temp != siz, errmsg,                              \
                 "Number of entries of '%s' (%d) does not match expected number (%d).", \
                 name, entries_read_temp, siz);                                 \
    }else{                                                                      \
      class_alloc(destination,siz*sizeof(int),errmsg);                          \
      for(n=0; n<siz; n++){destination[n] = val_default;}                       \
    }                                                                           \
  } while(0);

#define class_read_list_of_doubles(name,destination,siz)                        \
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_doubles(pfc,name,                            \
                                           &entries_read_temp,                  \
                                           &(destination),                      \
                                           &flag_temp,                          \
                                           errmsg),                             \
               errmsg,                                                          \
               errmsg);                                                         \
    class_test(flag_temp == _FALSE_, errmsg,                                    \
               "Entry '%s' is required but not found!", name)                   \
    class_test(entries_read_temp != siz, errmsg,                                \
               "Number of entries of '%s' (%d) does not match expected number (%d).", \
               name,entries_read_temp, siz);                                    \
  } while(0);

#define class_read_list_of_integers(name,destination,siz)                       \
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_integers(pfc,name,                           \
                                            &entries_read_temp,                 \
                                            &(destination),                     \
                                            &flag_temp,                         \
                                            errmsg),                            \
               errmsg,                                                          \
               errmsg);                                                         \
    class_test(flag_temp == _FALSE_, errmsg,                                    \
               "Entry '%s' is required but not found!", name)                   \
    class_test(entries_read_temp != siz, errmsg,                                \
               "Number of entries of '%s' (%d) does not match expected number (%d).", \
               name,entries_read_temp, siz);                                    \
  } while(0);

/**
 * temporary parameters for background fzero function
 */

#define _NUM_TARGETS_ 7 //Keep this number as number of target_names
enum target_names {theta_s, Omega_dcdmdr, omega_dcdmdr, Omega_scf, Omega_ini_dcdm, omega_ini_dcdm, sigma8};
enum computation_stage {cs_background, cs_thermodynamics, cs_perturbations, cs_primordial, cs_nonlinear, cs_transfer, cs_spectra};

struct input_pprpba{
  struct precision * ppr;
  struct background * pba;
};

struct fzerofun_workspace {
  int * unknown_parameters_index;
  struct file_content fc;
  enum target_names * target_name;
  double * target_value;
  int target_size;
  enum computation_stage required_computation_stage;
};

/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

  /* Main functions */
  int input_init(int argc,
                 char **argv,
                 struct precision * ppr,
                 struct background * pba,
                 struct thermo * pth,
                 struct perturbs * ppt,
                 struct transfers * ptr,
                 struct primordial * ppm,
                 struct spectra * psp,
                 struct nonlinear * pnl,
                 struct lensing *ple,
                 struct distortions *psd,
                 struct output *pop,
                 ErrorMsg errmsg);

  int input_find_file(int argc,
                      char ** argv,
                      struct file_content * fc,
                      ErrorMsg errmsg);

  int input_set_root(char* input_file,
                     struct file_content** ppfc_input,
                     struct file_content* pfc_setroot,
                     ErrorMsg errmsg);

  int file_exists(const char *fname);

  int input_read_from_file(struct file_content * pfc,
                           struct precision * ppr,
                           struct background *pba,
                           struct thermo *pth,
                           struct perturbs *ppt,
                           struct transfers *ptr,
                           struct primordial *ppm,
                           struct spectra *psp,
                           struct nonlinear *pnl,
                           struct lensing *ple,
                           struct distortions *psd,
                           struct output *pop,
                           ErrorMsg errmsg);

  /* Shooting */
  int input_shooting(struct file_content * pfc,
                     struct precision * ppr,
                     struct background * pba,
                     struct thermo * pth,
                     struct perturbs * ppt,
                     struct transfers * ptr,
                     struct primordial * ppm,
                     struct spectra * psp,
                     struct nonlinear * pnl,
                     struct lensing *ple,
                     struct distortions *psd,
                     struct output *pop,
                     int input_verbose,
                     int * has_shooting,
                     ErrorMsg errmsg);

  int input_needs_shooting_for_target(struct file_content * pfc,
                                      enum target_names target_name,
                                      double target_value,
                                      int * aux_flag,
                                      ErrorMsg errmsg);

  int input_find_root(double * xzero,
                      int * fevals,
                      struct fzerofun_workspace * pfzw,
                      ErrorMsg errmsg);

  int input_fzerofun_1d(double input,
                        void * fzerofun_workspace,
                        double * output,
                        ErrorMsg error_message);

  int class_fzero_ridder(int (*func)(double x,
                                     void * param,
                                     double * y,
                                     ErrorMsg error_message),
                         double x1,
                         double x2,
                         double xtol,
                         void * param,
                         double * Fx1,
                         double * Fx2,
                         double * xzero,
                         int * fevals,
                         ErrorMsg error_message);

  int input_get_guess(double * xguess,
                      double * dxdy,
                      struct fzerofun_workspace * pfzw,
                      ErrorMsg errmsg);

  int input_try_unknown_parameters(double * unknown_parameter,
                                   int unknown_parameters_size,
                                   void * pfzw,
                                   double * output,
                                   ErrorMsg errmsg);

  /* Read from precision.h */
  int input_read_precisions(struct file_content * pfc,
                            struct precision * ppr,
                            struct background * pba,
                            struct thermo * pth,
                            struct perturbs * ppt,
                            struct transfers * ptr,
                            struct primordial * ppm,
                            struct spectra * psp,
                            struct nonlinear * pnl,
                            struct lensing * ple,
                            struct distortions *psd,
                            struct output * pop,
                            ErrorMsg errmsg);

  /* Read from .ini file */
  int input_read_parameters(struct file_content * pfc,
                            struct precision * ppr,
                            struct background * pba,
                            struct thermo * pth,
                            struct perturbs * ppt,
                            struct transfers * ptr,
                            struct primordial * ppm,
                            struct spectra * psp,
                            struct nonlinear * pnl,
                            struct lensing * ple,
                            struct distortions *psd,
                            struct output * pop,
                            ErrorMsg errmsg);

  int input_read_parameters_general(struct file_content * pfc,
                                    struct background * pba,
                                    struct thermo * pth,
                                    struct perturbs * ppt,
                                    struct distortions * psd,
                                    ErrorMsg errmsg);

  int input_read_parameters_species(struct file_content * pfc,
                                    struct precision * ppr,
                                    struct background * pba,
                                    struct thermo * pth,
                                    struct perturbs * ppt,
                                    int input_verbose,
                                    ErrorMsg errmsg);

  int input_read_parameters_heating(struct file_content * pfc,
                                    struct thermo * pth,
                                    ErrorMsg errmsg);

  int input_read_parameters_nonlinear(struct file_content * pfc,
                                      struct precision * ppr,
                                      struct perturbs * ppt,
                                      struct nonlinear * pnl,
                                      int input_verbose,
                                      ErrorMsg errmsg);

  int input_prepare_pk_eq(struct precision * ppr,
                          struct background * pba,
                          struct thermo * pth,
                          struct nonlinear * pnl,
                          int input_verbose,
                          ErrorMsg errmsg);

  int input_read_parameters_primordial(struct file_content * pfc,
                                       struct perturbs * ppt,
                                       struct primordial * ppm,
                                       ErrorMsg errmsg);

  int input_read_parameters_spectra(struct file_content * pfc,
                                    struct precision * ppr,
                                    struct background * pba,
                                    struct perturbs * ppt,
                                    struct transfers * ptr,
                                    struct spectra * psp,
                                    struct output * pop,
                                    ErrorMsg errmsg);

  int input_read_parameters_lensing(struct file_content * pfc,
                                    struct precision * ppr,
                                    struct perturbs * ppt,
                                    struct transfers * ptr,
                                    struct lensing * ple,
                                    ErrorMsg errmsg);

  int input_read_parameters_distortions(struct file_content * pfc,
                                        struct precision * ppr,
                                        struct distortions * psd,
                                        ErrorMsg errmsg);

  int input_read_parameters_additional(struct file_content * pfc,
                                       struct precision * ppr,
                                       struct background * pba,
                                       struct thermo * pth,
                                       ErrorMsg errmsg);

  int input_read_parameters_output(struct file_content * pfc,
                                   struct background * pba,
                                   struct thermo * pth,
                                   struct perturbs * ppt,
                                   struct transfers * ptr,
                                   struct primordial * ppm,
                                   struct spectra * psp,
                                   struct nonlinear * pnl,
                                   struct lensing *ple,
                                   struct distortions *psd,
                                   struct output *pop,
                                   ErrorMsg errmsg);

  int compare_doubles(const void * a,
                      const void * b);

  /* Set default parameters */
  int input_default_params(struct background *pba,
                           struct thermo *pth,
                           struct perturbs *ppt,
                           struct transfers *ptr,
                           struct primordial *ppm,
                           struct spectra *psp,
                           struct nonlinear *pnl,
                           struct lensing *ple,
                           struct distortions *psd,
                           struct output *pop);


#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
/* @endcond */
