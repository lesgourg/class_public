#ifndef __HEATING__
#define __HEATING__

struct heating{

  double* z_table;
  int z_size;

  ErrorMsg error_message;
};



/**************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif
int heating_init(struct thermo* pth, struct heating* phe);

#ifdef __cplusplus
}
#endif

#endif
