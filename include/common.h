/** @file common.h Generic libraries, parameters and functions used in the whole code. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#ifndef __COMMON__
#define __COMMON__

#define _TRUE_ 1 /**< integer associated to true statement */
#define _FALSE_ 0 /**< integer associated to false statement */

#define _SUCCESS_ 0 /**< integer returned after successfull call of a function */
#define _FAILURE_ 1 /**< integer returned after failure in a function */

#define _ERRORMSGSIZE_ 2048 /**< generic error messages are cut beyond this number of characters */
typedef char ErrorMsg[_ERRORMSGSIZE_]; /**< Generic error messages (there is such a field in each structure) */

#define _PI_ 3.1415926535897932384626433832795e0 /**< The number pi */

#define _Mpc_over_m_ 3.085678e22 /**< conversion factor from meters to megaparsecs */
#define _Gyr_over_Mpc_ 306.61 /**< conversion factor from megaparsecs to gigayears (c=1 units) */

#define _MAX_IT_ 10000/**< default maximum number of iterations in conditional loops (to avoid infinite loops) */

#define min(a,b) (((a)<(b)) ? (a) : (b) ) /**< the usual "min" function */
#define max(a,b) (((a)<(b)) ? (b) : (a) ) /**< the usual "max" function */

#endif
