/***
 * A table that  grows automatically.
 */

#ifndef __GROWTABLE__
#define __GROWTABLE__

#include "common.h"

#define _GT_INITSIZE_ 4096 /**< Init size of a growTable (in bytes)*/
#define _GT_FACTOR_ 2      /**< inflating factor when current max size is reached */
#define _GT_END_ -1        /**< flag meaning the end of the current data in the growTable */

/**
 * growTable structure.
 */
typedef struct {
  void* buffer; /**< stack of data */
  long sz;      /**< total size */
  long csz;     /**< real size */
  int freeze;   /**< if set to _TRUE_ no data can be added */

  ErrorMsg error_message; /**< error message slot */
} growTable;

/** Boilerplate for C++ */
#ifdef __cplusplus
extern "C" {
#endif

int gt_init(growTable*);

int gt_add(growTable*, long idx, void* data, long sz);
int gt_retrieve(growTable *,long idx, long sz, void* data);
int gt_retrieveAll(growTable *,void* data);

int gt_getSize(growTable*, long *idx);

int gt_getPtr(growTable*, void** ptr);

int gt_free(growTable*);

/** Boilerplate for C++ */
#ifdef __cplusplus
}
#endif
#endif
