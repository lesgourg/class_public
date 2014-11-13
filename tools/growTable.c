/***
 * growTable provides automatically growing tables.
 */

#include "growTable.h"

/***
 * gt_init Initialize the growTable.
 * gt_init will initialize the growTable structure. It must be already allocated.
 *
 * Called by background_solve().
 */
int gt_init(
  growTable* self /***< a pointer on an empty growTable */
  ) {

  class_alloc(self->buffer,_GT_INITSIZE_,self->error_message);
  self->sz=_GT_INITSIZE_;
  self->csz=0;
  self->freeze=_FALSE_;  /**< This line added by JL */
  return _SUCCESS_;
}

/**
 * Add data to the growTable.
 *
 * Called by background_solve().
 */
int gt_add(
  growTable* self, /**< a growTable*/
  long idx,        /**< index at wich to add the data (in bytes). #_GT_END_ means the end of the currently written data*/
  void* data,      /**< data to be added*/
  long sz          /**< size of the data (in bytes)*/
  ) {
  long ridx;
  void *res;
  void *nbuffer;

  /** - assumes the growTable is correctly initialized */

  class_test(self->freeze == _TRUE_,
	     self->error_message,
	     "cannot add any more data in the growTable (freeze is on)");

  if (idx==_GT_END_) {
    ridx=self->csz;
  }
  else {
    ridx=idx;
  }
  class_test(ridx<0,
	     self->error_message,
	     "Don't know what to do with idx=%ld",ridx);

  if (ridx+sz>self->sz) {
    /** - test -> pass -> ok we need to grow */
    nbuffer=realloc(self->buffer,self->sz*_GT_FACTOR_);
    class_test(nbuffer==NULL,
	       self->error_message,
	       "Cannot grow growTable");
    self->buffer=nbuffer;
    self->sz=self->sz*_GT_FACTOR_;
  }

  res=memcpy((void*) (self->buffer+ridx),(void*) data,(size_t) sz);
  class_test(res!=self->buffer+ridx,
	     self->error_message,
	     "Cannot add data to growTable");
  self->csz=ridx+sz;

  return _SUCCESS_;
}

/**
 * Retrieve data from the growTable.
 *
 * Not called.
 */
int gt_retrieve(
  growTable *self, /**< a growTable*/
  long idx,        /**< index at wich to retrieve the data (in bytes).*/
  long sz,         /**< size of the data (in bytes)*/
  void* data       /**< OUTPUT : data must be allocated to ::sz bytes*/
  ) {
  void *res;

  class_test(idx<0,
	     self->error_message,
	     "don't know what to do with idx=%ld",idx);

  class_test((idx>self->csz) || (idx+sz>self->csz),
	     self->error_message,
	     "not enough data in growTable");

  res=memcpy(data,self->buffer+idx,sz);
  class_test(res!=self->buffer+idx,
	     self->error_message,
	     "cannot retrieve data from the growTable");

  return _SUCCESS_;
}

/**
 * Retrieve all data from the growTable.
 *
 * Not called.
 */
int gt_retrieveAll(
  growTable *self, /**< a growTable*/
  void* data       /**< OUTPUT : data must be allocated to the size of the growTable (see gt_getSize)*/
  ) {
  return gt_retrieve(self,0,self->csz,data);
}

/**
 * returns the size of the growTable
 *
 *  Not called.
 */
int gt_getSize(
  growTable* self,/**< a growTable*/
  long *idx /**< OUTPUT : the size of the growTable ::self*/
  ) {
  class_test(self->csz<0,
	     self->error_message,
	     "growTable does not make sense");
  *idx=self->csz;
  return _SUCCESS_;
}

/**
 * returns a pointer on the data contained in the growTable.
 * No Data can be added afterward !!!!! This is not for the faint of heart.
 *
 * Called by background_solve().
 */
int gt_getPtr(
  growTable* self, /**< a growTable*/
  void** ptr       /**< OUTPUT : pointer on the data */
  ) {
  self->freeze=_TRUE_;
  *ptr=self->buffer;

  return _SUCCESS_;
}


/**
 * free the growTable
 *
 * Called by background_solve().
 */
int gt_free(growTable* self) {
  free(self->buffer);
  self->csz=-1;
  self->sz=-1;
  self->freeze=_FALSE_;  /**< This line added by JL */

  return _SUCCESS_;
}

