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
    
  self->buffer = malloc(_GT_INITSIZE_);
  if (self->buffer==NULL) {
    sprintf(self->error_message,"%s(L:%d): Cannot allocate growTable->buffer \n",__func__,__LINE__);
    return _FAILURE_;
  }
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
  
  if (self->freeze == _TRUE_) {
    sprintf(self->error_message,"%s(L:%d): cannot add any more data in the growTable (freeze is on)\n",__func__,__LINE__);
    return _FAILURE_;
  }
  
  if (idx==_GT_END_) {
    ridx=self->csz;
  }
  else {
    ridx=idx;
  }
  if (ridx<0) {
    sprintf(self->error_message,"%s(L:%d): Don't know what to do with idx=%ld\n",__func__,__LINE__,ridx);
    return _FAILURE_;
  }
 
  if (ridx+sz>self->sz) {
    /** - test -> pass -> ok we need to grow */
    nbuffer=realloc(self->buffer,self->sz*_GT_FACTOR_);
    if (nbuffer==NULL) {
      sprintf(self->error_message,"%s(L:%d): Cannot grow growTable\n",__func__,__LINE__);
      return _FAILURE_;
    }
    self->buffer=nbuffer;
    self->sz=self->sz*_GT_FACTOR_;
  }
  
  res=memcpy((void*) (self->buffer+ridx),(void*) data,(size_t) sz);
  if (res!=self->buffer+ridx) {
    sprintf(self->error_message,"%s(L:%d): Cannot add data to growTable\n",__func__,__LINE__);
    return _FAILURE_;  
  }
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
  
  if (idx<0) {
    sprintf(self->error_message,"%s(L:%d): Don't know what to do with idx=%ld\n",__func__,__LINE__,idx);
    return _FAILURE_;
  }

  if ((idx>self->csz) || (idx+sz>self->csz)) {
    sprintf(self->error_message,"%s(L:%d): Not enough data in growTable (asked for [%ld,%ld[ got only  %ld data)\n",__func__,__LINE__);
    return _FAILURE_;  
  }
  
  res=memcpy(data,self->buffer+idx,sz);
  if (res!=self->buffer+idx) {
    sprintf(self->error_message,"%s(L:%d): Cannot retrieve data from the growTable\n",__func__,__LINE__);
    return _FAILURE_; 
  }
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
  if (self->csz<0) {
    sprintf(self->error_message,"%s(L:%d): growTable does not make sense\n",__func__,__LINE__);
    return _FAILURE_;
  }
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
}

