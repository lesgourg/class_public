#ifndef BASE_MODULE_H
#define BASE_MODULE_H

/* class modules */
#include "common.h"
#include "input.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "primordial.h"
#include "nonlinear.h"
#include "transfer.h"
#include "spectra.h"
#include "lensing.h"
#include "output.h"


class BaseModule {
public:
  BaseModule(const Input& input)
  : ppr(&input.precision_)
  , pba(&input.background_)
  , pth(&input.thermodynamics_)
  , ppt(&input.perturbations_)
  , ppm(&input.primordial_)
  , psp(&input.spectra_)
  , pnl(&input.nonlinear_)
  , ptr(&input.transfers_)
  , ple(&input.lensing_)
  , pop(&input.output_) {}
  BaseModule(const BaseModule&) = delete;
  
  mutable ErrorMsg error_message_;
protected:
  const precision * const ppr;
  const background * const pba;
  const thermo * const pth;
  const perturbs * const ppt;
  const primordial * const ppm;
  const nonlinear * const pnl;
  const transfers * const ptr;
  const spectra * const psp;
  const lensing * const ple;
  const output * const pop;
};


#endif //BASE_MODULE_H
