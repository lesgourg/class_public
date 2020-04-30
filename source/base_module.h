#ifndef BASE_MODULE_H
#define BASE_MODULE_H

/* class modules */
#include "common.h"
#include "input_module.h"
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
  BaseModule(InputModulePtr input_module)
  : ppr(&input_module->precision_)
  , pba(&input_module->background_)
  , pth(&input_module->thermodynamics_)
  , ppt(&input_module->perturbations_)
  , ppm(&input_module->primordial_)
  , psp(&input_module->spectra_)
  , pnl(&input_module->nonlinear_)
  , ptr(&input_module->transfers_)
  , ple(&input_module->lensing_)
  , pop(&input_module->output_) {
    input_module_ = std::move(input_module);
  }
  BaseModule(const BaseModule&) = delete;
  
  mutable ErrorMsg error_message_;
protected:
  InputModulePtr input_module_;

  const precision* const ppr;
  const background* const pba;
  const thermo* const pth;
  const perturbs* const ppt;
  const primordial* const ppm;
  const nonlinear* const pnl;
  const transfers* const ptr;
  const spectra* const psp;
  const lensing* const ple;
  const output* const pop;
};


#endif //BASE_MODULE_H
