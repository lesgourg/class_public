#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "input.h"

class Cosmology {
public:
  Cosmology(const Input& input)
  : input_(input) {}

  BackgroundModulePtr& GetBackgroundModule();
  ThermodynamicsModulePtr& GetThermodynamicsModule();
  PerturbationsModulePtr& GetPerturbationsModule();
  PrimordialModulePtr& GetPrimordialModule();
  NonlinearModulePtr& GetNonlinearModule();
  TransferModulePtr& GetTransferModule();
  SpectraModulePtr& GetSpectraModule();
  LensingModulePtr& GetLensingModule();

private:
  const Input& input_;
  BackgroundModulePtr background_module_ptr_ = nullptr;
  ThermodynamicsModulePtr thermodynamics_module_ptr_ = nullptr;
  PerturbationsModulePtr perturbations_module_ptr_ = nullptr;
  PrimordialModulePtr primordial_module_ptr_ = nullptr;
  NonlinearModulePtr nonlinear_module_ptr_ = nullptr;
  TransferModulePtr transfer_module_ptr_ = nullptr;
  SpectraModulePtr spectra_module_ptr_ = nullptr;
  LensingModulePtr lensing_module_ptr_ = nullptr;
};
#endif //COSMOLOGY_H
