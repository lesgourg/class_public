#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "input_module.h"

class Cosmology {
public:
  Cosmology(FileContent& fc)
  : input_module_ptr_(InputModulePtr(new InputModule(fc))) {}
  Cosmology(std::unique_ptr<InputModule> input_module)
  : input_module_ptr_(std::move(input_module)) {}

  InputModulePtr& GetInputModule();
  BackgroundModulePtr& GetBackgroundModule();
  ThermodynamicsModulePtr& GetThermodynamicsModule();
  PerturbationsModulePtr& GetPerturbationsModule();
  PrimordialModulePtr& GetPrimordialModule();
  NonlinearModulePtr& GetNonlinearModule();
  TransferModulePtr& GetTransferModule();
  SpectraModulePtr& GetSpectraModule();
  LensingModulePtr& GetLensingModule();

private:
  InputModulePtr input_module_ptr_;
  BackgroundModulePtr background_module_ptr_;
  ThermodynamicsModulePtr thermodynamics_module_ptr_;
  PerturbationsModulePtr perturbations_module_ptr_;
  PrimordialModulePtr primordial_module_ptr_;
  NonlinearModulePtr nonlinear_module_ptr_;
  TransferModulePtr transfer_module_ptr_;
  SpectraModulePtr spectra_module_ptr_;
  LensingModulePtr lensing_module_ptr_;
};
#endif //COSMOLOGY_H
