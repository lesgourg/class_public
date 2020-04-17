#include "cosmology.h"

const ThermodynamicsModule& Cosmology::GetThermodynamicsModule() {
  if (!thermodynamics_module_ptr_) {
    thermodynamics_module_ptr_ = std::unique_ptr<ThermodynamicsModule>(new ThermodynamicsModule(input_));
  }
  return *thermodynamics_module_ptr_;
}

const PerturbationsModule& Cosmology::GetPerturbationsModule() {
  if (!perturbations_module_ptr_) {
    perturbations_module_ptr_ = std::unique_ptr<PerturbationsModule>(new PerturbationsModule(input_, GetThermodynamicsModule()));
  }
  return *perturbations_module_ptr_;
}

const PrimordialModule& Cosmology::GetPrimordialModule() {
  if (!primordial_module_ptr_) {
    primordial_module_ptr_ = std::unique_ptr<PrimordialModule>(new PrimordialModule(input_, GetPerturbationsModule()));
  }
  return *primordial_module_ptr_;
}

const NonlinearModule& Cosmology::GetNonlinearModule() {
  if (!nonlinear_module_ptr_) {
    nonlinear_module_ptr_ = std::unique_ptr<NonlinearModule>(new NonlinearModule(input_, GetPerturbationsModule(), GetPrimordialModule()));
  }
  return *nonlinear_module_ptr_;
}

const TransferModule& Cosmology::GetTransferModule() {
  if (!transfer_module_ptr_) {
    transfer_module_ptr_ = std::unique_ptr<TransferModule>(new TransferModule(input_, GetThermodynamicsModule(), GetPerturbationsModule(), GetNonlinearModule()));
  }
  return *transfer_module_ptr_;
}

const SpectraModule& Cosmology::GetSpectraModule() {
  if (!spectra_module_ptr_) {
    spectra_module_ptr_ = std::unique_ptr<SpectraModule>(new SpectraModule(input_, GetPerturbationsModule(), GetPrimordialModule(), GetNonlinearModule(), GetTransferModule()));
  }
  return *spectra_module_ptr_;
}

const LensingModule& Cosmology::GetLensingModule() {
  if (!lensing_module_ptr_) {
    lensing_module_ptr_ = std::unique_ptr<LensingModule>(new LensingModule(input_, GetSpectraModule()));
  }
  return *lensing_module_ptr_;
}
