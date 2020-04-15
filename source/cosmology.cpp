#include "cosmology.h"

const PrimordialModule& Cosmology::GetPrimordialModule() {
  if (computation_stage_ < cs_primordial) {
    primordial_module_ptr_ = std::unique_ptr<PrimordialModule>(new PrimordialModule(input_));
    computation_stage_ = cs_primordial;
  }
  return *primordial_module_ptr_;
}

const NonlinearModule& Cosmology::GetNonlinearModule() {
  if (computation_stage_ < cs_nonlinear) {
    nonlinear_module_ptr_ = std::unique_ptr<NonlinearModule>(new NonlinearModule(input_, GetPrimordialModule()));
    computation_stage_ = cs_nonlinear;
  }
  return *nonlinear_module_ptr_;
}

const TransferModule& Cosmology::GetTransferModule() {
  if (computation_stage_ < cs_transfer) {
    transfer_module_ptr_ = std::unique_ptr<TransferModule>(new TransferModule(input_, GetNonlinearModule()));
    computation_stage_ = cs_transfer;
  }
  return *transfer_module_ptr_;
}

const SpectraModule& Cosmology::GetSpectraModule() {
  if (computation_stage_ < cs_spectra) {
    spectra_module_ptr_ = std::unique_ptr<SpectraModule>(new SpectraModule(input_, GetPrimordialModule(), GetNonlinearModule(), GetTransferModule()));
    computation_stage_ = cs_spectra;
  }
  return *spectra_module_ptr_;
}

const LensingModule& Cosmology::GetLensingModule() {
  if (computation_stage_ < cs_lensing) {
    lensing_module_ptr_ = std::unique_ptr<LensingModule>(new LensingModule(input_, GetSpectraModule()));
    computation_stage_ = cs_lensing;
  }
  return *lensing_module_ptr_;
}
