#include "cosmology.h"

const TransferModule& Cosmology::GetTransferModule() {
  if (computation_stage_ < cs_transfer) {
    transfer_module_ptr_ = std::unique_ptr<TransferModule>(new TransferModule(input_));
    computation_stage_ = cs_transfer;
  }
  return *transfer_module_ptr_;
}

const SpectraModule& Cosmology::GetSpectraModule() {
  if (computation_stage_ < cs_spectra) {
    spectra_module_ptr_ = std::unique_ptr<SpectraModule>(new SpectraModule(input_, GetTransferModule()));
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
