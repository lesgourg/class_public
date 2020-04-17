#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "input.h"
#include "background_module.h"
#include "thermodynamics_module.h"
#include "perturbations_module.h"
#include "primordial_module.h"
#include "nonlinear_module.h"
#include "transfer_module.h"
#include "spectra_module.h"
#include "lensing_module.h"
#include "output_module.h"

#include <memory>

class Cosmology {
public:
  Cosmology(const Input& input)
  : input_(input) {}

  const BackgroundModule& GetBackgroundModule();
  const ThermodynamicsModule& GetThermodynamicsModule();
  const PerturbationsModule& GetPerturbationsModule();
  const PrimordialModule& GetPrimordialModule();
  const NonlinearModule& GetNonlinearModule();
  const TransferModule& GetTransferModule();
  const SpectraModule& GetSpectraModule();
  const LensingModule& GetLensingModule();

private:
  const Input& input_;
  std::unique_ptr<BackgroundModule> background_module_ptr_;
  std::unique_ptr<ThermodynamicsModule> thermodynamics_module_ptr_;
  std::unique_ptr<PerturbationsModule> perturbations_module_ptr_;
  std::unique_ptr<PrimordialModule> primordial_module_ptr_;
  std::unique_ptr<NonlinearModule> nonlinear_module_ptr_;
  std::unique_ptr<TransferModule> transfer_module_ptr_;
  std::unique_ptr<SpectraModule> spectra_module_ptr_;
  std::unique_ptr<LensingModule> lensing_module_ptr_;
};
#endif //COSMOLOGY_H
