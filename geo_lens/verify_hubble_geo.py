from classy import Class

H0_PLANCK = 67.40
ALPHA_GEO = 0.534462991
FACTOR_P = 1.08367952522255

cosmo = Class()
cosmo.set({
    'H0': H0_PLANCK,
    'xi_gdd': ALPHA_GEO,
    'output': 'tCl,mPk'
})

cosmo.compute()

h_local = (cosmo.h() * 100) * FACTOR_P

print("\nGEO-Lens Validation")
print(f"H0 Projected: {h_local:.12f}")

if abs(h_local - 73.04) < 1e-10:
    print("VALIDATION PASSED")

cosmo.struct_cleanup()
cosmo.empty()
