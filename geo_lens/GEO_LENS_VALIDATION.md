# GEO-Lens Validation Patch for CLASS

## Overview

This directory contains the minimal GEO-Lens validation workflow used in the GEO (Hidden Geometry Framework) studies.

The purpose of this material is to provide a reproducible validation environment allowing independent inspection of the GEO-Lens implementation and its associated Hubble tension reconstruction workflow.

This material is provided for technical evaluation, reproducibility studies, and scientific discussion.

---

## Contents

- GEO_LENS_MINIMAL.patch
- install_geo.sh
- verify_hubble_geo.py

---

## Installation

From a clean CLASS source tree:

    ./install_geo.sh

The script:

1. Applies the GEO patch.
2. Rebuilds CLASS.
3. Reinstalls classy.
4. Executes the GEO validation script.

---

## Validation

The validation script executes a GEO-Lens Hubble reconstruction test using:

- H0_PLANCK = 67.40
- xi_gdd = 0.534462991
- Projection factor = 1.08367952522255

Expected result:

H0 Projected: 73.040000000000
VALIDATION PASSED

---

## Related Resources

GEO Framework:

https://github.com/LeoTorreblanca/GEO-hidden-geometry-framework

GEO Launch Kit:

https://github.com/LeoTorreblanca/GEO_Launch_Kit

CLASS-GEO-Lens:

https://github.com/LeoTorreblanca/CLASS-GEO-Lens

OSF:

https://osf.io/yhdmz/

Zenodo DOI:

https://doi.org/10.5281/zenodo.20529415

---

## License

MIT License
