#!/bin/bash

echo "Installing GEO-Lens validation patch"

patch -p1 --forward --batch --force < geo_lens/GEO_HUBBLE_FINAL.patch

echo "Building CLASS..."
make clean
make -j

echo "Installing classy..."
python3 -m pip install ./python --user

echo "Running GEO-Lens validation..."
python3 geo_lens/verify_hubble_geo.py
