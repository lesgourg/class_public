`GW_CLASS`: Cosmological Gravitational Wave Background in `CLASS`
=================================================================

Authors:
Florian Schulze,
Lorenzo Valbusa Dall’Armi

together with Julien Lesgourgues, Angelo Ricciardone, Nicola Bartolo,
Daniele Bertacca, Christian Fidler and Sabino Matarrese.

A Jupyter notebook is provided as an example in the `notebook` folder:
[GW_CLASS-Tutorial](notebooks/GW_CLASS-Tutorial.ipynb).


Compiling CLASS and getting started
-----------------------------------

Compiling the code is identical to the basic `CLASS` code.
For more details see documentation of `CLASS`.

 1. Download the `GW_CLASS` code:
    ```
    git clone git@github.com:lesgourg/class_public.git
    cd class_public
    git checkout GW_CLASS
    ```
 2. Compile `GW_CLASS`:
    ```
    make clean; make class -j
    ```
    OR: compile `GW_CLASS` with python wrapper `classy`:
    ```
    make clean; make -j
    ```


Main instructions of `GW_CLASS`
-------------------------------

To compute the CGWB anisotropies you have to give the input:
```
output = gwCl, OmGW
```
The background GW energy density \f$ \bar{\Omega}_{\rm GW}(f) \f$ (`OmGW`) and the
angular power spectrum \f$ C_{\ell}^{\rm CGWB \times CGWB} \f$ (`gwCl`)
always have to be computed in combination.

The background GW energy density \f$ \bar{\Omega}_{\rm GW}(f) \f$ is controlled by:
```
gwb_source_type = analytic_gwb
                  OR inflationary_gwb
                  OR external_gwb
                  OR PBH_gwb
                  OR PT_gwb
n_gwb = 0.
f_pivot = 1
f_min = 1e-3
f_max = 1e2
```
and source dependent parameter.

The contributions to the CGWB anisotropies are set with:
```
gravitational_wave_contributions = ad, tsw, pisw, eisw, lisw, ini
```
New parameter involved in the computation of the CGWB anisotropies are the
for example the fraction of relativistic decoupled particles at GW production
\f$ f_{\rm dec}(\eta_{\rm in}) \f$ and the spectrum of an additional non-adiabtic mode (`gwi`) \f$ P_\Gamma^\mathrm{NAD}(k) = A_\mathrm{gwi} \, \left( \frac{k}{k_*} \right)^{n_\mathrm{gwi}} \f$
```
f_dec_ini = 0.0

ic = ad, gwi
A_gwi = 1e-10.
n_gwi = 0.
```

An exhaustive list of all parameters and description of the code can be found in the appendix of [[2305.01602](https://arxiv.org/abs/2305.01602)].


Examples
--------

The parameter file `cgwb.ini` provides an example on how to calculate the CGWB using `GW_CLASS`.

A tutorial using the python wrapper is provided in the `notebook` folder:
[GW_CLASS-Tutorial](notebooks/GW_CLASS-Tutorial.ipynb).
It shows how to reproduce figures 1-7 of [[2305.01602](https://arxiv.org/abs/2305.01602)].


Calculate SNR
-------------

`GW_CLASS` provides an external python file `exteranal/GW_SNR/compute_GW_SNR.py`,
to compute the SNR of the CGWB energy density \f$ \bar{\Omega}_{\rm GW}(f) \f$.
It can be used, to directly calculate the SNR for a CGWB computed by `GW_CLASS`:
```
python external/GW_SNR/compute_GW_SNR.py <Omega_GW_file> [h]
```


Using the code
--------------

You can use CLASS freely, provided that in your publications, you cite at least the paper
`GW_CLASS: Cosmological Gravitational Wave Background in the Cosmic Linear Anisotropy Solving System` [[2305.01602](https://arxiv.org/abs/2305.01602)].
Feel free to cite also the original CLASS papers!
