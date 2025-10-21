#!/usr/bin/env bash
set -euo pipefail
awk 'NR==1{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\tD_L_Mpc\tE_z\tmu_mag";next} NR==2{H0=$5} NR>1{Dl=(1+$1)^2*$3; Ez=$5/H0; mu=5*log(Dl)/log(10)+25; printf "%s\t%.9f\t%.9f\t%.6f\n",$1"\t"$2"\t"$3"\t"$4"\t"$5,Dl,Ez,mu}' hs_distances.tsv > hs_distances_plus.tsv
