#!/usr/bin/env bash

BVp="/hpc/soft/brainvisa/brainvisa/bin/"
FSLPREFIX="fsl5.0-"
mkdir surfaces
OUTp="surfaces/"

IN="In-sub-Maga_ses-01_mp2rageT1w_denoised_debiased_in-inia19_sub-Maga_res-7WM_4CSF_0mvt_macaque_vs_all_t"
cIN=$OUTp$IN

"${FSLPREFIX}fslmaths" $IN $cIN -odt short
"${FSLPREFIX}fslmaths" $cIN -mul -1 $cIN -odt short
"${FSLPREFIX}fslmaths" $cIN -thr 5 -bin $cIN -odt short

"${BVp}AimsMesh" -i $cIN
SUFF="_1_0.gii" # suffix of the surface that you want to smooth
sIN=$cIN$SUFF
sOUT="${OUTp}smoothed_${IN}${SUFF}"
"${BVp}AimsMeshSmoothing" -i $sIN --algoType laplacian -o $sOUT

# rm "${OUTp}*.minf"