#!/usr/bin/env bash

BVp="/hpc/soft/brainvisa/brainvisa/bin/"
FSLPREFIX="fsl5.0-"
mkdir surfaces
OUTp="surfaces/"

IN="sub-Musly_ses-01_T1w_Mean_denoised_debiased_in-inia19_grey_white"
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