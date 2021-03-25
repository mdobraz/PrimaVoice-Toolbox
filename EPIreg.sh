#!/usr/bin/env bash

export ANTSPATH=${ANTSPATH:="/hpc/soft/ANTS/antsbin/bin/"}

# Registration to T1w seems to work a bit better
anat=/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Maga/anat/segmentation/inia19_sub-Maga_ses-01_T1w0p6mmDenoised_cropped_debiased/sub-Maga_ses-01_T1w0p6mmDenoised_cropped_debiased-rigid-in-template.nii.gz
anat_brain=/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Maga/anat/segmentation/inia19_sub-Maga_ses-01_T1w0p6mmDenoised_cropped_debiased/sub-Maga_ses-01_T1w0p6mmDenoised_cropped_debiased-brain-rigid-in-template.nii.gz
anat_brainmask=/hpc/banco/Primavoice_Data_and_Analysis/analysis_sub-Maga/anat/segmentation/inia19_sub-Maga_ses-01_T1w0p6mmDenoised_cropped_debiased/sub-Maga_ses-01_T1w0p6mmDenoised_cropped_debiased-brainmask-in-template.nii.gz

in=/hpc/banco/Primavoice_Data_and_Analysis/sub-Maga/ses-13/func/sub-Maga_ses-13_task-sparse_run-04_bold.nii
ref_scan=ref_scan

BETf=0.3
N4ITER=6

# average the time series
${ANTSPATH}antsMotionCorr  -d 3 -a $in -o ${ref_scan}.nii.gz 

# Debias with N4
ref_scan_N4=${ref_scan}_N4
cp ${ref_scan}.nii.gz ${ref_scan_N4}.nii.gz

for ((i = 1 ; i  <= $N4ITER ; i++)); do
  ${ANTSPATH}N4BiasFieldCorrection -i ${ref_scan_N4}.nii.gz -o ${ref_scan_N4}.nii.gz
done

# BET
bash /hpc/banco/Primavoice_Scripts/Pipeline3/functions/bash/T1xT2BET.sh \
-t1 $ref_scan_N4 -t2 $ref_scan_N4 -n 3 -p fsl5.0- -f $BETf

# IterREGBET
ref2anat=ref2anat
bash /hpc/banco/Primavoice_Scripts/Pipeline3/functions/bash/IterREGBET.sh \
-inw $ref_scan_N4 -inb ${ref_scan_N4}_BET -refb $anat_brain \
-n 3 -dof 6 -p fsl5.0- -xp $ref2anat

# Move anat to ref space
anat2ref=anat2ref # useless because already done in IterREGBET
fsl5.0-convert_xfm -omat ${anat2ref}.xfm -inverse ${ref2anat}.xfm

anat_in_ref=anat_in_ref.nii.gz
anat_brain_in_ref=anat_brain_in_ref.nii.gz
fsl5.0-flirt -in $anat -ref $ref_scan_N4 -out $anat_in_ref -applyxfm -init ${anat2ref}.xfm
fsl5.0-flirt -in $anat_brain -ref $ref_scan_N4 -out $anat_brain_in_ref -applyxfm -init ${anat2ref}.xfm

# SyN
/hpc/soft/ANTS/ANTs/Scripts/antsRegistrationSyNQuick.sh -d 3 \
-f $anat_brain_in_ref -f $anat_in_ref \
-m ${ref_scan_N4}_IRbrain.nii.gz -m ${ref_scan_N4}.nii.gz \
-o ${ref_scan_N4}_ -j 1

# Apply transforms to original ref_scan
${ANTSPATH}/antsApplyTransforms -i ${ref_scan}.nii.gz \
-r ${ref_scan_N4}_Warped.nii.gz \
-o ${ref_scan}_Warped.nii.gz \
-t ${ref_scan_N4}_1Warp.nii.gz -t ${ref_scan_N4}_0GenericAffine.mat




# Apply to other runs
in=/hpc/banco/Primavoice_Data_and_Analysis/sub-Maga/ses-14/func/sub-Maga_ses-14_task-sparse_run-01_bold.nii
run_avg=run_avg
# average the time series
${ANTSPATH}antsMotionCorr  -d 3 -a $in -o ${run_avg}.nii.gz

# SyN
run_avg_to_ref=${run_avg}_to_ref_
/hpc/soft/ANTS/ANTs/Scripts/antsRegistrationSyNQuick.sh -d 3 \
-f ${ref_scan}_Warped.nii.gz \
-m ${run_avg}.nii.gz \
-o $run_avg_to_ref -j 1 # -x ${ref_scan_N4}_Warped.nii.gz
# ----> does not give good enough results




# Other strategy woth additional steps (N4, brainmask,SyN to anat)

# N4
run_avg_N4=${run_avg}_N4
cp ${run_avg}.nii.gz ${run_avg_N4}.nii.gz
for ((i = 1 ; i  <= $N4ITER ; i++)); do
  ${ANTSPATH}N4BiasFieldCorrection -i ${run_avg_N4}.nii.gz -o ${run_avg_N4}.nii.gz
done

# Get brain mask
${ANTSPATH}/antsApplyTransforms -i $anat_brain_in_ref \
-r ${run_avg}.nii.gz \
-o ${run_avg}_brainmask.nii.gz \
-t ${run_avg_to_ref}1InverseWarp.nii.gz -t [${run_avg_to_ref}0GenericAffine.mat,1] \
-n NearestNeighbor
fsl5.0-fslmaths ${run_avg}_brainmask.nii.gz -bin ${run_avg}_brainmask.nii.gz
fsl5.0-fslmaths $run_avg_N4 -mas ${run_avg}_brainmask.nii.gz ${run_avg_N4}_BET

# IterREGBET
avg2anat=avg2anat
bash /hpc/banco/Primavoice_Scripts/Pipeline3/functions/bash/IterREGBET.sh \
-inw $run_avg_N4 -inb ${run_avg_N4}_BET -refb $anat_brain \
-n 3 -dof 12 -p fsl5.0- -xp $avg2anat

# SyN
/hpc/soft/ANTS/ANTs/Scripts/antsRegistrationSyNQuick.sh -d 3 \
-f $anat_brain_in_ref -f $anat_in_ref \
-m ${run_avg_N4}_IRbrain.nii.gz -m ${run_avg_N4}.nii.gz \
-o ${run_avg_N4}_ -j 1

# Apply transforms to original run_avg
${ANTSPATH}/antsApplyTransforms -i ${run_avg}.nii.gz \
-r ${run_avg_N4}_Warped.nii.gz \
-o ${run_avg}_Warped.nii.gz \
-t ${run_avg_N4}_1Warp.nii.gz -t ${run_avg_N4}_0GenericAffine.mat

# ----->> not better than direct SyN (even worse)





##### Trying to SyN 2 fieldmap corrrected volumes:

in1=/hpc/banco/Primavoice_Data_and_Analysis/tests_segmentation/realign_runs/sub-Apache_ses-03_task-sparse_run-01_bold_198.nii.gz
in2=/hpc/banco/Primavoice_Data_and_Analysis/tests_segmentation/realign_runs/sub-Apache_ses-04_task-sparse_run-01_bold_381.nii.gz

/hpc/soft/ANTS/ANTs/Scripts/antsRegistrationSyNQuick.sh -d 3 \
-f $in1 \
-m $in2 \
-o run4_ -j 1 


fsl5.0-flirt -in $in2 -ref $in1 -out flirt6 -dof 6 -cost normmi


# ----->> One solution might be to only replace the flirt stage by a SyN stage in the script

# Metric 0 - MeanSquareDifference, 1 - Cross-Correlation, 2-Mutual Information , 3-SMI 
${ANTSPATH}MeasureImageSimilarity 3 1 $in1 run4_Warped.nii.gz




run_avg_aff=run_avg_aff
/hpc/soft/ANTS/ANTs/Scripts/antsRegistrationSyNQuick.sh -d 3 \
-f ${ref_scan}.nii.gz \
-m ${run_avg}.nii.gz \
-o ${run_avg_aff}_ -t a

run_avg_SyN=run_avg_SyN
${ANTSPATH}/antsApplyTransforms -i ${run_avg}.nii.gz \
-r ${ref_scan_N4}_Warped.nii.gz \
-o ${run_avg_SyN}_Warped.nii.gz \
-t ${ref_scan_N4}_1Warp.nii.gz -t ${run_avg_aff}_0GenericAffine.mat














# Motion correction 
# in: /hpc/banco/Primavoice_Data_and_Analysis/tests_segmentation/realign_runs/

export ANTSPATH=${ANTSPATH:="/opt/ANTs/bin/"}

in2=sub-Apache_ses-02_task-sparse_run-03_roi.nii.gz
in5=sub-Apache_ses-05_task-sparse_run-04_roi.nii.gz
in12=sub-Apache_ses-12_task-sparse_run-03_roi.nii.gz

avg2=avg-02.nii.gz
avg5=avg-05.nii.gz
avg12=avg-12.nii.gz

out2=ses-02
out5=ses-05
out12=ses-12

fixed=Reference_scan.nii.gz

in=$in12
out=$out12

${ANTSPATH}antsMotionCorr  -d 3 -a $in2 -o $avg2
${ANTSPATH}antsMotionCorr  -d 3 -a $in5 -o $avg5
${ANTSPATH}antsMotionCorr  -d 3 -a $in12 -o $avg12

${ANTSPATH}antsMotionCorr  -d 3 -o [${out},${out}.nii.gz,${out}_avg.nii.gz] \
-m MI[ $fixed , $in , 1 , 1 , Random, 0.05  ] \
-t Rigid[ 0.005 ] -i 20 -u 1 -e 1 -s 0 -f 1  \
-m MI[ $fixed , $in , 1 , 1 , Random, 0.05  ] \
-t Affine[ 0.005 ] -i 20 -u 1 -e 1 -s 0 -f 1 \
-m MI[  $fixed , $in , 1 , 2 ] \
-t SyN[0.15,3,0.5] -i 20 -u 1 -e 1 -s 0 -f 1 -n 10







${ANTSPATH}antsMotionCorr  -d 3 -o [${out},${out}.nii.gz,${out}_avg.nii.gz] \
-m MI[ $fixed , $in , 1 , 1 , Random, 0.05  ] \
-t Affine[ 0.005 ] -i 20 -u 1 -e 1 -s 0 -f 1 -n 10






bash antsRegistrationSyNQuick.sh -d 3 \
-f $avg12 \
-m $avg5 \
-o SyN_ -j 1












# Registering two func volumes
out=avg_sub-Maga_ses-13_task-sparse_run-04_bold

${ANTSPATH}antsMotionCorr  -d 3 -a $in -o ${out}.nii.gz # average the time series

/hpc/soft/ANTS/ANTs/Scripts/antsRegistrationSyNQuick.sh -d 3 \
-f $fixed -m ${out}.nii.gz \
-o ${out}_SyN_ -j 1 -t so # so: deformable syn only






# c3d_affine_tool tests

${ANTSPATH}/antsApplyTransforms -i $in \
-r $ref \
-o u.nii.gz \
-t ants.mat \
-t fsl.mat

${ANTSPATH}/antsApplyTransforms -i $in \
-r $ref \
-o v.nii.gz \
-t both.mat

${ANTSPATH}/antsApplyTransforms -i $ref \
-r $in \
-o x.nii.gz \
-t both_inv.mat