#!/bin/tcsh

# This script aligns the Human Connectome Project extended atlas (Huang et al., 2021) with the MNI 2009c asymmetric template, then combines it with the spatially unbiased atlas template of the cerebellum and brainstem (Diedrichsen, 2006). It will create the custom atlases in the atlas materials folder and make a CustomAtlases directory if not already created.

3dAllineate -base atlas_materials/MNI152_2009_template.nii.gz -source atlas_materials/HCPex.nii.gz -prefix atlas_materials/HCPex_aligned.nii.gz -final NN -overwrite

mkdir ~/CustomAtlases

@AfniEnv -set AFNI_SUPP_ATLAS_DIR ~/CustomAtlases

3dresample -master atlas_materials/HCPex_aligned.nii.gz -inset atlas_materials/Cerebellum-MNIsegment.nii.gz -prefix atlas_materials/Cerebellum-MNIsegment_resampled.nii.gz

3dcalc -a atlas_materials/HCPex_aligned.nii.gz -b atlas_materials/Cerebellum-MNIsegment_resampled.nii.gz -expr 'a + ispositive(b)*(b+1000)' -prefix atlas_materials/HCPex_SUIT_atlas.nii.gz

@Atlasize -dset atlas_materials/HCPex_aligned.nii.gz -lab_file atlas_materials/HCPex_labels.txt 1 0 -atlas_name HCPex -atlas_type G -auto_backup

@Atlasize -dset atlas_materials/HCPex_SUIT_atlas.nii.gz -lab_file atlas_materials/HCPex_SUIT_labels.txt 1 0 -atlas_name HCPex_SUIT -atlas_type G -auto_backup