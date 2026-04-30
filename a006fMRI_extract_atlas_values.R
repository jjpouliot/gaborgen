# data_location <- '/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf/results'
# where_to_copy <- '/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder'
data_location <- '/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/fmri_preprocessed'
where_to_copy <- '/home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/roi_timeseries'

use_unblurred_and_unscaled_func <- T

# copy design matrices; I don't want this to be necessary, it should happen in a006prep_design_matrix.R
# full_design_matrix_paths <- list.files(
#   path = data_location,
#   pattern = 'X.nocensor.xmat.1D$',
#   recursive = T
# )

# full_design_matrix_new_names <- sub(
#   replacement = "_",
#   pattern = "/",
#   x = full_design_matrix_paths
# )

# file.copy(
#   from = paste0(data_location, "/", full_design_matrix_paths),
#   to = paste0(where_to_copy, "/", full_design_matrix_new_names)
# )

# # copy censor info
# censor_info_paths <- list.files(
#   path = data_location,
#   pattern = 'combined_2.1D$',
#   recursive = T
# )

# file.copy(
#   from = paste0(data_location, "/", censor_info_paths),
#   to = paste0(where_to_copy, "/", basename(censor_info_paths))
# )

# copy functional data
# can also look for my pb files in /home/andrewfarkas/Research_data/EEG/Gaborgen24_EEG_fMRI/pb_files
processed_functional_paths <- list.files(
  path = data_location,
  pattern = "pb05",
  recursive = T
)

file.copy(
  from = paste0(data_location, "/", processed_functional_paths),
  to = paste0(where_to_copy, "/", basename(processed_functional_paths))
)

if (use_unblurred_and_unscaled_func) {
  processed_functional_paths <- list.files(
    path = data_location,
    pattern = "pb03",
    recursive = T
  )

  file.copy(
    from = paste0(data_location, "/", processed_functional_paths),
    to = paste0(where_to_copy, "/", basename(processed_functional_paths))
  )
}

# add what is necessary here to make proper function HCP SUIT mask for ROI timeseries extractions
# From emoclips, not so sure about every part of this code. The problem is that the mask may
# extend beyond some ares, and those areas become new ROIs on the bottom of the atlas list.
# The reason I want to try this new atlas is that I think it may better fit the MNI template.
# A manual check suggests it aligns so much better than the previous. I am using HCPex_resam, the same
# one used for emoclips pilots.

# #remaking atlas
# # 1) Align HCPex with NN, then force integer type
# 3dAllineate -base atlas_materials/MNI152_2009_template.nii.gz \
#   -source atlas_materials/HCPex.nii.gz \
#   -prefix atlas_materials/HCPex_aligned_tmp.nii.gz \
#   -final NN -overwrite

# 3dcalc \
#   -a atlas_materials/HCPex_aligned_tmp.nii.gz \
#   -expr 'a' \
#   -datum short \
#   -prefix atlas_materials/HCPex_aligned_2.nii.gz

# # 2) Resample cerebellum segmentation with NN and integer type
# 3dresample -rmode NN \
#   -master atlas_materials/HCPex_aligned.nii.gz \
#   -inset atlas_materials/Cerebellum-MNIsegment.nii.gz \
#   -prefix atlas_materials/Cerebellum-MNIsegment_resampled_tmp.nii.gz

# 3dcalc \
#   -a atlas_materials/Cerebellum-MNIsegment_resampled_tmp.nii.gz \
#   -expr 'a' \
#   -datum short \
#   -prefix atlas_materials/Cerebellum-MNIsegment_resampled_2.nii.gz

# # 3) Merge by replacement (no summed labels)
# 3dcalc -a atlas_materials/HCPex_aligned_2.nii.gz \
#   -b atlas_materials/Cerebellum-MNIsegment_resampled_2.nii.gz \
#   -expr 'a*(1-step(b)) + step(b)*(b+1000)' \
#   -datum short \
#   -prefix atlas_materials/HCPex_SUIT_atlas_2.nii.gz

# @Atlasize -dset atlas_materials/HCPex_aligned_2.nii.gz -lab_file atlas_materials/HCPex_labels.txt 1 0 -atlas_name HCPex_2 -atlas_type G -auto_backup

# @Atlasize -dset atlas_materials/HCPex_SUIT_atlas_2.nii.gz -lab_file atlas_materials/HCPex_SUIT_labels.txt 1 0 -atlas_name HCPex_SUIT_2 -atlas_type G -auto_backup

# 3dmaskdump -noijk -mask atlas_materials/HCPex_SUIT_atlas_2.nii.gz \
#   atlas_materials/HCPex_SUIT_atlas_2.nii.gz | sort -n -u | tail

# cd /home/andrewfarkas/Research_data/multimodal/emoclips/processed

# 3dresample \
#   -master /home/andrewfarkas/Research_data/multimodal/emoclips/processed/EMOCLIPS001.results/pb04.EMOCLIPS001.r01.scale+tlrc \
#   -inset  HCPex_SUIT_atlas_2.nii.gz \
#   -prefix HCPex_resam+tlrc \
#   -rmode NN

# calculate ROI times series per participant, just day 1 right now

participants_to_process <- c(
  "101", # made
  "102", # 22% censored# made
  "103", # made
  "106", # made
  "107", # made
  "108", # made
  "109", # made
  "113", # made
  "114", # 23% censored # made
  "115", # made
  "116", # made
  "117", # made
  "119", # made
  "120", # 16% censored # moved over 2 mm # made
  "121", # made
  "122", # made
  "123", # made
  #"124", # 34% censored moved head to every cue
  "125", # made
  "126", # made
  "127", # made
  "128", # made
  "129", # made
  "131", # made
  "132", # made
  "133", # made
  "134", # made
  "135", # made
  # "136", # 50% censored
  "137", # made
  "138", # made
  # "139", # 31% censored, movement over 3mm
  "140", # 23.7% censored # made
  "141", # made
  # double-check these
  #"142", #probably asleep # made
  "143", #probably asleep, not convinced enough to keep them out # made
  #"144", # 29% censored severe TSNR warnings, large pitch shifts above 5 degrees on way out of acquisition
  #
  "145", # made
  #"146", # messed up alignment
  #"147", #asleep , not aligned right?
  #"148", #messed up alignment
  "149", # made
  "150", # made
  "151", # made
  "152", # made
  "153", # made
  "154", # made
  "155", # made
  "158", # made
  "159", # made
  "160", # made
  "161" # made
)

skip_all_but_v1_bun <- F
if (use_unblurred_and_unscaled_func) {
  for (i in 1:length(participants_to_process)) {
    where_to_copy
    setwd(where_to_copy) # go to ROI folder

    # run this terminal command, expect per participant, gets timeseries per ROI into text file
    #   3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_149.r01.scale+tlrc > 149_roi_stats.txt

    if (participants_to_process[i] < 123) {
      if (!skip_all_but_v1_bun) {
        system2(
          'tcsh',
          args = c(
            '-c',
            shQuote(paste(
              # '3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_',
              '3dROIstats -mask HCPex_resam+tlrc -nzmean -median pb03.GABORGEN24_',
              participants_to_process[i],
              ".r01.volreg+tlrc > ",
              participants_to_process[i],
              "_roi_stats_pb03.txt",
              sep = ""
            ))
          )
        )
      }
      # V1 fovea L HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_L.nii.gz -nzmean -median pb03.GABORGEN24_',
            participants_to_process[i],
            ".r01.volreg+tlrc > ",
            participants_to_process[i],
            "_roi_stats_pb03_HCPex_resam_V1_fovea_L.txt",
            sep = ""
          ))
        )
      )
      # V1 fovea R HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_R.nii.gz -nzmean -median pb03.GABORGEN24_',
            participants_to_process[i],
            ".r01.volreg+tlrc > ",
            participants_to_process[i],
            "_roi_stats_pb03_HCPex_resam_V1_fovea_R.txt",
            sep = ""
          ))
        )
      )

      # # V1_julich_fovea_L.nii.gz
      #     system2(
      #   'tcsh',
      #   args = c(
      #     '-c',
      #     shQuote(paste(
      #       '3dROIstats -mask V1_julich_fovea_L.nii.gz -nzmean -median pb05.GABORGEN24_',
      #       participants_to_process[i],
      #       ".r01.scale+tlrc > ",
      #       participants_to_process[i],
      #       "_roi_stats_V1_julich_fovea_L.txt",
      #       sep = ""
      #     ))
      #   )
      # )

      # # V1_julich_fovea_R.nii.gz
      #     system2(
      #   'tcsh',
      #   args = c(
      #     '-c',
      #     shQuote(paste(
      #       '3dROIstats -mask V1_julich_fovea_R.nii.gz -nzmean -median pb05.GABORGEN24_',
      #       participants_to_process[i],
      #       ".r01.scale+tlrc > ",
      #       participants_to_process[i],
      #       "_roi_stats_V1_julich_fovea_R.txt",
      #       sep = ""
      #     ))
      #   )
      # )
    } else {
      if (!skip_all_but_v1_bun) {
        system2(
          'tcsh',
          args = c(
            '-c',
            shQuote(paste(
              # '3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_',
              '3dROIstats -mask HCPex_resam+tlrc -nzmean -median pb03.GABORGEN24_DAY1_',
              participants_to_process[i],
              ".r01.volreg+tlrc > ",
              participants_to_process[i],
              "_roi_stats_pb03.txt",
              sep = ""
            ))
          )
        )
      }

      # V1 fovea L HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_L.nii.gz -nzmean -median pb03.GABORGEN24_DAY1_',
            participants_to_process[i],
            ".r01.volreg+tlrc > ",
            participants_to_process[i],
            "_roi_stats_pb03_HCPex_resam_V1_fovea_L.txt",
            sep = ""
          ))
        )
      )
      # V1 fovea R HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_R.nii.gz -nzmean -median pb03.GABORGEN24_DAY1_',
            participants_to_process[i],
            ".r01.volreg+tlrc > ",
            participants_to_process[i],
            "_roi_stats_pb03_HCPex_resam_V1_fovea_R.txt",
            sep = ""
          ))
        )
      )
    }
  }
} else {
  for (i in 1:length(participants_to_process)) {
    where_to_copy
    setwd(where_to_copy) # go to ROI folder

    # run this terminal command, expect per participant, gets timeseries per ROI into text file
    #   3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_149.r01.scale+tlrc > 149_roi_stats.txt

    if (participants_to_process[i] < 123) {
      if (!skip_all_but_v1_bun) {
        system2(
          'tcsh',
          args = c(
            '-c',
            shQuote(paste(
              # '3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_',
              '3dROIstats -mask HCPex_resam+tlrc -nzmean -median pb05.GABORGEN24_',
              participants_to_process[i],
              ".r01.scale+tlrc > ",
              participants_to_process[i],
              "_roi_stats.txt",
              sep = ""
            ))
          )
        )
      }
      # V1 fovea L HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_L.nii.gz -nzmean -median pb05.GABORGEN24_',
            participants_to_process[i],
            ".r01.scale+tlrc > ",
            participants_to_process[i],
            "_roi_stats_HCPex_resam_V1_fovea_L.txt",
            sep = ""
          ))
        )
      )
      # V1 fovea R HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_R.nii.gz -nzmean -median pb05.GABORGEN24_',
            participants_to_process[i],
            ".r01.scale+tlrc > ",
            participants_to_process[i],
            "_roi_stats_HCPex_resam_V1_fovea_R.txt",
            sep = ""
          ))
        )
      )

      # # V1_julich_fovea_L.nii.gz
      #     system2(
      #   'tcsh',
      #   args = c(
      #     '-c',
      #     shQuote(paste(
      #       '3dROIstats -mask V1_julich_fovea_L.nii.gz -nzmean -median pb05.GABORGEN24_',
      #       participants_to_process[i],
      #       ".r01.scale+tlrc > ",
      #       participants_to_process[i],
      #       "_roi_stats_V1_julich_fovea_L.txt",
      #       sep = ""
      #     ))
      #   )
      # )

      # # V1_julich_fovea_R.nii.gz
      #     system2(
      #   'tcsh',
      #   args = c(
      #     '-c',
      #     shQuote(paste(
      #       '3dROIstats -mask V1_julich_fovea_R.nii.gz -nzmean -median pb05.GABORGEN24_',
      #       participants_to_process[i],
      #       ".r01.scale+tlrc > ",
      #       participants_to_process[i],
      #       "_roi_stats_V1_julich_fovea_R.txt",
      #       sep = ""
      #     ))
      #   )
      # )
    } else {
      if (!skip_all_but_v1_bun) {
        system2(
          'tcsh',
          args = c(
            '-c',
            shQuote(paste(
              # '3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_',
              '3dROIstats -mask HCPex_resam+tlrc -nzmean -median pb05.GABORGEN24_DAY1_',
              participants_to_process[i],
              ".r01.scale+tlrc > ",
              participants_to_process[i],
              "_roi_stats.txt",
              sep = ""
            ))
          )
        )
      }

      # V1 fovea L HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_L.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_',
            participants_to_process[i],
            ".r01.scale+tlrc > ",
            participants_to_process[i],
            "_roi_stats_HCPex_resam_V1_fovea_L.txt",
            sep = ""
          ))
        )
      )
      # V1 fovea R HCPex
      system2(
        'tcsh',
        args = c(
          '-c',
          shQuote(paste(
            '3dROIstats -mask HCPex_resam_V1_fovea_R.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_',
            participants_to_process[i],
            ".r01.scale+tlrc > ",
            participants_to_process[i],
            "_roi_stats_HCPex_resam_V1_fovea_R.txt",
            sep = ""
          ))
        )
      )
    }
  }
}

# Make design matrices with shock, BLOCK

participants_to_find <- c(159:161)

# data_directory <- "~/research_data/gaborgen/raw_data"
data_directory <- "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/raw_data"

log_file_paths <- list.files(
  path = data_directory,
  pattern = "logfile.dat$",
  recursive = T,
  full.names = T
)

# where_to_save <- "/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder"
where_to_save <- "/home/andrewfarkas/tmp/restore3/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI/roi_data_and_info/shock_design_matrix_column_to_add/"

# write shock times

for (i in 1:length(log_file_paths)) {
  log_file <- read.delim(log_file_paths[i], header = T, sep = ",")

  name_to_save <- sub(
    basename(log_file_paths[i]),
    pattern = "logfile.dat$",
    replacement = "shock_onsets.1D"
  )

  # approximate shock onset,
  shock_onsets <- subset(
    log_file,
    paired == 1,
    select = timeSinceFirstTR,
    drop = T
  ) +
    2 # shock comes two seconds after the paired stimulus (approximately, double check needed)

  write(
    x = shock_onsets,
    file = paste0(where_to_save, name_to_save),
    ncolumns = 1
  )
}

# below i don't think is necessary because of other 006 script
#BLOCK(d, p) ; d = duration seconds, p = amplitude, 1 is easiest to understand and should be default unless you know what you're doing

shock_onset_stim_files <- list.files(
  path = where_to_save,
  pattern = "shock_onsets.1D$",
  recursive = T,
  full.names = T
)

setwd(where_to_save)


#enters(\n) will cause command to enter early

for (i in 1:length(shock_onset_stim_files)) {
  sub_name <- sub(
    basename(shock_onset_stim_files[i]),
    pattern = "shock_onsets.1D$",
    replacement = ""
  )

  system2(
    'tcsh',
    args = c(
      '-c',
      shQuote(paste(
        '3dDeconvolve -nodata 1070 2.0 -num_stimts 1 -stim_times 1',
        paste0(shock_onset_stim_files[i], " "),
        '"BLOCK(.1,1)" -x1D',
        paste0(sub_name, "design_matrix.xmat.1D "),
        '-xjpeg',
        paste0(sub_name, "design_matrix.jpg "),
        '-x1D_stop'
      ))
    )
  )
}

# 3dDeconvolve -nodata 1070 2.0 \
# -num_stimts 1 \
# -stim_times 1 /Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder/gaborgen24_fMRI_Day1_121_shock_onsets.1D "BLOCK(.1,1)" \
# -x1D design_matrix.xmat.1D \
# -xjpeg design_matrix.jpg \
# -x1D_stop
#
#
#
#
#
# 3dDeconvolve -nodata 200 2.0 \
# -num_stimts 1 \
# -stim_times 1 stimfile.1D 'BLOCK(2,1)' \
# -x1D design_matrix.xmat.1D \
# -xjpeg design_matrix.jpg \
# -x1D_stop
