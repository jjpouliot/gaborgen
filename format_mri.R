# Convert files to NIFTI format and deobliques
# If there are multiple of the same scan, it processes the last one recorded by
# looking for the largest number at the end of file
# legacy version did this by directory, now everything in the directory gets converted to dicoms first
for (directory_index in 1:length(participant_directories)) {
  print(paste(
    "Dicom 2 nifti MRI conversion for participant ",
    basename(participant_directories[directory_index])
  ))

  setwd(paste0(participant_directories[directory_index], "/fMRI"))

  # This terminal command dcm2niix needs to be installed on your system, homebrew is
  # a good option on a mac
  system2(
    'tcsh',
    args = c(
      '-c',
      shQuote(paste(
        'dcm2niix -w 0',
        paste0(participant_directories[directory_index], "/fMRI")
      ))
    )
  )

  current_participant_found_structural_volume <- list.files(
    paste0(participant_directories[directory_index], "/fMRI"),
    pattern = 't1_mprage_sag_p2_iso.*.nii',
    full.names = T
  )

  # assumes the last recorded volume is the correct one if there are duplicates
  current_structural_volume <-
    current_participant_found_structural_volume[
      which.max(as.numeric(
        sub(
          ".*_(\\d+)\\.nii$",
          "\\1",
          current_participant_found_structural_volume
        )
      ))
    ]

  current_participant_found_functional_volumes <- list.files(
    paste0(participant_directories[directory_index], "/fMRI"),
    pattern = "EPI-cmrr-2s_.*.nii",
    full.names = TRUE
  )

  # Filter out those containing "MoCo", this is a motion correction version that we don't need becauase afni will do motion correction for us
  current_participant_found_functional_volumes <- current_participant_found_functional_volumes[
    !grepl("MoCo", current_participant_found_functional_volumes)
  ]

  current_functional_volumes <-
    current_participant_found_functional_volumes[
      which.max(as.numeric(
        sub(
          ".*_(\\d+)\\.nii$",
          "\\1",
          current_participant_found_functional_volumes
        )
      ))
    ]

  ## For participants where the last function is not the correct one
  if (
    basename(participant_directories[directory_index]) == "GABORGEN24_DAY1_146"
  ) {
    current_functional_volumes <- paste0(
      participant_directories[directory_index],
      "/fMRI/fMRI_BOLD-EPI-cmrr-2s_20250224110209_10.nii"
    )
  }

  current_participant_found_blip_forward <- list.files(
    paste0(participant_directories[directory_index], '/fMRI'),
    pattern = 'fMRI_cmmr-DistMap_AP.*.nii',
    full.names = T
  )

  current_participant_found_blip_forward <- current_participant_found_blip_forward[
    !grepl("MoCo", current_participant_found_blip_forward)
  ]

  current_blip_forward <- current_participant_found_blip_forward[
    which.max(as.numeric(
      sub(".*_(\\d+)\\.nii$", "\\1", current_participant_found_blip_forward)
    ))
  ]

  current_participant_found_blip_reverse <- list.files(
    paste0(participant_directories[directory_index], '/fMRI'),
    pattern = 'fMRI_cmmr-DistMap_PA.*.nii',
    full.names = T
  )

  current_participant_found_blip_reverse <- current_participant_found_blip_reverse[
    !grepl("MoCo", current_participant_found_blip_reverse)
  ]

  current_blip_reverse <- current_participant_found_blip_reverse[
    which.max(as.numeric(
      sub(".*_(\\d+)\\.nii$", "\\1", current_participant_found_blip_reverse)
    ))
  ]

  # Structural deobliqued
  if (
    length(list.files(
      paste0(participant_directories[directory_index], "/fMRI"),
      pattern = '^ANAT_DEOB'
    )) >
      0
  ) {
    print(
      'Deobliqued structural file already appears to exist, skipping conversion'
    )
  } else {
    system2(
      'tcsh',
      args = c(
        '-c',
        shQuote(paste(
          '3dWarp -deoblique -prefix ANAT_DEOB',
          current_structural_volume
        ))
      )
    )
  }

  # Functional deobliqued
  if (
    length(list.files(
      paste0(participant_directories[directory_index], "/fMRI"),
      pattern = '^FUNC_DEOB'
    )) >
      0
  ) {
    print(
      'Deobliqued functional file already appears to exist, skipping conversion'
    )
  } else {
    system2(
      'tcsh',
      args = c(
        '-c',
        shQuote(paste(
          '3dWarp -deoblique -prefix FUNC_DEOB',
          current_functional_volumes
        ))
      )
    )
  }

  # Blip forward deobliqued
  if (
    length(list.files(
      paste0(participant_directories[directory_index], "/fMRI"),
      pattern = '^BLIP_FOR_DEOB'
    )) >
      0
  ) {
    print(
      'Deobliqued blip forward file already appears to exist, skipping conversion'
    )
  } else {
    system2(
      'tcsh',
      args = c(
        '-c',
        shQuote(paste(
          '3dWarp -deoblique -prefix BLIP_FOR_DEOB',
          current_blip_forward
        ))
      )
    )
  }

  # Blip reverse deobliqued
  if (
    length(list.files(
      paste0(participant_directories[directory_index], "/fMRI"),
      pattern = '^BLIP_REV_DEOB'
    )) >
      0
  ) {
    print(
      'Deobliqued blip forward file already appears to exist, skipping conversion'
    )
  } else {
    system2(
      'tcsh',
      args = c(
        '-c',
        shQuote(paste(
          '3dWarp -deoblique -prefix BLIP_REV_DEOB',
          current_blip_reverse
        ))
      )
    )
  }
  setwd(local_git_directory)
}
