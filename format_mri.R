# Convert structural and functional images to nifti fomat, deoblique structural
# functional should not need to be deobliqued
# If there are multiple of the same scan, it processes the last one recorded by
# looking for the largest number at the end of folder
for (participant_index in 1:length(participants_to_preprocess)) {
  
  print(paste("Reformat MRI volumnes for participant ", 
              participants_to_preprocess[participant_index]))
  
  current_participant_found_structural_directories <- list.files(
    paste0(participant_directories[participant_index],'/fMRI'),
    pattern = 'T1_MPRAGE_SAG_P2_ISO', full.names = T)
  
  current_structural_path <- 
    current_participant_found_structural_directories[
      which.max(as.numeric(
        sub(".*?(\\d+)$", "\\1", 
            current_participant_found_structural_directories)))]
  
  current_participant_found_functional_directories <- list.files(
    paste0(participant_directories[participant_index],'/fMRI'),
    pattern = 'BOLD-EPI-CMRR-2S', full.names = T)
  
  current_functional_path <- current_participant_found_functional_directories[
    which.max(as.numeric(
      sub(".*?(\\d+)$", "\\1", 
          current_participant_found_functional_directories)))]
  
  current_participant_found_blip_forward <- list.files(
    paste0(participant_directories[participant_index],'/fMRI'),
    pattern = 'CMMR-DISTMAP_AP', full.names = T)
  
  current_blip_forward_path <- current_participant_found_blip_forward[
    which.max(as.numeric(
      sub(".*?(\\d+)$", "\\1", 
          current_participant_found_blip_forward)))]
  
  current_participant_found_blip_reverse <- list.files(
    paste0(participant_directories[participant_index],'/fMRI'),
    pattern = 'CMMR-DISTMAP_PA', full.names = T)
  
  current_blip_reverse_path <- current_participant_found_blip_reverse[
    which.max(as.numeric(
      sub(".*?(\\d+)$", "\\1", 
          current_participant_found_blip_reverse)))]
  
  
  # This terminal command (dcm2niix) may need to be installed on your system, homebrew is
  # a good option on a mac
  
  # Structural to NIFTI
  if (length(list.files(current_structural_path, pattern = 'nii$')) == 1 &
      length(list.files(current_structural_path, pattern = 'json$')) == 1) {
    
    print('Structural NIFTI file already appears to exist, skipping conversion')
    
  } else if (length(list.files(current_structural_path, pattern = 'nii$')) > 1 &
             length(list.files(current_structural_path, pattern = 'json$')) > 1) {
    
    stop('Multiple NIFTI files found in the current structural path. There should be 1 or 0.')
    
  } else {
    system2('tcsh', 
            args = c('-c', shQuote(paste('dcm2niix', current_structural_path))))
  }
  
  
  # Functional to NIFTI
  if (length(list.files(current_functional_path, pattern = 'nii$')) == 1 &
      length(list.files(current_functional_path, pattern = 'json$')) == 1) {
    
    print('Functional NIFTI file already appears to exist, skipping conversion')
    
  } else if (length(list.files(current_functional_path, pattern = 'nii$')) > 1 &
             length(list.files(current_functional_path, pattern = 'json$')) > 1) {
    
    stop('Multiple NIFTI files found in the current functional path. There should be 1 or 0.')
    
  }  else {
    system2('tcsh', 
            args = c('-c', shQuote(paste('dcm2niix', current_functional_path))))
  }

  # Blip distortion maps to NIFTI
  if (length(list.files(current_blip_forward_path, pattern = 'nii$')) == 1 &
      length(list.files(current_blip_forward_path, pattern = 'json$')) == 1) {
    print('Blip forward NIFTI file already appears to exist, skipping conversion')
  } else if (length(list.files(current_blip_forward_path, pattern = 'nii$')) > 1 &
             length(list.files(current_blip_forward_path, pattern = 'json$')) > 1) {
    stop('Multiple NIFTI files found in the current blip forward path. There should be 1 or 0.')
  }  else {
    system2('tcsh', 
            args = c('-c', shQuote(paste('dcm2niix', current_blip_forward_path))))
  }

  if (length(list.files(current_blip_reverse_path, pattern = 'nii$')) == 1 &
      length(list.files(current_blip_reverse_path, pattern = 'json$')) == 1) {
    print('Blip reverse NIFTI file already appears to exist, skipping conversion')
  } else if (length(list.files(current_blip_reverse_path, pattern = 'nii$')) > 1 &
             length(list.files(current_blip_reverse_path, pattern = 'json$')) > 1) {
    stop('Multiple NIFTI files found in the current blip reverse path. There should be 1 or 0.')
  }  else {
    system2('tcsh', 
            args = c('-c', shQuote(paste('dcm2niix', current_blip_reverse_path))))
  }
  
  
  # Structural deobliqued
  if (length(list.files(current_structural_path, pattern = '^ANAT_DEOB')) > 0) {
    print('Deobliqued structural file already appears to exist, skipping conversion')
  } else {
    
    # necessary so the obliqued data is written to the right place
    setwd(current_structural_path)
    
    current_structural_NIFTI_path <- list.files(current_structural_path, 
                                                pattern = '.nii$', 
                                                full.names = T)
    system2('tcsh', 
            args = c('-c', 
                     shQuote(paste('3dWarp -deoblique -prefix ANAT_DEOB', 
                                   current_structural_NIFTI_path))))
    
    setwd(local_git_directory)
  }
}
