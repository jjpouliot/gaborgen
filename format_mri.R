# Convert structural and functional images to nifti fomat, deoblique structural
# functional should not need to be deobliqued
for (participant_index in 1:length(participants_to_preprocess)) {
  
  current_structural_path <- paste0(participant_directories[participant_index],
                                    '/fMRI/T1_MPRAGE_SAG_P2_ISO_0005')
  
  current_functional_path <- paste0(participant_directories[participant_index],
                                    '/fMRI/BOLD-EPI-CMRR-2S_0006')
  
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
