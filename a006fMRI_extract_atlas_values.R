data_location <- '/Users/andrewfarkas/UFL Dropbox/Andrew Farkas/TwoDayFearConditioning/AFNI_results_andrewf/results'
where_to_copy <- '/Users/andrewfarkas/research_data/gaborgen/temp_analyses/ROI_folder'

# copy design matrices
full_design_matrix_paths <- list.files(path = data_location, pattern = 'X.nocensor.xmat.1D$', recursive = T)

full_design_matrix_new_names <- sub(replacement = "_", pattern = "/", x = full_design_matrix_paths)

file.copy(from = paste0(data_location, "/", full_design_matrix_paths), to = paste0(where_to_copy, "/", full_design_matrix_new_names))

# copy censor info
censor_info_paths <- list.files(path = data_location, pattern = 'combined_2.1D$', recursive = T)

file.copy(from = paste0(data_location, "/", censor_info_paths), to = paste0(where_to_copy, "/", basename(censor_info_paths)))


# copy functional data

processed_functional_paths <- list.files(path = data_location, pattern = "pb05", recursive = T)

file.copy(from = paste0(data_location, "/", processed_functional_paths), to = paste0(where_to_copy, "/", basename(processed_functional_paths)))


# calculate ROI times series per participant, just day 1 right now

participants_to_process <- c(101:154)

for (i in 1:length(participants_to_process)) {
  setwd(where_to_copy) # go to ROI folder
  
  # run this terminal command, expect per participant, gets timeseries per ROI into text file
  #   3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_149.r01.scale+tlrc > 149_roi_stats.txt
  
  if (participants_to_process[i] < 123) {
    system2('tcsh', 
            args = c('-c', shQuote(paste('3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_', 
                                         participants_to_process[i], ".r01.scale+tlrc > ", participants_to_process[i], "_roi_stats.txt",sep = ""))))
    
  } else {
    system2('tcsh', 
             args = c('-c', shQuote(paste('3dROIstats -mask HCPex_SUIT_func_mask.nii.gz -nzmean -median pb05.GABORGEN24_DAY1_', 
                                    participants_to_process[i], ".r01.scale+tlrc > ", participants_to_process[i], "_roi_stats.txt",sep = ""))))
  }
}


