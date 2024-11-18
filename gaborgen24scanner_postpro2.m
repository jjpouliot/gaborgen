%% Script for analyzing all gaborgen24 data after initial preprocessing
% this script opens the folders with EEGlab data, one for each condition
% then it runs the sliding window code and saves the condition trialwise
% amplitide at 15 hz or whatever other frequency
% Get a list of all files and folders in the current directory
clear 

files = dir("gaborgen*");

% Filter out the non-folder entries
dirFlags = [files.isdir];

% Extract the names of the folders
folderNames = {files(dirFlags).name};

% Remove the '.' and '..' folders
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

% Display the folder names
disp('Folders in the current working directory:');
disp(folderNames);

% loop over subjects
for subindex = 6:size(folderNames,2)

    eval(['cd ' folderNames{subindex} '/EEG'])

 amp3dwin = [];

     conditioncounter = 1; 

    for condition = 21:24

        datacon = []; 

        basename = folderNames{subindex};
        currentFile = sprintf('%s_%d_EEG.mat.slidwin.mat', basename(end-2:end), condition);

        load(currentFile);


        [amp3dwin(:,:,conditioncounter), freqs3dwin, fftcompwin] = freqtag_FFT3D(outmat.winmat, 600);
        
       
 conditioncounter = conditioncounter+1; 

    end % condition

    outname = [basename(end-2:end) 'spectrum.mat']

    save(outname,'amp3dwin','-mat')


    cd ..

    pause(.10)

    cd ..

    pause(.5)

    fclose('all');
    close('all');

end % subject
