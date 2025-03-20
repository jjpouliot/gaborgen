%% Script for analyzing all gaborgen24 data after initial preprocessing
% this script opens the folders with EEGlab data, one for each condition
% then it runs the sliding window code and saves the condition trialwise
% amplitide at 15 hz or whatever other frequency
% Get a list of all files and folders in the current directory
clear 
temp99 = eeglab; 

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
for subindex = 1:size(folderNames,2)

    eval(['cd ' folderNames{subindex} '/EEG'])

    ssVEP_ampmat = []; 

    all_amplitudes = [];  % To store the amplitudes across conditions %% intentar
    
    for condition = 21:24

        conditioncounter = 1; 

        datacon = []; 

        basename = folderNames{subindex}; 

        if condition == 21 
            filenametoload_121 = [basename(end-2:end) '_121_EEG.mat'];
            load(filenametoload_121);
            data121 = double(EEG.data);
        end

        filenametoload = [basename(end-2:end) '_' num2str(condition) '_EEG.mat'];

        load(filenametoload);

        datacon = double(EEG.data);

        if condition == 21 
            datacon = cat(3, datacon, data121);
        end

        [trialamp15Hz,winmat3d15Hz,phasestabmat15Hz,trialSNR15Hz] = freqtag_slidewin(datacon, 0, 1:121, 961:2158, 15, 600, 600, filenametoload);

        [amp3dwin, freqs3dwin, fftcompwin] = freqtag_FFT3D(winmat3d15Hz, 600);
        
       %  plot(freqs3dwin(1:20),amp3dwin(:,1:20)'); pause(1)

   %    topomap(amp3dwin(:, 5), [0 max(amp3dwin(:, 5))]); pause(20)
       % topomap(amp3dwin(:, 5), [0 0.19]); pause(20)

         
        ssVEP_ampmat(:,conditioncounter) = amp3dwin(:, 5); 

        conditioncounter = conditioncounter+1;


    end % condition

    outname = [basename(end-2:end) '_ssVEP_amp.mat']

    save(outname,'ssVEP_ampmat','-mat')


    cd ..

    pause(.10)

    cd ..

    pause(.5)

    fclose('all');
    close('all');

end % subject
