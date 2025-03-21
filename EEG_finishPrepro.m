function EEG_finishPrepro(partID, excludeICs, interpolateChans, parentFolder, day1, day2, CSDtransform)

if nargin < 5
    day1 = 1;
end

if nargin < 6
    day2 = 1;
end

if nargin < 7
    CSDtransform = 1;
end

for partI = 1:length(partID)
    % set random number generator (just in case)
    rng(1);

    [currentParticipantDirectories, dataFolder, gaborgenCodeRepository] = ...
        gaborgenMriReturnDirs(partID(partI), parentFolder, day1, day2);

    for j = 1:length(currentParticipantDirectories)

        %% initialize eeglab
        [AllEEG, ~, ~, ~] = eeglab;

        %% load dataset
        disp('Step 1/7 - load EEG data');

        currentDirectory =  [dataFolder '/' currentParticipantDirectories{j} '/EEG'];

        currentFilenames = {dir(currentDirectory).name};
        EEGICAIndex = find(endsWith(currentFilenames, '_03_ICA.set'));
        if length(EEGICAIndex) == 1
            EEGICAFileName = currentFilenames{EEGICAIndex};
        elseif EEGICAIndex > 1
            error(['More than one 03_ICA.set file found in ' currentDirectory]);
        else
            error(['No 03_ICA.set file found in ' currentDirectory]);
        end
        [~, EEGFileName, ~] = fileparts(currentFilenames{EEGICAIndex});


        EEG = pop_loadset('filename', EEGICAFileName, 'filepath', currentDirectory);


        [AllEEG, EEG, ~] = eeg_store(AllEEG, EEG, 0);

        %% remove ICs
        disp('Step 2/7 - removing independent components from data');
        EEG = eeg_checkset(EEG);
        EEG = pop_subcomp(EEG, excludeICs{partI}, 0);

        %% list indices of channels to interpolate
        chanList = struct2cell(EEG.chanlocs);
        chanList = chanList(1,1,:);

        interpolateChans_ind = find(ismember(squeeze(chanList),interpolateChans{partI}));

        %% interpolate channels
        disp('Step 3/7 - interpolating channels');
        EEG = eeg_checkset(EEG);
        EEG = pop_interp(EEG, interpolateChans_ind, 'spherical');

        %% apply average reference
        disp('Step 4/7 - apply average reference');
        EEG = eeg_checkset(EEG);
        EEG = pop_reref( EEG, [],'exclude',32);

        %% remove artifacts
        disp('Step 5/7 - remove artifacts');
        EEG = eeg_checkset(EEG);
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
        [~, EEG, ~] = pop_newset(AllEEG, EEG, 1,'gui','off');

        %% transform into CSD
        disp('Step 6/7 - CSD transformation');
        EEG = eeg_checkset(EEG);
        % we have to remove non-EEG channels for the ERPlab plugin to work
        EEG = pop_select(EEG, 'rmchannel',{'ECG'});
        % do the CSD transformation
        if CSDtransform
            EEG = csdFromErplabAutomated(EEG); % using defaults: m = 4 splines, lambda = 1e-5, 10 cm radius
        end
        [~, EEG, ~] = pop_newset(AllEEG, EEG, 1,'gui','off');

        %% save data in eeglab format
        disp('Step 7/7 - save preprocessed data');
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename',[EEGFileName(1:(end - length('_03_ICA'))) '_04_preprocessed.set'], ...
            'filepath',currentDirectory);

        %% generate logfile
        logText = strcat('logfile for gaborgen_mri_eeg: finish PrePro\n', ...
            'date_time: ', string(datetime()), '\n', ...
            'participant: ', int2str(partID(partI)), '\n', ...
            'removed ICs: ', int2str(excludeICs{partI}), '\n', ...
            'interpolated channels: ', sprintf('%s ',string(interpolateChans{partI})));
        fID = fopen([currentDirectory '/log03_finishPrepro_' int2str(partID(partI)) '.txt'], 'w');
        fprintf(fID, logText);
        fclose(fID);
    end
end
end