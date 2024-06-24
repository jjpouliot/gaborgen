function EEG_03_finishPrepro(partID,excludeICs,interpolateChans,parentFolder)
    %% setting paths
    if nargin < 3
        cd ..;
        parentFolder = pwd;
    end

    for partI = 1:length(partID)
        %% set random number generator (just in case)
        rng(1);
        
        %% initialize eeglab
        [ALLEEG, ~, ~, ~] = eeglab;
        
        %% load dataset
        disp('Step 1/7 - load EEG data');
        dataFolder = [parentFolder '/data/' int2str(partID(partI)) '/EEG/'];
        EEG = pop_loadset('filename',['ssv4att_MRI_' int2str(partID(partI)) '_03_ICA.set'], ...
                          'filepath',dataFolder);
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG, 0);

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
        [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
        
        %% transform into CSD
        disp('Step 6/7 - CSD transformation');
        EEG = eeg_checkset(EEG);
        % we have to remove non-EEG channels for the ERPlab plugin to work
        EEG = pop_select(EEG, 'rmchannel',{'ECG'});
        % do the CSD transformation
        EEG = csdFromErplabAutomated(EEG); % using defaults: m = 4 splines, lambda = 1e-5, 10 cm radius
        [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'gui','off');

        %% save data in eeglab format
        disp('Step 7/7 - save preprocessed data');
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename',['ssv4att_MRI_' int2str(partID(partI)) '_04_preprocessed.set'], ...
                         'filepath',dataFolder);

        %% generate logfile
        logText = strcat('logfile for ssv4att_MRI: finish preprocessing\n', ...
                   'date_time: ', string(datetime()), '\n', ...
                   'participant: ', int2str(partID(partI)), '\n', ...
                   'removed ICs: ', int2str(excludeICs{partI}), '\n', ...                   
                   'interpolated channels: ', sprintf('%s ',string(interpolateChans{partI})));
        fID = fopen([dataFolder '/log03_finishPrepro_' int2str(partID(partI)) '.txt'], 'w');
        fprintf(fID, logText);
        fclose(fID);
    end
end