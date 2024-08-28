function EEG_finishPrepro(partID, excludeICs, interpolateChans, parentFolder, day1, day2)
%% setting paths
if nargin < 4
    error('Set a parentFolder')
end

dataFolder = [parentFolder '/raw_data'];
if ~(exist(dataFolder) == 7)
    error('raw_data folder not inside parentFolder')
end

if nargin < 5
    day1 = 1;
end

if nargin < 6
    day2 = 1;
end

if day1 == 0 && day2 == 0
    error('day1 and day2 arguments are zero which means neither day should be processed')
end

gaborgenCodeRepository = fileparts(mfilename('fullpath'));

participantDirectories = dir(dataFolder);
participantDirectories = participantDirectories(~ismember({participantDirectories.name}, {'.', '..'}));


for partI = 1:length(partID)
    %% set random number generator (just in case)
    rng(1);

    matchingDirs = {};
    for i = 1:length(participantDirectories)
        dirname = participantDirectories(i).name;
        if contains(dirname,num2str(partID(partI)))
            matchingDirs{end+1} = dirname;
        end
    end

    %day 1 or 2 exclusion
    if day1 == 0
        finalDirs = {};
        for i = 1:length(matchingDirs)
            if contains(matchingDirs{i}, 'DAY2')
                finalDirs{end+1} = matchingDirs{i};
            end
        end
        matchingDirs = finalDirs;
    end

    if day2 == 0
        finalDirs = {};
        for i = 1:length(matchingDirs)
            if ~contains(matchingDirs{i}, 'DAY2')
                finalDirs{end+1} = matchingDirs{i};
            end
        end
        matchingDirs = finalDirs;
    end

    if length(matchingDirs) > 2
        error('There are more than two directories with the same subject number which should be impossible.')
    end
    for j = 1:length(matchingDirs)

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
end