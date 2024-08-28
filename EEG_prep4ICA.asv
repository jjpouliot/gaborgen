function EEG_prep4ICA(partID, newSamplingRate, filterCutOff, parentFolder, day1, day2)
% set up and error checking
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
    error('day1 and day2 are zero which means neither day should be processed')
end

gaborgenCodeRepository = fileparts(mfilename('fullpath'));

participantDirectories = dir(dataFolder);
participantDirectories = participantDirectories(~ismember({participantDirectories.name}, {'.', '..'}));

for partI = 1:length(partID)
    % set random number generator (just in case)
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

        % initialize eeg lab
        [ALLEEG, ~, ~, ~] = eeglab;

        % load dataset
        disp('Step 1/10 - import BVA data');
        currentDir =  [dataFolder '/' matchingDirs{j} '/EEG'];

        currentFilenames = {dir(currentDir).name};
        vhdrIndex = find(endsWith(currentFilenames, '.vhdr'));
        if ~isempty(vhdrIndex)
            vhdrFileName = currentFilenames{vhdrIndex};
        elseif vhdrIndex > 1
            error(['More than one .vhdr file found in ' currentDir]);
        else
            error(['No .vhdr file found in ' currentDir]);
        end
        [~, EEGFileName, ~] = fileparts(currentFilenames{vhdrIndex});

        EEG = pop_loadbv(currentDir, vhdrFileName); % If line doesn't work, maybe get bva-io eeglab plugin
        [ALLEEG, EEG, ~] = pop_newset(ALLEEG, EEG, 0, 'setname', EEGFileName, 'gui', 'off');

        % add channel locations
        disp('Step 2/10 - add channel coordinates');
        EEG = eeg_checkset(EEG);
        EEG = pop_chanedit(EEG, 'lookup',[gaborgenCodeRepository '/standard_1005.elc']);
        %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

        % save raw data in eeglab format
        disp('Step 3/10 - save raw data in eeglab format');
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG, 'filename',[EEGFileName '_01_raw.set'], ...
            'filepath',currentDir);

        % rename TR markers
        disp('Step 4/10 - rename TR markers');
        for i = 1:length(EEG.event)
            if strcmp(EEG.event(i).type(1:4), 'T  1')
                EEG.event(i).type = 'T1';
            end
        end

        % remove MRI gradient artifacts, need FMRIB plugin
        disp('Step 5/10 = remove gradient artifacts (this may take a while...)');
        EEG = eeg_checkset(EEG);
        [EEG, ~] = pop_fmrib_fastr(EEG, 70, 4, 31, 'T1', 0, 0, 0, [], [], [], 32, 'auto');

        % downsample
        disp('Step 6/10 - downsample data');
        EEG = eeg_checkset(EEG);
        EEG = pop_resample(EEG, newSamplingRate);

        % detect QRS
        disp('Step 7/10 - detect R spikes');
        EEG = eeg_checkset(EEG);
        EEG = pop_fmrib_qrsdetect(EEG, 32, 'qrs', 'no');

        % remove cardioballistic artifacts
        disp('Step 8/10 - remove cardioballistic artifacts');
        EEG = eeg_checkset(EEG);
        EEG = pop_fmrib_pas(EEG,'qrs','obs',3);

        % filtering
        disp('Step 9/10 - filter data');
        EEG = eeg_checkset(EEG);
        EEG = pop_eegfiltnew(EEG, 'locutoff',filterCutOff(1),'hicutoff',filterCutOff(2),'plotfreqz',0,'channels',{'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','Oz','FC1','FC2','CP1','CP2','FC5','FC6','CP5','CP6','TP9','TP10','POz'});
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[EEGFileName '_prepped4ICA'],'gui','off');
        [~, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

        % save prepped data
        disp('Step 10/10 - save data prepped 4 ICA');
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename',[EEGFileName '_02_prepped4ICA.set'],'filepath',currentDir);

        % generate logfile
        logText = strcat('logfile for gaborgen_mri_eeg: prepping 4 ICA\n', ...
            'date_time: ', string(datetime()), '\n', ...
            'participant: ', EEGFileName, '\n', ...
            'new sampling rate: ', int2str(newSamplingRate), '\n', ...
            'filter cut-offs: ', int2str(filterCutOff));
        fID = fopen([currentDir '/log01_prep4ICA_' EEGFileName '.txt'], 'w');
        fprintf(fID, logText);
        fclose(fID);
    end
end
end