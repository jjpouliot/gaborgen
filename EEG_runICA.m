function EEG_runICA(partID, exclude4ICA, parentFolder, day1, day2)
% setting paths
if nargin < 3
    error('Set a parentFolder')
end

dataFolder = [parentFolder '/raw_data'];
if ~(exist(dataFolder) == 7)
    error('raw_data folder not inside parentFolder')
end

if nargin < 4
    day1 = 1;
end

if nargin < 5
    day2 = 1;
end

if day1 == 0 && day2 == 0
    error('day1 and day2 arguments are zero which means neither day should be processed')
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

        % initialize eeglab
        [ALLEEG, ~, ~, ~] = eeglab;

        % load dataset
        disp('Step 1/3 - load EEG data');
        currentDir =  [dataFolder '/' matchingDirs{j} '/EEG'];

        currentFilenames = {dir(currentDir).name};
        EEGpreICAIndex = find(endsWith(currentFilenames, '_02_prepped4ICA.set'));
        if ~isempty(EEGpreICAIndex)
            EEGpreICAFileName = currentFilenames{EEGpreICAIndex};
        elseif EEGIndex > 1
            error(['More than one 02_prepped4ICA.set file found in ' currentDir]);
        else
            error(['No 02_prepped4ICA.set file found in ' currentDir]);
        end
        [~, EEGFileName, ~] = fileparts(currentFilenames{EEGpreICAIndex});


        EEG = pop_loadset('filename', EEGpreICAFileName, 'filepath', currentDir);
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG, 0);

        % list indices of channels to include
        chanList = struct2cell(EEG.chanlocs);
        chanList = chanList(1,1,:);

        exclude4ICA_ind = find(ismember(squeeze(chanList), exclude4ICA{partI}));

        chans2include = 1:31;
        exclude4ICA_ind = unique([exclude4ICA_ind; 32]);
        chans2include = chans2include(~ismember(chans2include, exclude4ICA_ind));

        % run ICA
        disp('Step 2/3 - run ICA');
        EEG = pop_runica(EEG,'icatype','sobi','chanind',chans2include);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', ...
            [EEGFileName '_ICA'],'gui','off');
        [~, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

        % save data in eeglab format
        disp('Step 3/3 - save data with IC weights');
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename', [EEGFileName '_03_ICA.set'], ...
            'filepath', currentDir);

        % generate logfile
        logText = strcat('logfile for gaborgen_mri_eeg: run ICA\n', ...
            'date_time: ', string(datetime()), '\n', ...
            'participant: ', EEGFileName, '\n', ...
            'channels excluded from ICA: ', sprintf('%s ',string(exclude4ICA{partI})));
        fID = fopen([currentDir '/log02_runICA_' EEGFileName '.txt'], 'w');
        fprintf(fID, logText);
        fclose(fID);
    end
end
end