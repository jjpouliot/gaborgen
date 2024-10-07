    function EEG_02_runICA(partID,exclude4ICA,parentFolder)
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
        disp('Step 1/3 - load EEG data');
        dataFolder = [parentFolder '/data/' int2str(partID(partI)) '/EEG/'];
        EEG = pop_loadset('filename',['ssv4att_MRI_' int2str(partID(partI)) '_02_prepped4ICA.set'], ...
                          'filepath',dataFolder);
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG, 0);
        
        %% list indices of channels to include
        chanList = struct2cell(EEG.chanlocs);
        chanList = chanList(1,1,:);

        exclude4ICA_ind = find(ismember(squeeze(chanList),exclude4ICA{partI}));

        chans2include = 1:31;
        exclude4ICA_ind = unique([exclude4ICA_ind; 32]);
        chans2include = chans2include(~ismember(chans2include,exclude4ICA_ind));
        
        %% run ICA
        disp('Step 2/3 - run ICA');
        EEG = pop_runica(EEG,'icatype','sobi','chanind',chans2include);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',['ssv4att_MRI_' int2str(partID(partI)) ' - ICA'],'gui','off');         
        [~, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %% save data in eeglab format
        disp('Step 3/3 - save data with IC weights');
        EEG = eeg_checkset(EEG);
        pop_saveset(EEG, 'filename',['ssv4att_MRI_' int2str(partID(partI)) '_03_ICA.set'], ...
                         'filepath',dataFolder);

        %% generate logfile
        logText = strcat('logfile for ssv4att_MRI: run ICA\n', ...
                   'date_time: ', string(datetime()), '\n', ...
                   'participant: ', int2str(partID(partI)), '\n', ...
                   'channels excluded from ICA: ', sprintf('%s ',string(exclude4ICA{partI})));
        fID = fopen([dataFolder '/log02_runICA_' int2str(partID(partI)) '.txt'], 'w');
        fprintf(fID, logText);
        fclose(fID);
    end
end