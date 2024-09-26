%% This will clear and edit your matlab path
% possibly necessary if you don't have eeglab set up correctly
restoredefaultpath
gaborgenCodeRepository = '/home/andrewf/Repositories/gaborgen';
eeglabDirectory = '/home/andrewf/Repositories/eeglab2024.0';
cd(eeglabDirectory)
[AllEEG, ~, ~, ~] = eeglab;
cd(gaborgenCodeRepository)

% Add EMGS directory and all subdirectories to path
emegs28path = '/home/andrewf/Repositories/emegs2.8';
addpath(genpath(emegs28path), '-end');

addpath(genpath('/home/andrewf/Repositories/freqTag'), '-end');


%% Load participant EEG, find good trials, extract 15Hz ssVEP

participantIDs = [110,111,112,113,114,115,116,117,118,119,120,121,122,123,125,126,127,128,129,130];
epochMs = [0 2000];
sampleRateHz = 500;
resampledRateHz = 510; % Check if this makes sense
epochSamplePoints = 1:(epochMs(2)*(sampleRateHz/1000));
rawDataPath = '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI';
startOfStimMarkerRegEx = '^S\s[1234]|^S121';



for participantIndex = 1:length(participantIDs)

    [currentParticipantDirectories, dataFolder, gaborgenCodeRepository] = ...
        gaborgenMriReturnDirs(participantIDs(participantIndex), rawDataPath, 1, 0);

    for j = 1:length(currentParticipantDirectories)

        %% set random number generator (just in case)
        rng(1);

        %% initialize eeglab
        [ALLEEG, ~, ~, ~] = eeglab;

        %% load dataset
        currentDirectory =  [dataFolder '/' currentParticipantDirectories{j} '/EEG/'];

        currentFilenames = {dir(currentDirectory).name};
        EEGIndex = find(endsWith(currentFilenames, '_04_preprocessed.set'));
        if length(EEGIndex) == 1
            EEGpreproFileName = currentFilenames{EEGIndex};
        elseif EEGcurrentDirectoryIndex > 1
            error(['More than one 04_preprocessed.set file found in ' currentDirectory]);
        else
            error(['No 04_preprocessed.set file found in ' currentDirectory]);
        end
        [~, EEGFileName, ~] = fileparts(currentFilenames{EEGIndex});


        EEG = pop_loadset('filename', EEGpreproFileName, 'filepath', currentDirectory);


        [ALLEEG, EEG, ~] = eeg_store(AllEEG, EEG, 0);

        chanList = struct2cell(EEG.chanlocs);
        chanList = chanList(1,1,:);

        EEG = eeg_checkset(EEG);
        EEG = pop_epoch(EEG, startOfStimMarkerRegEx, ...
            epochMs ./ 1000, 'newname', 'segmented', 'epochinfo', 'yes');
        [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'gui','off');

        %% list the original trial indices from urevents for each trial category
        eventTable = struct2table(EEG.event);
        ureventTable = struct2table(EEG.urevent);
        %cleanSegmentsBool = cell(1);

        eventInd = [];
        ureventInd = [];
        for i = 1:length(ureventTable.type)
            if regexp(ureventTable.type{i}, startOfStimMarkerRegEx, 'once')
                ureventInd = [ureventInd, i];
            end
        end
        for i = 1:length(eventTable.type)
            if regexp(eventTable.type{i}, startOfStimMarkerRegEx, 'once')
                eventInd = [eventInd, eventTable.urevent(i)];
            end
        end
        cleanSegmentsBool = ismember(ureventInd, eventInd);

        trialNumber = 1:length(ureventInd);

        cleanTrialNumber = trialNumber(cleanSegmentsBool);

        % sliding window post-cue - 15 Hz
        [trialamp15Hz, winmat3d15Hz, phasestabmat15Hz, trialSNR15Hz] = ...
            freqtag_slidewin(EEG.data, 0, epochSamplePoints, epochSamplePoints, ...
            15, resampledRateHz, EEG.srate, 'whatever.txt');

        % FFTs of raw and sliding window corrected data
        [rawAmp, rawFreqs, rawFFTcomp] = ...
            freqtag_FFT3D(EEG.data, 500);

        [slidingWindowAmp, slidingWindowFreqs, slidingWindowFFTcomp] = ...
            freqtag_FFT3D(winmat3d15Hz, resampledRateHz);



        % Save the raw and sliding window timeseries and FFT by trial so
        % that I know which trials are missing
         mkdir([rawDataPath '/single_trial_timeseries_FFTs'])

         matSavePath = [rawDataPath '/single_trial_timeseries_FFTs/'];

        for trialIndex = 1:length(cleanTrialNumber)
            participantID = participantIDs(participantIndex);
            currentTrialNumber = cleanTrialNumber(trialIndex);

            rawTimeseriesName = ...
                ['rawChanTime_participant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{eventInd(trialIndex)} ...
                '_sampleRate' num2str(EEG.srate) 'Hz.mat'];

            rawTimeseriesName = regexprep(rawTimeseriesName, '\s', '');

            SlidingWindowTimeseriesName = ...
                ['slidingWindowChanTime_participant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{eventInd(trialIndex)} ...
                '_sampleRate' num2str(resampledRateHz) 'Hz.mat'];

            SlidingWindowTimeseriesName = regexprep(SlidingWindowTimeseriesName, '\s', '');

            fullPath = [matSavePath rawTimeseriesName];
            currentRawTimeSeries = EEG.data(:, :, trialIndex) ;

            save(fullPath, 'currentRawTimeSeries');

            fullPath = [matSavePath SlidingWindowTimeseriesName];
            currentSlidingWindowTimeSeries = winmat3d15Hz(:, :, trialIndex) ;

            save(fullPath, 'currentSlidingWindowTimeSeries');






            rawFFTName = ...
                ['rawChanFrequency_particpant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{eventInd(trialIndex)} ...
                '_sampleRate' num2str(EEG.srate) 'Hz.mat'];

            rawFFTName = regexprep(rawFFTName, '\s', '');

            SlidingWindowFFTName = ...
                ['slidingWindow_ChanFrequency_particpant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{eventInd(trialIndex)} ...
                '_sampleRate' num2str(resampledRateHz) 'Hz.mat'];

            SlidingWindowFFTName = regexprep(SlidingWindowFFTName, '\s', '');

            fullPath = [matSavePath rawFFTName];
            currentRawFFT = abs(rawFFTcomp(:, :, trialIndex));

            save(fullPath, 'currentRawFFT', 'rawAmp', 'rawFreqs');

            fullPath = [matSavePath SlidingWindowFFTName];
            currentSlidingWindowFFT = abs(slidingWindowFFTcomp(:, :, trialIndex));

            save(fullPath, 'currentSlidingWindowFFT', 'slidingWindowAmp', 'slidingWindowFreqs');
        end
    end
end