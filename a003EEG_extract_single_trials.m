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

% participantIDs = [110,111,112,113,114,115,116,117,118,119,120,121,122,123,125,126,127,128,129,130];
% participantIDs = [131,132,133,134,135,136,139,140,141,142,143];
participantIDs = [110,111,112,113,114,115,116,117,118,119,120,121,122,123,125,126,...
                            127,128,129,130,131,132,133,134,135,136,139,140,141,142,143];
epochMs = [0 2000];
sampleRateHz = 500;
resampledRateHz = 600; % Check if this makes sense
postEpochedSamplePoints = 1:((epochMs(2) - epochMs(1))*(sampleRateHz/1000));
% epochSamplePoints = 1:(epochMs(2)*(sampleRateHz/1000));
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
        elseif EEGIndex > 1
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

        % urevent has all possible event triggers sent, this will have all
        % the stimuli triggers sent. The event field will have events that
        % are retained after preprocessing. After epoching, even more will
        % events will be lost; and importantly, if the epoch extends
        % infront of the trigger, then the stimulus event will be lost. So
        % the following code saves the index for the urevent and then uses
        % the retained epoch indices to find the final trials kept.
        ureventTable = struct2table(EEG.urevent);
        eventTable = struct2table(EEG.event);

        allPossibleEvents = {EEG.urevent.type};
        allPossibleEventIndices = 1:length(allPossibleEvents);

        allEventsIndicesLeft = {EEG.event.urevent};
        allEventsTypesLeft = {EEG.event.type};
        allEventsIndicesLeft = {EEG.event.urevent};
        
        % Use regular expressions to find stimuli events
        matchIdx = ~cellfun('isempty', regexp(allPossibleEvents, startOfStimMarkerRegEx, 'once'));
        allPossibleStimuliIndices = allPossibleEventIndices(matchIdx);

        allStimuliTriggersSent = allPossibleEvents(matchIdx);

        trialNumber = 1:length(allStimuliTriggersSent);

        matchIdx = ~cellfun('isempty', regexp(allEventsTypesLeft, startOfStimMarkerRegEx, 'once'));
        currentStimuliIndicesLeft = allEventsIndicesLeft(matchIdx);


        [EEG, acceptedEventIndices] = pop_epoch(EEG, startOfStimMarkerRegEx, ...
            epochMs ./ 1000, 'newname', 'segmented', 'epochinfo', 'yes');
        [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'gui','off');

        finalStimuliIndicesLeft = currentStimuliIndicesLeft(acceptedEventIndices);
        finalStimuliIndicesLeft = cell2mat(finalStimuliIndicesLeft);

        cleanSegmentsBool = ismember(allPossibleStimuliIndices, finalStimuliIndicesLeft);

        cleanTrialNumber = trialNumber(cleanSegmentsBool);


%         %% list the original trial indices from urevents for each trial category
%         eventTable = struct2table(EEG.event);
%         ureventTable = struct2table(EEG.urevent);
%         %cleanSegmentsBool = cell(1);
% 
%         eventInd = [];
%         ureventInd = [];
%         for i = 1:length(ureventTable.type)
%             if regexp(ureventTable.type{i}, startOfStimMarkerRegEx, 'once')
%                 ureventInd = [ureventInd, i];
%             end
%         end
%         for i = 1:length(eventTable.type)
%             if regexp(eventTable.type{i}, startOfStimMarkerRegEx, 'once')
%                 eventInd = [eventInd, eventTable.urevent(i)];
%             end
%         end
%         cleanSegmentsBool = ismember(ureventInd, eventInd);
% 
%         trialNumber = 1:length(ureventInd);
% 
%         cleanTrialNumber = trialNumber(cleanSegmentsBool);

        % sliding window post-cue - 15 Hz
        % [trialamp,winmat3d,phasestabmat,trialSNR] = freqtag_slidewin(data, plotflag, bslvec, ssvepvec, foi, sampnew, fsamp, outname)
        % bslvec and ssvepvec have to be the indices of the already
        % extracted epoch.
        [trialamp15Hz, winmat3d15Hz, phasestabmat15Hz, trialSNR15Hz] = ...
            freqtag_slidewin(EEG.data, 0, postEpochedSamplePoints, postEpochedSamplePoints, ...
            15, resampledRateHz, EEG.srate, 'whatever.txt');





        % Save the raw and sliding window timeseries and FFT by trial so
        % that I know which trials are missing
         mkdir([rawDataPath '/single_trial_timeseries_FFTs_CSD'])

         matSavePath = [rawDataPath '/single_trial_timeseries_FFTs_CSD/'];

        for trialIndex = 1:length(cleanTrialNumber)
            participantID = participantIDs(participantIndex);
            currentTrialNumber = cleanTrialNumber(trialIndex);

            % FFTs of raw and sliding window corrected data
            [rawAmp, rawFreqs, rawFFTcomp] = ...
                freqtag_FFT3D(EEG.data(:,:,trialIndex), EEG.srate);

            [slidingWindowAmp, slidingWindowFreqs, slidingWindowFFTcomp] = ...
                freqtag_FFT3D(winmat3d15Hz(:,:,trialIndex), resampledRateHz);

            rawTimeseriesName = ...
                ['rawChanTime_participant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{finalStimuliIndicesLeft(trialIndex)} ...
                '_sampleRate' num2str(EEG.srate) 'Hz.mat'];

            rawTimeseriesName = regexprep(rawTimeseriesName, '\s', '');

            SlidingWindowTimeseriesName = ...
                ['slidingWindowChanTime_participant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{finalStimuliIndicesLeft(trialIndex)} ...
                '_sampleRate' num2str(resampledRateHz) 'Hz.mat'];

            SlidingWindowTimeseriesName = regexprep(SlidingWindowTimeseriesName, '\s', '');

            fullPath = [matSavePath rawTimeseriesName];
            currentRawTimeSeries = EEG.data(:, :, trialIndex) ;

            save(fullPath, 'currentRawTimeSeries');

            fullPath = [matSavePath SlidingWindowTimeseriesName];
            currentSlidingWindowTimeSeries = winmat3d15Hz(:, :, trialIndex) ;

            save(fullPath, 'currentSlidingWindowTimeSeries');






            rawFFTName = ...
                ['rawChanFrequency_participant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{finalStimuliIndicesLeft(trialIndex)} ...
                '_sampleRate' num2str(EEG.srate) 'Hz.mat'];

            rawFFTName = regexprep(rawFFTName, '\s', '');

            SlidingWindowFFTName = ...
                ['slidingWindow_ChanFrequency_participant' num2str(participantID) '_trial' ...
                num2str(currentTrialNumber) '_condition'...
                ureventTable.type{finalStimuliIndicesLeft(trialIndex)} ...
                '_sampleRate' num2str(resampledRateHz) 'Hz.mat'];

            SlidingWindowFFTName = regexprep(SlidingWindowFFTName, '\s', '');

            fullPath = [matSavePath rawFFTName];
            %This was likely the wrong way to do this, more complex than
            %just switching the complex components to be absolute values
            %Now just do fft on single trials and keep the rawAmp
%             currentRawFFT = abs(rawFFTcomp(:, :, trialIndex));

            save(fullPath, 'rawAmp', 'rawFreqs');

            fullPath = [matSavePath SlidingWindowFFTName];
%             currentSlidingWindowFFT = abs(slidingWindowFFTcomp(:, :, trialIndex));

            save(fullPath, 'slidingWindowAmp', 'slidingWindowFreqs');
        end
    end
end