function EEG_04_slidingWindow(partID,markerStrings,condStrings,resampleTo,segTimesMs,blTimesMs_bl,ssvepTimesMs_bl,blTimesMs_cue,ssvepTimesMs_cue,avgChannels,parentFolder)
    %% setting paths
    if nargin < 11
        cd ..;
        parentFolder = pwd;
    end

    for partI = 1:length(partID)
         for markI = 1:length(markerStrings)
 
            disp(['now processing participant ' int2str(partID(partI)) ', condition ' condStrings{markI}])
            %% set random number generator (just in case)
            rng(1);
            
            %% initialize eeglab
            [ALLEEG, ~, ~, ~] = eeglab;
            
            %% load dataset
            dataFolder = [parentFolder '/data/' int2str(partID(partI)) '/EEG/'];
            EEG = pop_loadset('filename',['ssv4att_MRI_' int2str(partID(partI)) '_04_preprocessed.set'], ...
                              'filepath',dataFolder);
            [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG, 0);
    
            %% list indices of channels to average across
            chanList = struct2cell(EEG.chanlocs);
            chanList = chanList(1,1,:);

            avgChannels_ind = find(ismember(squeeze(chanList),avgChannels));

            %% resample
            EEG = pop_resample(EEG, resampleTo);

            %% extract epochs and remove baseline (whole epoch)
            EEG = eeg_checkset(EEG);
            EEG = pop_epoch(EEG, markerStrings{markI}, segTimesMs ./ 1000, 'newname', 'segmented', 'epochinfo', 'yes');
            [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 


%             %% list the original trial indices from urevents for each trial category
%             %disp('Step 2/6 - removing independent components from data');
%             eventTable = struct2table(EEG.event);
%             ureventTable = struct2table(EEG.urevent);
%             %cleanSegmentsBool = cell(1);
% 
%             eventInd = [];
%             ureventInd = [];
%             for i = 1:length(ureventTable.type)
%                 if ismember(ureventTable.type{i}, markerStrings{markI})
%                     ureventInd = [ureventInd, i];
%                 end
%             end
%             for i = 1:length(eventTable.type)
%                 if ismember(eventTable.type{i}, markerStrings{markI})
%                      eventInd = [eventInd, eventTable.urevent(i)];
%                 end
%             end
%             cleanSegmentsBool = ismember(ureventInd, eventInd);
% 
% 
%             %% determine baseline and data samples
%             blSamples_bl = (floor((blTimesMs_bl(1) - segTimesMs(1)) .* (EEG.srate/1000)) : ceil((blTimesMs_bl(2) - segTimesMs(1)) .* (EEG.srate/1000))) + 1;
%             ssvepSamples_bl = (floor((ssvepTimesMs_bl(1) - segTimesMs(1)) .* (EEG.srate/1000)) : ceil((ssvepTimesMs_bl(2) - segTimesMs(1)) .* (EEG.srate/1000))) + 1;
%             blSamples_cue = (floor((blTimesMs_cue(1) - segTimesMs(1)) .* (EEG.srate/1000)) : ceil((blTimesMs_cue(2) - segTimesMs(1)) .* (EEG.srate/1000))) + 1;
%             ssvepSamples_cue = (floor((ssvepTimesMs_cue(1) - segTimesMs(1)) .* (EEG.srate/1000)) : ceil((ssvepTimesMs_cue(2) - segTimesMs(1)) .* (EEG.srate/1000))) + 1;
% 
%             %% sliding window baseline - 8.57 Hz
%             [powmat_bl_857,winmat3d_bl_857,~,~] = freqtag_slidewin(EEG.data, 0, blSamples_bl, ssvepSamples_bl, 60/7, EEG.srate, EEG.srate, 'whatever.txt');
% 
%             %% sliding window baseline - 15 Hz
%             [powmat_bl_15,winmat3d_bl_15,~,~] = freqtag_slidewin(EEG.data, 0, blSamples_bl, ssvepSamples_bl, 15, EEG.srate, EEG.srate, 'whatever.txt');
% 
%             %% sliding window post-cue - 8.57 Hz
%             [powmat_cue_857,winmat3d_cue_857,~,~] = freqtag_slidewin(EEG.data, 0, blSamples_cue, ssvepSamples_cue, 60/7, EEG.srate, EEG.srate, 'whatever.txt');
% 
%             %% sliding window post-cue - 15 Hz
%             [powmat_cue_15,winmat3d_cue_15,~,~] = freqtag_slidewin(EEG.data, 0, blSamples_cue, ssvepSamples_cue, 15, EEG.srate, EEG.srate, 'whatever.txt');
% 
% 
%             %% create and save vectors for parametric modulation
%             powmat_blcorr_857 = powmat_cue_857 ./ powmat_bl_857;
%             powmat_blcorr_15 = powmat_cue_15 ./ powmat_bl_15;
%     
%             parmodvec_857_abs = NaN(length(ureventInd),1);
%             parmodvec_15_abs = NaN(length(ureventInd),1);
%             parmodvec_857_blcorr = NaN(length(ureventInd),1);
%             parmodvec_15_blcorr = NaN(length(ureventInd),1);
%     
%             parmodvec_857_abs(cleanSegmentsBool) = mean(powmat_cue_857(avgChannels_ind,:), 1);
%             parmodvec_15_abs(cleanSegmentsBool) = mean(powmat_cue_15(avgChannels_ind,:), 1);
%             parmodvec_857_blcorr(cleanSegmentsBool) = mean(powmat_blcorr_857(avgChannels_ind,:), 1);
%             parmodvec_15_blcorr(cleanSegmentsBool) = mean(powmat_blcorr_15(avgChannels_ind,:), 1);
%     
%             filename = [dataFolder, int2str(partID(partI)), '_parMod_', condStrings{markI}, '_857_abs.txt'];
%             writematrix(parmodvec_857_abs,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_parMod_', condStrings{markI}, '_15_abs.txt'];
%             writematrix(parmodvec_15_abs,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_parMod_', condStrings{markI}, '_857_blcorr.txt'];
%             writematrix(parmodvec_857_blcorr,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_parMod_', condStrings{markI}, '_15_blcorr.txt'];
%             writematrix(parmodvec_15_blcorr,filename);        
%     
%             
%             %% save single-trial amplitudes at driving frequencies
% 
%             %powmat_blcorr_857 = powmat_cue_857 ./ powmat_bl_857;
%             %powmat_blcorr_15 = powmat_cue_15 ./ powmat_bl_15;
%             
%             filename = [dataFolder, int2str(partID(partI)), '_fftST_', condStrings{markI}, '_857_abs.txt'];
%             writematrix(powmat_cue_857,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_fftST_', condStrings{markI}, '_15_abs.txt'];
%             writematrix(powmat_cue_15,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_fftST_', condStrings{markI}, '_857_blcorr.txt'];
%             writematrix(powmat_blcorr_857,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_fftST_', condStrings{markI}, '_15_blcorr.txt'];
%             writematrix(powmat_blcorr_15,filename);
% 
% 
%             %% compute and save ERPs (from sliding window)
%             erp_bl857 = mean(winmat3d_bl_857,3);
%             erp_bl15 = mean(winmat3d_bl_15,3);
%             erp_cue857 = mean(winmat3d_cue_857,3);
%             erp_cue15 = mean(winmat3d_cue_15,3);
%             
%             filename = [dataFolder, int2str(partID(partI)), '_erp_', condStrings{markI}, '_857_bl.txt'];
%             writematrix(erp_bl857,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_erp_', condStrings{markI}, '_15_bl.txt'];
%             writematrix(erp_bl15,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_erp_', condStrings{markI}, '_857_cue.txt'];
%             writematrix(erp_cue857,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_erp_', condStrings{markI}, '_15_cue.txt'];
%             writematrix(erp_cue15,filename);
% 
% 
% 
%             %% compute and save FFT on ERP (from sliding window)
%             fft_bl857 = freqtag_FFT(erp_bl857,EEG.srate);
%             fft_bl15 = freqtag_FFT(erp_bl15,EEG.srate);
%             fft_cue857 = freqtag_FFT(erp_cue857,EEG.srate);
%             fft_cue15 = freqtag_FFT(erp_cue15,EEG.srate);
%             fft_blcorr857 = fft_cue857 ./ fft_bl857;
%             fft_blcorr15 = fft_cue15 ./ fft_bl15;
% 
%             filename = [dataFolder, int2str(partID(partI)), '_erpfft_', condStrings{markI}, '_857_abs.txt'];
%             writematrix(fft_cue857,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_erpfft_', condStrings{markI}, '_15_abs.txt'];
%             writematrix(fft_cue15,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_erpfft_', condStrings{markI}, '_857_blcorr.txt'];
%             writematrix(fft_blcorr857,filename);
%             filename = [dataFolder, int2str(partID(partI)), '_erpfft_', condStrings{markI}, '_15_blcorr.txt'];
%             writematrix(fft_blcorr15,filename);


            %% compute and save Hilbert
            hilb_857 = freqtag_HILB_chris(double(mean(EEG.data,3)),60/7,8,20,0,EEG.srate);
            hilb_15 = freqtag_HILB_chris(double(mean(EEG.data,3)),15,8,20,0,EEG.srate);

            filename = [dataFolder, int2str(partID(partI)), '_hilbert_', condStrings{markI}, '_857.txt'];
            writematrix(hilb_857,filename);
            filename = [dataFolder, int2str(partID(partI)), '_hilbert_', condStrings{markI}, '_15.txt'];
            writematrix(hilb_15,filename);

         end


        %% generate logfile
        logText = strcat('logfile for ssv4att_MRI: sliding window\n', ...
                   'date_time: ', string(datetime()), '\n', ...
                   'participant: ', int2str(partID(partI)), '\n', ...
                   'new sampling rate: ', int2str(resampleTo), '\n', ...                   
                   'segment boundaries in ms: ', int2str(segTimesMs), '\n', ...
                   'baseline window for baseline segment (ms): ', int2str(blTimesMs_bl), '\n', ...
                   'ssVEP window for baseline segment (ms): ', int2str(ssvepTimesMs_bl), '\n', ...
                   'baseline window for post-cue segment (ms): ', int2str(blTimesMs_cue), '\n', ...
                   'ssVEP window for post-cue segment (ms): ', int2str(ssvepTimesMs_cue), '\n', ...
                   'channels averaged for parametric modulators: ', sprintf('%s ',string(avgChannels)));
        fID = fopen([dataFolder '/log04_slidingWindow_' int2str(partID(partI)) '.txt'], 'w');
        fprintf(fID, logText);
        fclose(fID);
    end
end