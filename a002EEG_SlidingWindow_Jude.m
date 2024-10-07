

%% This will clear and edit your matlab path
% possibly necessary if you don't have eeglab set up correctly
restoredefaultpath
gaborgenCodeRepository = 'C:\Users\jcedi\OneDrive\Documents\GitHub\gaborgen';
eeglabDirectory = 'C:\Users\jcedi\OneDrive\Documents\eeglab2024.2';
cd(eeglabDirectory)
[AllEEG, ~, ~, ~] = eeglab;
cd(gaborgenCodeRepository)

% Add EMGS directory and all subdirectories to path
emegs28path = 'C:\Users\jcedi\OneDrive\Documents\GitHub\EMEGShelper'; 
addpath(genpath(emegs28path), '-end'); 

addpath(genpath('C:\Users\jcedi\OneDrive\Documents\GitHub\freqTag'), '-end');



%% Sliding Window

% function EEG_slidingWindow(partID, 
%     markerStrings, ... 
%     condStrings, ...
%     resampleTo, segTimesMs, blTimesMs_bl, ...
%     ssvepTimesMs_bl, blTimesMs_cue, ssvepTimesMs_cue, ...
%     avgChannels, parentFolder, day1, day2)

EEG_slidingWindowJ([110 111], ...
                     {{'S 11'},{'S 12'},{'S 13'},{'S 14'},{'S121'},{'S 21'},{'S 22'},{'S 23'},{'S 24'},{'S 31'},{'S 32'},{'S 33'},{'S 34'}}, ...
                     {'11','12','13','14','121','21','22','23','24','31','32','33','34'}, ...
                     500,[-1600 2000],[-1600 -1400], ...
                     [-1400 -466],[1 200],[0 1800], ...
                     {'Oz','O1','O2'},'C:\Users\jcedi\OneDrive\Documents\EEGStudies\Gaborgen24_EEG_fMRI',1,0)


