%% This will clear and edit your matlab path
% possibly necessary if you don't have eeglab set up correctly
clear
restoredefaultpath
gaborgenCodeRepository = '/home/andrewf/Repositories/gaborgen';
eeglabDirectory = '/home/andrewf/Repositories/eeglab2024.0';
% gaborgenCodeRepository = 'C:\Users\jcedi\OneDrive\Documents\GitHub\gaborgen';
% eeglabDirectory = 'C:\Users\jcedi\OneDrive\Documents\eeglab2024.2';
cd(eeglabDirectory)
[AllEEG, ~, ~, ~] = eeglab;
cd(gaborgenCodeRepository)

% Add EMGS directory and all subdirectories to path
emegs28path = '/home/andrewf/Repositories/emegs2.8'; 
% emegs28path = 'C:\Users\jcedi\OneDrive\Documents\GitHub\EMEGShelper'; This should not be the EMEGShelper package, that is different than emegs2.8 
addpath(genpath(emegs28path), '-end'); 

addpath(genpath('/home/andrewf/Repositories/freqTag'), '-end');
% addpath(genpath('C:\Users\jcedi\OneDrive\Documents\GitHub\freqTag'), '-end');


%% Preprocessing, Step 1
% correct for scanner & cardioballistic artifacts, downsample, and filter
% EEG_prep4ICA(partID,newSamplingRate,filterCutOff,parentFolder)
% Data should be dropbox organization: parentFolder/raw_data/...partID/EEG
%124 bad, may not have T1 triggers, not sure though
EEG_prep4ICA([131,132,133,134,135,136,139,140,141,142,143], ...
    500,[3 40], '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI',1,0);

%%% manual step: check which channels should not be included in ICA



%% Preprocessing, Step 2
% run ICA, while excluding bad channels
% Example:
% EEG_runICA([101:102 104:107], ...
%     {{'CP1','FC2'}, {'CP1','FC2'}, {}, {'F7','Cz'}, {}, {'T7'}}, ...
%     '/Volumes/TOSHIBA_4TB/ssV4att_MRI')
% 

%EEG_runICA(partID, exclude4ICA, parentFolder, day1, day2)
% day1 and day2 split here to exclude different channels
%112 all bad after 600 seconds
%115 most channels look identical so there could be a reference issue
%116 seem to lose  F4, C4, Fz at 1862 seconds in, but only until 1905 or so
%119 ECG seems to glitch out often, Fz out at 1501   hold:
%121 C3 out at 380ish back at 934 and then gone again at 1090 haha, FC6 fading at 1390, may be possible to save channels with ICA
%123 T7, T8 impossibly large and noisy, maybe P and F are too
%128 ECG is there but drifts away for some reason
% EEG_runICA([110,111,    112,                                  113,   114,115,  116,117,118,119,   120,121,122,123,         125,126,127,128,129,130], ...
%                      {{},   {'FC1'},{'FC3', 'FC1', 'Cz', 'C3'},    {'C4'}, {},    {'F8'},{},    {},   {},    {'Fz'}, {},    {},    {},   {'T8','T7'},{},    {},    {},   {},    {},    {}}, ...
%     '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI', 1, 0)

% Alot of 142 is bad
EEG_runICA([131,132,133,134,135,                                                          136,                                                                            139,      140,141,142,      143], ...
                     {{},    {},    {},   {},    {'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8'}, {'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'FC2', 'Cz'}, {'FC1'}, {},    {},    {'CP1'}, {}}, ... #not sure if 139 channel is actually bad
    '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI', 1, 0)

%%% manual step: check which ICs to remove and which channels to interpolate
% How do you know which ICs to manuall remove? Andrew


%% Preprocessing, Step 3
% exclude artifact ICs, interpolate bad channels, apply average reference,
% run artifact rejection, apply CSD

%EEG_finishPrepro(partID, excludeICs, interpolateChans, parentFolder, day1, day2,CSDtransform)

EEG_finishPrepro([110,111,    112,                                  113,   114,115,  116,117,118,119,   120,121,122,123,         125,126,127,128,129,130], ...
                              {[],    [],        [],                                     [],       [],    [],     [],    [],    [],   [],       [],    [],   [],    [],            [],    [],    [],   [],    [],    []}, ...
                              {{},   {'FC1'},{'FC3', 'FC1', 'Cz', 'C3'},    {'C4'}, {},    {'F8'},{},    {},   {},    {'Fz'}, {},    {},    {},   {'T8','T7'},{},    {},    {},   {},    {},    {}}, ...
    '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI',1,0,1)

EEG_finishPrepro([131,132,133,134,135,                                                          136,                                                                            139,      140,141,142,      143], ...
                              {[],    [],   [],    [],    [],                                                             [],                                                                                [],          [],   [],    [],          []},...
                              {{},    {},    {},   {},    {'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8'}, {'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'FC2', 'Cz'}, {'FC1'}, {},    {},    {'CP1'}, {}}, ... #not sure if 139 channel is actually bad
    '/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI', 1, 0,1)


% % 101 ICs - Gradients: 1,2,6,9,15; eyes: 3,5; CBA: 4
% EEG_03_finishPrepro(101,{[1 2 3 4 5 6 9 15]},{{'CP1','FC2'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')
% 
% % 102 ICs - Gradients: 1,6; eyes: 3; CBA: 2,7
% EEG_03_finishPrepro(102,{[1 2 3 6 7]},{{'CP1','FC2'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')
% 
% % 104 ICs - Gradients: 3,5; eyes: 2,10; CBA: 6,7
% EEG_03_finishPrepro(104,{[2,3,5,6,7,10]},{{}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')
% 
% % 105 ICs - Gradients: 1,2,3; eyes: 4,10; CBA: 5
% EEG_03_finishPrepro(105,{[1 2 3 4 5 10]},{{'F7','Cz'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')
% 
% % 106 ICs - Gradients: 1,2,4,7,8,14,18,24; eyes: 3; CBA: 9,11
% EEG_03_finishPrepro(106,{[1 2 3 4 7 8 9 11 14 18 24]},{{}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')
% 
% % 107 ICs - Gradients: 1,8; eyes: 2; CBA: 4
% EEG_03_finishPrepro(107,{[1 2 4 8]},{{'T7'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')



%% Sliding Window
% This finds category averages, Andrew uses the next script a003 to find
% single trials

% function EEG_slidingWindow(partID, 
%     markerStrings, ... 
%     condStrings, ...
%     resampleTo, segTimesMs, blTimesMs_bl, ...
%     ssvepTimesMs_bl, blTimesMs_cue, ssvepTimesMs_cue, ...
%     avgChannels, parentFolder, day1, day2)

EEG_slidingWindow([110], ...
                     {{'S 11'},{'S 12'},{'S 13'},{'S 14'},{'S121'},{'S 21'},{'S 22'},{'S 23'},{'S 24'}}, ...
                     {'11','12','13','14','121','21','22','23','24'}, ...
                     510,[-1600 2000],[-1600 -1400], ...
                     [-1400 -466],[1 200],[0 1800], ...
                     {'Oz','O1','O2'},'/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI',1,0)

EEG_slidingWindow([111,112,113,114,115,116,117,118,119120,121,122,123,125,126,127,128,129,130], ...
                     {{'S 11'},{'S 12'},{'S 13'},{'S 14'},{'S121'},{'S 21'},{'S 22'},{'S 23'},{'S 24'},{'S 31'},{'S 32'},{'S 33'},{'S 34'}}, ...
                     {'11','12','13','14','121','21','22','23','24','31','32','33','34'}, ...
                     510,[-1600 2000],[-1600 -1400], ...
                     [-1400 -466],[1 200],[0 1800], ...
                     {'Oz','O1','O2'},'/home/andrewf/Research_data/EEG/Gaborgen24_EEG_fMRI',1,0)


% Example:
% EEG_slidingWindow([101:102 104:107], ...
%                      {{'S 11', 'S111'},{'S 12', 'S112'},{'S 21', 'S121'},{'S 22', 'S122'}},{'11','12','21','22'}, ...
%                      480,[-1600 4000],[-1600 -1400],[-1400 -466],[-200 0],[933 3734], ...
%                      {'Oz','O1','O2'},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')