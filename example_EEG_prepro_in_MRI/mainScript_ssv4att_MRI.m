%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAIN SCRIPT SSV4ATT_MRI %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preprocessing, Step 1
% correct for scanner & cardioballistic artifacts, downsample, and filter
EEG_01_prep4ICA([101:102 104:107],500,[3 40],'/Volumes/TOSHIBA_4TB/ssV4att_MRI');

%%% manual step: check which channels should not be included in ICA



%% Preprocessing, Step 2
% run ICA, while excluding bad channels
EEG_02_runICA([101:102 104:107], ...
    {{'CP1','FC2'}, {'CP1','FC2'}, {}, {'F7','Cz'}, {}, {'T7'}}, ...
    '/Volumes/TOSHIBA_4TB/ssV4att_MRI')

% NOTES: 105 still riddled with gradient artifacts

%%% manual step: check which ICs to remove and which channels to interpolate



%% Preprocessing, Step 3
% exclude artifact ICs, interpolate bad channels, apply average reference,
% run artifact rejection, apply CSD

% 101 ICs - Gradients: 1,2,6,9,15; eyes: 3,5; CBA: 4
EEG_03_finishPrepro(101,{[1 2 3 4 5 6 9 15]},{{'CP1','FC2'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')

% 102 ICs - Gradients: 1,6; eyes: 3; CBA: 2,7
EEG_03_finishPrepro(102,{[1 2 3 6 7]},{{'CP1','FC2'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')

% 104 ICs - Gradients: 3,5; eyes: 2,10; CBA: 6,7
EEG_03_finishPrepro(104,{[2,3,5,6,7,10]},{{}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')

% 105 ICs - Gradients: 1,2,3; eyes: 4,10; CBA: 5
EEG_03_finishPrepro(105,{[1 2 3 4 5 10]},{{'F7','Cz'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')

% 106 ICs - Gradients: 1,2,4,7,8,14,18,24; eyes: 3; CBA: 9,11
EEG_03_finishPrepro(106,{[1 2 3 4 7 8 9 11 14 18 24]},{{}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')

% 107 ICs - Gradients: 1,8; eyes: 2; CBA: 4
EEG_03_finishPrepro(107,{[1 2 4 8]},{{'T7'}},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')



%% Sliding Window
EEG_04_slidingWindow([101:102 104:107], ...
                     {{'S 11', 'S111'},{'S 12', 'S112'},{'S 21', 'S121'},{'S 22', 'S122'}},{'11','12','21','22'}, ...
                     480,[-1600 4000],[-1600 -1400],[-1400 -466],[-200 0],[933 3734], ...
                     {'Oz','O1','O2'},'/Volumes/TOSHIBA_4TB/ssV4att_MRI')