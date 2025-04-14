%% This will clear and edit your matlab path
% possibly necessary if you don't have eeglab set up correctly
restoredefaultpath;
clc
clear
close all

% gaborgenCodeRepository = '/Users/jcedielescobar/Documents/GitHub';
% eeglabDirectory = '/Users/jcedielescobar/Documents/MATLAB/eeglab2024.2';
gaborgenCodeRepository = '/home/andrewf/Repositories/gaborgen';
eeglabDirectory = '/home/andrewf/Repositories/eeglab2024.0';

cd(eeglabDirectory)
[AllEEG, ~, ~, ~] = eeglab;
cd(gaborgenCodeRepository)
% Add EMGS directory and all subdirectories to path
emegs28path = '/home/andrewf/Repositories/emegs2.8'; 
% emegs28path = 'C:\Users\jcedi\OneDrive\Documents\GitHub\EMEGShelper'; 
addpath(genpath(emegs28path), '-end'); 

% addpath(genpath('/Users/jcedielescobar/Documents/GitHub'), '-end');
%% Format transformation

% Define the main folder path
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';
% main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Define the specific folder to search
vhdr_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(vhdr_folder)
    % Get the list of subject folders within the data folder
    subjects = dir(vhdr_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 1:length(subjects) % Changed to iterate over all subjects
        subject = subjects(i).name;
        eeg_folder = fullfile(vhdr_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Read all .vhdr files in the EEG folder
            vhdr_files = dir(fullfile(eeg_folder, '*.vhdr'));
            
            for j = 1:length(vhdr_files)
                vhdr_file = vhdr_files(j).name;
                setFilePath = fullfile(eeg_folder, vhdr_file);
                
                % Load the .vhdr file
                EEG = pop_loadbv(eeg_folder, vhdr_file);
                
                % Define the output file name and path
                outputFileName = strrep(vhdr_file, '.vhdr', '_01_raw.set'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the saved file name to the results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The data folder does not exist.');
end


%% Chan Location, renaming T1, removing MRI artifacts, resampling
% Define the main folder path
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Define the specific folder to search
vhdr_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(vhdr_folder)
    % Get the list of subject folders within the raw_data folder
    subjects = dir(vhdr_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 1:length(subjects)
        subject = subjects(i).name;
        eeg_folder = fullfile(vhdr_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Read all previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_01_raw.set'));
            
            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);
                
                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
                
                % Load channel locations
                EEG = pop_chanedit(EEG, 'lookup', '/Users/jcedielescobar/Documents/GitHub/gaborgen1/standard_1005.elc');
                EEG.allchan = EEG.chanlocs; % Save the original locations
                
                % Rename events
                for indextrial = 1:length(EEG.event)
                    if strcmp(EEG.event(indextrial).type(1:4), 'T  1')
                        EEG.event(indextrial).type = 'T1';
                    end
                end
                
                % Remove MRI gradient artifacts
                disp('Remove gradient artifacts (this may take a while...)');
                EEG = eeg_checkset(EEG);
                [EEG, ~] = pop_fmrib_fastr(EEG, 70, 4, 31, 'T1', 0, 0, 0, [], [], [], 32, 'auto');
                
                % Resample
                EEG = eeg_checkset(EEG);
                EEG = pop_resample(EEG, 500);
                
                % Define the output file name and path
                outputFileName = strrep(set_file, '_01_raw.set', '_02_Ch_MRI_RS.set'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the saved file name to the results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The raw_data folder does not exist.');
end

%% Cardiobalistic artifacts

% Define the main folder path
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Define the specific folder to search
vhdr_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(vhdr_folder)
    % Get the list of subject folders within the raw_data folder
    subjects = dir(vhdr_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 1:length(subjects)
        subject = subjects(i).name;
        eeg_folder = fullfile(vhdr_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Read all previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_02_Ch_MRI_RS.set'));
            
            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);
                
                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
                
                % Detect QRS
                EEG = eeg_checkset(EEG);
                EEG = pop_fmrib_qrsdetect(EEG, 32, 'qrs', 'no');
                    
                % Remove cardiobalistic artifacts
                EEG = eeg_checkset(EEG);
                EEG = pop_fmrib_pas(EEG, 'qrs', 'obs', 3);
        
                % Define the output file name and path
                outputFileName = strrep(set_file, '_02_Ch_MRI_RS.set', '_03_CB.set'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the saved file name to the results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The raw_data folder does not exist.');
end

%% Filtering

main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Define the specific folder to search
vhdr_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(vhdr_folder)
    % Get the list of subject folders within the raw_data folder
    subjects = dir(vhdr_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));
 for i = 1:length(subjects)
        subject = subjects(i).name;
        eeg_folder = fullfile(vhdr_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Read all previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_03_CB.set'));
            
            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);
                
                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
                    
                % Filtering
                EEG = pop_eegfiltnew(EEG, 'locutoff', 3, 'hicutoff', 40, 'plotfreqz', 0, 'channels', {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'Fz', 'Cz', 'Pz', 'Oz', 'FC1', 'FC2', 'CP1', 'CP2', 'FC5', 'FC6', 'CP5', 'CP6', 'TP9', 'TP10', 'POz'});
 
                % Define the output file name and path
                outputFileName = strrep(set_file, '_03_CB.set', '_04_FL.set'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the saved file name to the results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The raw_data folder does not exist.');
end

%% bad channeles (check manually) Just for the first time. The
%interpolaton has ready the struct with the bad channels
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Initialize an array to store bad channels info
badChannelsStruct = struct('subjectName', {}, 'numBadChannels', {});

% Define the specific folder to search
vhdr_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(vhdr_folder)
    % Get the list of subject folders within the raw_data folder
    subjects = dir(vhdr_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 21:length(subjects)  % You might want to change the loop range
        subject = subjects(i).name;
        eeg_folder = fullfile(vhdr_folder, subject, 'EEG');

        if isfolder(eeg_folder)
            % Read all previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_04_FL.set'));

            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);

                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);

                % Remove bad channels (example: remove channel 32)
                EEG = pop_select(EEG, 'nochannel', [32]);
                eeglab redraw

                % pop_eegplot(EEG, 1, 1, 1);

                % Mark bad trials
                R1 = input('Highlight bad trials, update marks and then press enter');

                % Save bad channels
                answer = inputdlg('Enter bad channels', 'Bad channel removal', [1 31]);
                str = answer{1};
                badChannels = str2double(strsplit(str));

                % Save the number of bad channels and subject name in the structure
                badChannelsStruct(i).subjectName = subject;
                badChannelsStruct(i).numBadChannels = length(badChannels);

                close all;

                % Define the output file name and path
                outputFileName = strrep(set_file, '_04_FL.set', '_05_BC.set'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);

                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);

                % Add the saved file name to the results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The raw_data folder does not exist.');
end

%% ICA

ruta_principal = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
resultados = {};

% Define the specific folder to search
carpeta_vhdr = fullfile(ruta_principal, 'raw_data');

% Check if the folder exists
if isfolder(carpeta_vhdr)
    % Get the list of subject folders within the data folder
    sujetos = dir(carpeta_vhdr);
    sujetos = sujetos([sujetos.isdir] & ~ismember({sujetos.name}, {'.', '..'}));

    for i = 7:length(sujetos) 
        sujeto = sujetos(i).name;
        eeg_carpeta = fullfile(carpeta_vhdr, sujeto, 'EEG');
        
        if isfolder(eeg_carpeta)
            % Read all the previously generated .set files
            archivos_set = dir(fullfile(eeg_carpeta, '*_05_BC.set'));
            
            for j = 1:length(archivos_set)
                archivo_set = archivos_set(j).name;
                setFilePath = fullfile(eeg_carpeta, archivo_set);
                
                % Load the .set file
                EEG = pop_loadset('filename', archivo_set, 'filepath', eeg_carpeta);
                
                % Run ICA
                EEG = pop_runica(EEG,'icatype','sobi');
                
                % Define the name and path of the output file
                outputFileName = strrep(archivo_set, '_05_BC.set', '_06_ICA.set'); % Change the name
                outputFilePath = fullfile(eeg_carpeta, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_carpeta);
                
                % Add the name of the saved file to results
                resultados{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The data folder does not exist.');
end

%% Eye components (this label the components but you should check manually)
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Define the specific folder to search
data_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(data_folder)
    % Get the list of subject folders within the data folder
    subjects = dir(data_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 1:length(subjects) % Starting index for subjects (adjust as needed)
        subject = subjects(i).name;
        eeg_folder = fullfile(data_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Get all previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_06_ICA.set'));
            
            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);
                
                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
                
                %When the first time you want to label the components
                % EEG = pop_iclabel(EEG, 'default');

                % Load the eye_components.mat structure
                eyeComponentsPath = fullfile(main_path, 'eye_components.mat');
                
                if isfile(eyeComponentsPath)
                    % Load the structure from eye_components.mat
                    load(eyeComponentsPath, 'eye_components'); % Assuming the structure is called 'eye_components'
                    
                    % Find the subject index in the eye_components structure
                    subjectIndex = find(strcmp({eye_components.subjectName}, subject));
                    
                    if ~isempty(subjectIndex)
                        % Get the list of components to remove
                        componentsToRemove = eye_components(subjectIndex).removedComponents;
                        
                        % Remove the components from EEG
                        EEG = pop_subcomp(EEG, componentsToRemove, 0);
                    else
                        disp(['No eye components found for subject ', subject, ' in eye_components.mat']);
                    end
                else
                    disp('The file eye_components.mat is not found.');
                end
                
                % Define the name and path of the output file
                outputFileName = strrep(set_file, '_06_ICA', '_07_comp'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the name of the saved file to results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The data folder does not exist.');
end


%% Interpolation
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Define the main data folder path
data_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(data_folder)
    % Get the list of subject folders inside the data folder
    subjects = dir(data_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 1:length(subjects) % Change the starting index as necessary
        subject = subjects(i).name;
        eeg_folder = fullfile(data_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Get all the previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_07_comp.set'));
 
            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);
                
                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
                
                % Load the badChannelsStruct from badChannelsStruct.mat
                badChannelsPath = fullfile(main_path, 'badChannelsStruct.mat');
                
                if isfile(badChannelsPath)
                    % Load the structure from badChannelsStruct.mat
                    load(badChannelsPath, 'badChannelsStruct'); % Load the full struct
                    
                    % Find the subject's bad channels using subjectName
                    subjectIndex = find(strcmp({badChannelsStruct.subjectName}, subject));
                    
                    if ~isempty(subjectIndex)
                        % Get the list of bad channels for the current subject
                        badChannelsList = badChannelsStruct(subjectIndex).badChannels;
                        
                        % Perform the interpolation
                        EEG = pop_interp(EEG, badChannelsList, 'spherical');
                    else
                        disp(['No defective channels found for subject ', subject, ' in badChannelsStruct']);
                    end
                else
                    disp('The file badChannelsStruct.mat is not found.');
                end
                
                % Define the output file name and path
                outputFileName = strrep(set_file, '_07_comp.', '_08_inter'); 
                outputFilePath = fullfile(eeg_folder, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the saved file name to results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The data folder does not exist.');
end


%% Reference
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell to store results
results = {};

% Define the specific folder to search
data_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(data_folder)
    % Get the list of subject folders within the data folder
    subjects = dir(data_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 1:length(subjects)
        subject = subjects(i).name;
        eeg_folder = fullfile(data_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Read all the previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_08_inter.set'));
            
            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);
                
                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
                
                % Reference
                EEG = eeg_checkset(EEG);
                EEG = pop_reref(EEG, [], 'exclude', 32);
               
                % Define the name and path of the output file
                outputFileName = strrep(set_file, '_08_inter', '_09_ref'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the name of the saved file to results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The data folder does not exist.');
end
%%
HARD REMOVAL_NOT FOR THIS TIME 'CUZ ELIMINATES A LOT OF TRIALS WE NEED
% %% Remove Artifacts
% main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';
% 
% % Initialize a cell to store results
% results = {};
% 
% % Define the specific folder to search
% data_folder = fullfile(main_path, 'raw_data');
% 
% % Check if the folder exists
% if isfolder(data_folder)
%     % Get the list of subject folders within the data folder
%     subjects = dir(data_folder);
%     subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));
% 
%     for i = 1:length(subjects)
%         subject = subjects(i).name;
%         eeg_folder = fullfile(data_folder, subject, 'EEG');
% 
%         if isfolder(eeg_folder)
%             % Read all the previously generated .set files
%             set_files = dir(fullfile(eeg_folder, '*_09_ref.set'));
% 
%             for j = 1:length(set_files)
%                 set_file = set_files(j).name;
%                 setFilePath = fullfile(eeg_folder, set_file);
% 
%                 % Load the .set file
%                 EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
% 
%                 % Remove artifacts
%                 EEG = eeg_checkset(EEG);
%                 EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 'off', 'ChannelCriterion', 'off', 'LineNoiseCriterion', 'off', 'Highpass', 'off', 'BurstCriterion', 20, 'WindowCriterion', 0.25, 'BurstRejection', 'on', 'Distance', 'Euclidean', 'WindowCriterionTolerances', [-Inf 7]);
%                 [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'gui', 'off');
% 
%                 % Define the name and path of the output file
%                 outputFileName = strrep(set_file, '_09_ref', '_10_remArt'); % Change the name
%                 outputFilePath = fullfile(eeg_folder, outputFileName);
% 
%                 % Save the dataset in .set format
%                 EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
% 
%                 % Add the name of the saved file to results
%                 results{end+1} = outputFilePath; %#ok<AGROW>
%             end
%         end
%     end
% else
%     disp('The data folder does not exist.');
% end


%% Transformation
main_path = '/Users/jcedielescobar/Documents/Prepro/Day1';

% Initialize a cell array to store results
results = {};

% Define the specific folder to search
data_folder = fullfile(main_path, 'raw_data');

% Check if the folder exists
if isfolder(data_folder)
    % Get the list of subject folders within the data folder
    subjects = dir(data_folder);
    subjects = subjects([subjects.isdir] & ~ismember({subjects.name}, {'.', '..'}));

    for i = 1:length(subject) %Adjust the range if needed
        subject = subjects(i).name;
        eeg_folder = fullfile(data_folder, subject, 'EEG');
        
        if isfolder(eeg_folder)
            % Read all the previously generated .set files
            set_files = dir(fullfile(eeg_folder, '*_09_ref.set'));
            
            for j = 1:length(set_files)
                set_file = set_files(j).name;
                setFilePath = fullfile(eeg_folder, set_file);
                
                % Load the .set file
                EEG = pop_loadset('filename', set_file, 'filepath', eeg_folder);
                EEG = eeg_checkset(EEG);
             
                % Remove non-EEG channels for ERPlab plugin to work
                EEG = pop_select(EEG, 'rmchannel', {'ECG'});
                EEG = csdFromErplabAutomated(EEG); % Using defaults: m = 4 splines, lambda = 1e-5, 10 cm radius
                [~, EEG, ~] = pop_newset(ALLEEG, EEG, 1, 'gui', 'off');
    
                % Define the name and path of the output file
                outputFileName = strrep(set_file, '_09_ref', '_11_V2_CDS'); % Change the name
                outputFilePath = fullfile(eeg_folder, outputFileName);
                
                % Save the dataset in .set format
                EEG = pop_saveset(EEG, 'filename', outputFileName, 'filepath', eeg_folder);
                
                % Add the name of the saved file to results
                results{end+1} = outputFilePath; %#ok<AGROW>
            end
        end
    end
else
    disp('The data folder does not exist.');
end
