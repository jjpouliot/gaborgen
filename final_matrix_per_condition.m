clear;
%temp99 = eeglab;

files = dir('gaborgen*');  % List all folders starting with 'gaborgen'

% Filter out non-folder entries
dirFlags = [files.isdir];

% Extract the names of the folders
folderNames = {files(dirFlags).name};

% Remove '.' and '..' folders
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

% Display the folder names
disp('Folders in the current working directory:');
disp(folderNames);

% Initialize final matrix as a cell array
final_matrix = {}; 

% Loop over all subjects
for subindex = 1:length(folderNames)  % Loop over all subjects

    % Full path to the EEG directory of the current subject
    eegDir = fullfile(folderNames{subindex}, 'EEG');
    
    % Loop over the conditions (21 to 24)
    for condition = 21:24
        % Create the filename for the current condition
        basename = folderNames{subindex};
        currentFile = sprintf('%s_%d_EEG.mat.slidwin.mat', basename(end-2:end), condition);
        
        % Build the full path to the file
        currentFilePath = fullfile(eegDir, currentFile);
        
        % Check if the file exists
        if exist(currentFilePath, 'file')
            % Load the .mat file
            data = load(currentFilePath);
            
            % Create a struct to store the data with participant and condition info
            dataStruct = struct();
            dataStruct.participant = basename;    % Store participant's name
            dataStruct.condition = condition;     % Store condition number
            dataStruct.data = data;               % Store the slid wind output
            
            % Append the struct to final_matrix
            final_matrix{end+1} = dataStruct;  % Store the loaded data structure
            
        else
            warning(['File not found: ', currentFilePath]);
        end
    end
end

