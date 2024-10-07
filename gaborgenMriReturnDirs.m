function [participantDirectories, dataFolder, gaborgenCodeRepository] = gaborgenMriReturnDirs(partID, parentFolder, day1, day2)

dataFolder = [parentFolder '/raw_data'];
if ~(exist(dataFolder) == 7)
    error('raw_data folder not inside parentFolder')
end

if nargin < 3
    day1 = 1;
end

if nargin < 4
    day2 = 1;
end

if day1 == 0 && day2 == 0
    error('day1 and day2 arguments are zero which means neither day should be processed')
end

gaborgenCodeRepository = fileparts(mfilename('fullpath'));

allParticipantsDirectories = dir(dataFolder);
allParticipantsDirectories = allParticipantsDirectories(~ismember({allParticipantsDirectories.name}, {'.', '..'}));


%% set random number generator (just in case)
rng(1);

participantDirectories = {};
for i = 1:length(allParticipantsDirectories)
    dirname = allParticipantsDirectories(i).name;
    if contains(dirname, num2str(partID))
        participantDirectories{end+1} = dirname;
    end
end

%day 1 or 2 exclusion
finalDirs = {};
if day1 == 1
    for i = 1:length(participantDirectories)
        if ~contains(participantDirectories{i}, 'DAY2')
            finalDirs{end+1} = participantDirectories{i};
        end
    end

end

if day2 == 1
    for i = 1:length(participantDirectories)
        if contains(participantDirectories{i}, 'DAY2')
            finalDirs{end+1} = participantDirectories{i};
        end
    end
end

participantDirectories = finalDirs;

if length(participantDirectories) > 2
    error('There are more than two directories with the same subject number which should be impossible.')
end
